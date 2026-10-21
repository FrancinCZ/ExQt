# ROI drawing usability (feedback 2026-10-17): drawing nuclei in Napari was hard to follow.
# The ROI step now shows instructions with a nucleus counter, numbers each shape with the
# cell_id it will get, previews the final per-nucleus ROI (with a warning when an ROI covers
# slices where the mask has no objects, e.g. a stack that starts below the cells), and offers
# the ROI saved by a previous run so it does not have to be redrawn.
import os

import numpy as np
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from roi_tools import load_saved_roi, roi_coverage_warnings, roi_instructions, shape_nucleus_ids


def _two_nuclei_drawn():
    # Shape 1 and 3 belong to the right nucleus (drawn on z=1 and z=3), shape 2 to the left one.
    drawn = np.zeros((5, 30, 30), dtype=np.int32)
    drawn[1, 5:15, 18:28] = 1
    drawn[2, 5:15, 2:12] = 2
    drawn[3, 6:14, 19:27] = 3
    return drawn


def test_shape_ids_follow_final_cell_ids():
    from Batch import _interpolate_or_extrude_roi

    drawn = _two_nuclei_drawn()
    cells = _interpolate_or_extrude_roi(drawn)
    ids = shape_nucleus_ids(drawn, cells, n_shapes=3)
    assert ids[0] == ids[2] != ids[1]
    assert ids[0] == int(cells[1, 10, 20]) and ids[1] == int(cells[2, 10, 5])


def test_shape_ids_in_2d_are_the_shape_labels():
    drawn = np.zeros((20, 20), dtype=np.int32)
    drawn[2:6, 2:6] = 1
    drawn[10:15, 10:15] = 2
    assert shape_nucleus_ids(drawn, drawn, n_shapes=2) == [1, 2]


def test_shape_missing_from_labels_gets_none():
    drawn = np.zeros((20, 20), dtype=np.int32)
    drawn[2:6, 2:6] = 1
    assert shape_nucleus_ids(drawn, drawn, n_shapes=2) == [1, None]


def test_warns_when_roi_covers_slices_without_objects():
    cells = np.zeros((25, 20, 20), dtype=np.int32)
    cells[:, 2:10, 2:10] = 1
    cells[:, 12:18, 12:18] = 2
    mask = np.zeros_like(cells, dtype=np.uint8)
    mask[14:, 4:6, 4:6] = 1
    mask[:, 14:16, 14:16] = 1
    warnings = roi_coverage_warnings(cells, mask)
    assert len(warnings) == 1
    assert "Nucleus 1" in warnings[0] and "14 of 25" in warnings[0]


def test_no_warning_for_a_few_empty_edge_slices():
    cells = np.zeros((25, 20, 20), dtype=np.int32)
    cells[:, 2:10, 2:10] = 1
    mask = np.zeros_like(cells, dtype=np.uint8)
    mask[2:24, 4:6, 4:6] = 1
    assert roi_coverage_warnings(cells, mask) == []


def test_no_coverage_warning_in_2d():
    cells = np.ones((20, 20), dtype=np.int32)
    assert roi_coverage_warnings(cells, np.zeros((20, 20))) == []


def test_load_saved_roi(tmp_path):
    roi = np.zeros((5, 20, 20), dtype=np.uint8)
    roi[:, 2:8, 2:8] = 1
    roi[:, 10:15, 10:15] = 2
    path = tmp_path / "img_ROI.tif"
    tifffile.imwrite(path, roi)
    loaded = load_saved_roi(path, (5, 20, 20))
    assert loaded is not None and loaded.max() == 2
    assert load_saved_roi(path, (6, 20, 20)) is None
    assert load_saved_roi(tmp_path / "missing_ROI.tif", (5, 20, 20)) is None


def test_whole_field_roi_is_not_offered(tmp_path):
    # Auto-ROI and the empty fallback save an all-ones ROI: that is not a drawn ROI.
    path = tmp_path / "img_ROI.tif"
    tifffile.imwrite(path, np.ones((5, 20, 20), dtype=np.uint8))
    assert load_saved_roi(path, (5, 20, 20)) is None


def test_instructions_show_counter_and_keys():
    text = roi_instructions(is_3d=True, n_shapes=3, has_saved_roi=True)
    assert "3" in text and "E" in text and "P" in text and "Backspace" in text
    assert "saved" in text.lower()
    assert "Z" in text
    text_2d = roi_instructions(is_3d=False, n_shapes=0, has_saved_roi=False)
    assert "saved" not in text_2d.lower()


@pytest.fixture
def window(tmp_path):
    from unittest.mock import patch
    import napari
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QWidget

    app = QApplication.instance() or QApplication([])
    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        from App import ExQt
        win = ExQt()
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)
    win.viewer = napari.components.ViewerModel()
    yield win
    win.close()
    app.processEvents()


def _start_roi(window, shape=(5, 30, 30), saved=None):
    from unittest.mock import MagicMock

    window.viewer.add_image(np.zeros(shape), name="Signal")
    mask = np.zeros(shape, dtype=int)
    mask[:, 8:10, 22:24] = 1
    mask[:, 8:10, 5:7] = 1
    window.viewer.add_labels(mask, name="Mask Condensates")
    window.worker = MagicMock()
    window.prepare_manual_roi({"shape": shape, "is_3d": len(shape) == 3, "saved_roi_path": saved})
    return window.viewer.layers["Draw ROI"]


def _ellipse(z, y0, x0, y1, x1):
    return np.array([[z, y0, x0], [z, y0, x1], [z, y1, x1], [z, y1, x0]], dtype=float)


def test_roi_panel_counts_nuclei_and_preview_numbers_shapes(window):
    layer = _start_roi(window)
    assert not window.roi_help_label.isHidden()
    layer.add_ellipses([_ellipse(1, 5, 18, 15, 28), _ellipse(2, 5, 2, 15, 12), _ellipse(3, 6, 19, 14, 27)])
    assert "3" in window.roi_help_label.text()
    window.preview_roi()
    assert "ROI preview" in window.viewer.layers
    numbers = list(layer.features["nucleus"])
    assert numbers[0] == numbers[2] != numbers[1]
    assert "2 nuclei" in window.roi_help_label.text()


def test_confirm_hides_panel_and_removes_preview(window):
    layer = _start_roi(window)
    layer.add_ellipses([_ellipse(1, 5, 18, 15, 28)])
    window.preview_roi()
    window.confirm_roi()
    assert window.roi_help_label.isHidden()
    assert "ROI preview" not in window.viewer.layers
    assert window.worker.user_roi_source == "manual"
    window.worker.roi_event.set.assert_called_once()


def test_saved_roi_can_be_reused(window, tmp_path):
    roi = np.zeros((5, 30, 30), dtype=np.uint8)
    roi[:, 5:15, 18:28] = 1
    path = tmp_path / "img_ROI.tif"
    tifffile.imwrite(path, roi)
    _start_roi(window, saved=str(path))
    assert not window.btn_use_saved_roi.isHidden()
    window.use_saved_roi()
    assert window.worker.user_roi_source == "saved_roi"
    assert np.array_equal(window.worker.user_roi_data, roi)
    assert window.btn_use_saved_roi.isHidden()


def test_saved_roi_button_hidden_without_file(window):
    _start_roi(window)
    assert window.btn_use_saved_roi.isHidden()


def test_worker_offers_roi_saved_by_previous_run(tmp_path):
    from App import AnalysisWorker

    raw = np.full((5, 20, 20), 100, dtype=np.uint16)
    raw[1:3, 5:9, 5:9] = 600
    mask = np.zeros((5, 20, 20), dtype=np.uint8)
    mask[1:3, 5:9, 5:9] = 1
    tifffile.imwrite(tmp_path / "img1.tif", raw)
    tifffile.imwrite(tmp_path / "img1_Mask.tif", mask)
    (tmp_path / "out").mkdir()
    saved = np.zeros((5, 20, 20), dtype=np.uint8)
    saved[:, 2:12, 2:12] = 1
    tifffile.imwrite(tmp_path / "out" / "img1_ROI.tif", saved)

    worker = AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": False, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
    })
    infos = []

    def roi_slot(info):
        infos.append(info)
        worker.user_roi_data = load_saved_roi(info["saved_roi_path"], info["shape"])
        worker.user_roi_source = "saved_roi"
        worker.roi_event.set()

    worker.request_roi_signal.connect(roi_slot)
    worker.run()
    assert infos and infos[0]["saved_roi_path"].endswith("img1_ROI.tif")
    import pandas as pd
    df = pd.read_csv(next((tmp_path / "out").glob("*_Output_Batch_3d.csv")))
    assert set(df["roi_source"]) == {"saved_roi"}


# Follow-up 2026-10-17: in 3D every shape showed "?" until Preview ROI, which looked broken.
# Numbers are now computed live from the drawn slices only and match the final cell_ids.
def test_live_numbers_match_final_cell_ids_in_3d():
    from roi_tools import live_shape_numbers
    from Batch import _interpolate_or_extrude_roi
    import napari

    shapes = [_ellipse(3, 5, 18, 15, 28), _ellipse(1, 5, 2, 15, 12), _ellipse(1, 6, 19, 14, 27)]
    layer = napari.layers.Shapes(shapes, shape_type="ellipse")
    numbers, error = live_shape_numbers(layer.data, layer.shape_type, (5, 30, 30))
    assert error is None
    cells = _interpolate_or_extrude_roi(layer.to_labels(labels_shape=(5, 30, 30)))
    assert numbers == [str(n) for n in shape_nucleus_ids(layer.to_labels(labels_shape=(5, 30, 30)), cells, 3)]
    assert numbers[0] == numbers[2] != numbers[1]


def test_live_numbers_report_touching_shapes():
    from roi_tools import live_shape_numbers
    import napari

    layer = napari.layers.Shapes([_ellipse(1, 5, 2, 15, 12), _ellipse(1, 5, 11, 15, 21)], shape_type="ellipse")
    numbers, error = live_shape_numbers(layer.data, layer.shape_type, (5, 30, 30))
    assert error is not None and "touch" in error
    assert numbers == ["!", "!"]


def test_shapes_are_numbered_without_preview(window):
    layer = _start_roi(window)
    layer.add_ellipses([_ellipse(1, 5, 18, 15, 28), _ellipse(2, 5, 2, 15, 12), _ellipse(3, 6, 19, 14, 27)])
    numbers = list(layer.features["nucleus"])
    assert "?" not in numbers
    assert numbers[0] == numbers[2] != numbers[1]
    assert "2 nuclei" in window.roi_help_label.text()
