# K4 (review 2026-09-23): "Stop and Discard" while drawing a manual ROI must not
# save the current image, and a stale ROI from the previous image must never be reused.
import os

import numpy as np
import pytest
import pandas as pd
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _write_pair(folder, name):
    raw = np.full((5, 20, 20), 100, dtype=np.uint16)
    raw[1:3, 5:9, 5:9] = 600
    mask = np.zeros((5, 20, 20), dtype=np.uint8)
    mask[1:3, 5:9, 5:9] = 1
    tifffile.imwrite(folder / f"{name}.tif", raw)
    tifffile.imwrite(folder / f"{name}_Mask.tif", mask)


def _run(tmp_path, on_roi_request, review=False, on_review=None):
    from App import AnalysisWorker

    _write_pair(tmp_path, "img1")
    _write_pair(tmp_path, "img2")
    worker = AnalysisWorker({
        "input_folder": str(tmp_path),
        "output_folder": str(tmp_path / "out"),
        "mode": "3d",
        "expansion_factor": 1.0,
        "min_voxels": 1,
        "auto_roi": False,
        "review_each_image": review,
        "show_napari": False,
        "generate_reports": False,
        "pixel_size_nm": 100.0,
        "z_step_nm": 300.0,
        "calibration_source": "gui",
        "signal_channel": 0,
        "dapi_channel": 0,
    })
    messages = []
    calls = {"roi": 0, "review": 0}
    worker.progress.connect(messages.append)

    # Same-thread emit → direct call, so the "GUI" answers before the worker waits.
    def roi_slot(info):
        calls["roi"] += 1
        on_roi_request(worker, info, calls["roi"])

    def review_slot():
        calls["review"] += 1
        (on_review or (lambda w, n: w.review_event.set()))(worker, calls["review"])

    worker.request_roi_signal.connect(roi_slot)
    worker.request_review_signal.connect(review_slot)
    worker.run()
    csvs = list((tmp_path / "out").glob("*_Output_Batch_3d.csv"))
    df = pd.read_csv(csvs[0]) if csvs else pd.DataFrame(columns=["filename"])
    return df, messages, calls


def _full_roi(worker, info):
    worker.user_roi_data = np.ones(info["shape"], dtype=int)
    worker.roi_event.set()


def test_abort_during_second_roi_keeps_only_first_image(tmp_path):
    def gui(worker, info, n):
        if n == 1:
            _full_roi(worker, info)
        else:
            worker.request_abort()

    df, messages, calls = _run(tmp_path, gui)
    assert calls["roi"] == 2
    assert set(df["filename"]) == {"img1.tif"}
    assert "Discarded and stopped at: img2.tif" in messages
    assert not (tmp_path / "img2_ROI.tif").exists()
    assert not (tmp_path / "out" / "img2_ROI.tif").exists()      # ROI now lives in the run folder
    assert (tmp_path / "out" / "img1_ROI.tif").exists()


def test_abort_during_second_roi_with_review_enabled(tmp_path):
    def gui(worker, info, n):
        if n == 1:
            _full_roi(worker, info)
        else:
            worker.request_abort()

    df, messages, calls = _run(tmp_path, gui, review=True)
    assert set(df["filename"]) == {"img1.tif"}
    assert calls["review"] == 1  # img2 never reaches review
    assert "Discarded and stopped at: img2.tif" in messages


def test_stale_roi_from_previous_image_is_never_reused(tmp_path):
    # The GUI wakes the worker for img2 without providing a new ROI.
    def gui(worker, info, n):
        if n == 1:
            _full_roi(worker, info)
        else:
            worker.roi_event.set()

    df, messages, _ = _run(tmp_path, gui)
    assert set(df["filename"]) == {"img1.tif"}
    assert any(m.startswith("Error") and "img2.tif" in m for m in messages)


def test_abort_during_review_still_discards_current_image(tmp_path):
    def review(worker, n):
        if n == 1:
            worker.review_event.set()
        else:
            worker.request_abort()

    df, messages, _ = _run(tmp_path, _full_roi_first_arg, review=True, on_review=review)
    assert set(df["filename"]) == {"img1.tif"}
    assert "Discarded and stopped at: img2.tif" in messages


def _full_roi_first_arg(worker, info, n):
    _full_roi(worker, info)


# --- K4.4: empty ROI → explicit dialog + roi_source column ---------------------

def test_roi_source_column_manual_and_auto(tmp_path):
    df, _, _ = _run(tmp_path, _full_roi_first_arg)
    assert set(df["roi_source"]) == {"manual"}

    from Batch import process_condensates
    auto = process_condensates(
        tif_path=tmp_path / "img1.tif", mask_path=tmp_path / "img1_Mask.tif",
        mode="3d", min_voxels=1, auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0,
        signal_channel=0, dapi_channel=0,
    )
    assert set(auto["roi_source"]) == {"auto_fov"}


def test_empty_roi_fallback_is_recorded(tmp_path):
    def gui(worker, info, n):
        _full_roi(worker, info)
        if n == 2:
            worker.user_roi_source = "empty_fallback_fov"

    df, _, _ = _run(tmp_path, gui)
    sources = df.groupby("filename")["roi_source"].unique().to_dict()
    assert list(sources["img1.tif"]) == ["manual"]
    assert list(sources["img2.tif"]) == ["empty_fallback_fov"]


class _FakeWorker:
    def __init__(self):
        import threading
        self.roi_event = threading.Event()
        self.user_roi_data = None
        self.user_roi_source = "manual"


@pytest.fixture
def window(tmp_path):
    from unittest.mock import patch
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QWidget

    app = QApplication.instance() or QApplication([])
    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        viewer.return_value.layers.__iter__.return_value = iter([])
        viewer.return_value.layers.__contains__.return_value = False
        from App import ExQt
        win = ExQt()
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)
    win.worker = _FakeWorker()
    win._current_roi_info = {"shape": (5, 20, 20), "is_3d": True}
    yield win
    win.close()
    app.processEvents()


def test_empty_roi_asks_and_can_return_to_drawing(window):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    with patch.object(QMessageBox, "question", return_value=QMessageBox.No) as question:
        window.confirm_roi()
    assert question.called
    assert not window.worker.roi_event.is_set()
    assert window.worker.user_roi_data is None


def test_empty_roi_confirmed_uses_fov_and_marks_source(window):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    with patch.object(QMessageBox, "question", return_value=QMessageBox.Yes):
        window.confirm_roi()
    assert window.worker.roi_event.is_set()
    assert window.worker.user_roi_data.shape == (5, 20, 20)
    assert window.worker.user_roi_data.all()
    assert window.worker.user_roi_source == "empty_fallback_fov"


# --- K3: ambiguous multi-nucleus ROI is rejected in the GUI, shapes kept -------

def _set_drawn_labels(window, labels):
    from unittest.mock import MagicMock

    layer = MagicMock()
    layer.to_labels.return_value = labels
    window.viewer.layers.__contains__.side_effect = lambda name: name == "Draw ROI"
    window.viewer.layers.__getitem__.return_value = layer
    window.viewer.layers.remove.reset_mock()


def test_ambiguous_roi_is_rejected_and_drawing_kept(window):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    labels = np.zeros((5, 20, 20), dtype=np.int64)
    labels[2, 2:8, 2:8] = 1
    labels[2, 8:14, 2:8] = 2  # touching shape on the same slice
    _set_drawn_labels(window, labels)
    with patch.object(QMessageBox, "warning") as warning:
        window.confirm_roi()
    assert warning.called and "same Z-slice" in warning.call_args.args[2]
    assert not window.worker.roi_event.is_set()
    window.viewer.layers.remove.assert_not_called()


def test_valid_multi_nucleus_roi_is_sent(window):
    labels = np.zeros((5, 20, 20), dtype=np.int64)
    labels[1, 2:6, 2:6] = 1
    labels[3, 12:16, 12:16] = 2
    _set_drawn_labels(window, labels)
    window.confirm_roi()
    assert window.worker.roi_event.is_set()
    np.testing.assert_array_equal(window.worker.user_roi_data, labels)
    assert window.worker.user_roi_source == "manual"
    window.viewer.layers.remove.assert_called_once_with("Draw ROI")
