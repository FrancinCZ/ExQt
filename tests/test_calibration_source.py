# K2 (review 2026-09-23): the user's calibration choice must reach the worker,
# implausible TIFF calibration must not pass silently, and metadata.json must
# record the calibration actually applied to each file.
import json
import os

import numpy as np
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from Batch import get_metadata_from_tif
from calibration_policy import (
    looks_pre_expansion,
    resolve_file_calibration,
    validate_calibration,
)


def _write_imagej_pair(folder, name, pixel_um, z_um):
    raw = np.full((5, 20, 20), 100, dtype=np.uint16)
    raw[1:3, 5:9, 5:9] = 600
    mask = np.zeros((5, 20, 20), dtype=np.uint8)
    mask[1:3, 5:9, 5:9] = 1
    tifffile.imwrite(
        folder / f"{name}.tif", raw, imagej=True,
        resolution=(1.0 / pixel_um, 1.0 / pixel_um),
        metadata={"spacing": z_um, "unit": "um", "axes": "ZYX"},
    )
    tifffile.imwrite(folder / f"{name}_Mask.tif", mask)


# --- metadata reader ---------------------------------------------------------

def test_72_dpi_tiff_is_not_read_as_calibration(tmp_path):
    path = tmp_path / "plain.tif"
    tifffile.imwrite(
        path, np.zeros((3, 8, 8), dtype=np.uint16),
        resolution=(72.0, 72.0), resolutionunit="INCH", photometric="minisblack",
    )
    meta = get_metadata_from_tif(path)
    assert meta is None or "pixel_size" not in meta


def test_imagej_length_unit_calibration_still_detected(tmp_path):
    _write_imagej_pair(tmp_path, "fiji", 0.0145, 0.0625)
    meta = get_metadata_from_tif(tmp_path / "fiji.tif")
    assert meta["pixel_size"] == pytest.approx(14.5, rel=1e-4)
    assert meta["z_step"] == pytest.approx(62.5, rel=1e-4)


# --- pure policy -------------------------------------------------------------

DETECTED = {"a.tif": {"pixel_size_nm": 14.5, "z_step_nm": 62.5}}


def test_gui_source_ignores_detected_tiff_calibration():
    applied = resolve_file_calibration("a.tif", "gui", 58.0, 250.0, DETECTED)
    assert applied == {"pixel_size_nm": 58.0, "z_step_nm": 250.0, "source": "gui"}


def test_tiff_source_uses_detected_and_falls_back_to_gui():
    applied = resolve_file_calibration("a.tif", "tiff_per_file", 58.0, 250.0, DETECTED)
    assert applied == {"pixel_size_nm": 14.5, "z_step_nm": 62.5, "source": "tiff"}
    fallback = resolve_file_calibration("b.tif", "tiff_per_file", 58.0, 250.0, DETECTED)
    assert fallback == {
        "pixel_size_nm": 58.0, "z_step_nm": 250.0, "source": "gui_fallback_no_tiff_metadata",
    }


def test_unknown_source_is_rejected():
    with pytest.raises(ValueError):
        resolve_file_calibration("a.tif", "auto", 58.0, 250.0, DETECTED)


@pytest.mark.parametrize("xy, z", [(352778.0, 250.0), (0.5, 250.0), (58.0, 6000.0), (58.0, float("nan"))])
def test_implausible_calibration_is_an_error(xy, z):
    with pytest.raises(ValueError):
        validate_calibration(xy, z)
    with pytest.raises(ValueError):
        resolve_file_calibration("a.tif", "gui", xy, z, {})


def test_pre_expansion_calibration_is_flagged():
    assert looks_pre_expansion(14.5, 62.5, 58.0, 250.0, 4.0)
    assert looks_pre_expansion(15.5, 60.0, 58.0, 250.0, 4.0)  # within 10 %
    assert not looks_pre_expansion(58.0, 250.0, 58.0, 250.0, 4.0)
    assert not looks_pre_expansion(14.5, 62.5, 58.0, 250.0, 1.0)


def test_size_preview_follows_gui_choice(tmp_path):
    from size_preview import collect_size_preview

    _write_imagej_pair(tmp_path, "fiji", 0.0145, 0.0625)
    kwargs = dict(input_folder=tmp_path, mode="3d", expansion_factor=4.0,
                  pixel_size_nm=58.0, z_step_nm=250.0, min_voxels=1, signal_channel=0)
    gui = collect_size_preview(**kwargs, calibration_source="gui")
    tiff = collect_size_preview(**kwargs, calibration_source="tiff_per_file")
    # 32 voxels × (58/4)² × (250/4) nm³ vs the same with 14.5/62.5 → ExF³ = 64× smaller.
    assert gui["biological_size"].iloc[0] == pytest.approx(32 * 14.5**2 * 62.5 / 1e9)
    assert tiff["biological_size"].iloc[0] == pytest.approx(gui["biological_size"].iloc[0] / 64)


# --- GUI choice --------------------------------------------------------------

@pytest.fixture
def window(tmp_path):
    from unittest.mock import MagicMock, patch
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QWidget

    app = QApplication.instance() or QApplication([])
    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        from App import ExQt
        win = ExQt()
    # Never touch the user's real QSettings (registry) from tests.
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)
    win.settings.setValue("adv_pixel_size", 58.0)
    win.settings.setValue("adv_z_step", 250.0)
    win.exp_factor_spin.setValue(4.0)
    yield win
    win.close()
    app.processEvents()


@pytest.mark.parametrize("answer, expected", [("No", "gui"), ("Yes", "tiff_per_file")])
def test_mixed_calibration_choice_sets_source(window, tmp_path, answer, expected):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    data = tmp_path / "data"
    data.mkdir()
    _write_imagej_pair(data, "a", 0.058, 0.25)
    _write_imagej_pair(data, "b", 0.06, 0.3)
    with patch.object(QMessageBox, "question", return_value=getattr(QMessageBox, answer)):
        window._try_auto_fill_metadata(str(data))
    assert window.calibration_source == expected


def test_single_calibration_keeps_gui_source_and_flags_pre_expansion(window, tmp_path):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    data = tmp_path / "data"
    data.mkdir()
    _write_imagej_pair(data, "fiji", 0.0145, 0.0625)
    with patch.object(QMessageBox, "question", return_value=QMessageBox.No) as question:
        window._try_auto_fill_metadata(str(data))
    text = question.call_args.args[2]
    default_button = question.call_args.args[4]
    assert "PRE-EXPANSION" in text
    assert default_button == QMessageBox.No
    assert window.calibration_source == "gui"
    assert float(window.settings.value("adv_pixel_size")) == 58.0


def test_implausible_tiff_calibration_is_reported_and_ignored(window, tmp_path):
    from unittest.mock import patch
    from PySide6.QtWidgets import QMessageBox

    data = tmp_path / "data"
    data.mkdir()
    _write_imagej_pair(data, "huge", 5.0, 0.25)  # 5000 nm XY pixel
    with patch.object(QMessageBox, "warning") as warning, patch.object(QMessageBox, "question") as question:
        window._try_auto_fill_metadata(str(data))
    assert warning.called and "huge.tif" in warning.call_args.args[2]
    assert not question.called
    assert window.detected_metadata_by_file == {}


# --- worker end-to-end -------------------------------------------------------

def _run_worker(tmp_path, calibration_source, detected):
    from App import AnalysisWorker

    params = {
        "input_folder": str(tmp_path),
        "output_folder": str(tmp_path / "out"),
        "mode": "3d",
        "expansion_factor": 4.0,
        "min_voxels": 1,
        "auto_roi": True,
        "review_each_image": False,
        "show_napari": False,
        "generate_reports": False,
        "pixel_size_nm": 58.0,
        "z_step_nm": 250.0,
        "calibration_source": calibration_source,
        "detected_metadata_by_file": detected,
        "signal_channel": 0,
        "dapi_channel": 0,
    }
    worker = AnalysisWorker(params)
    messages = []
    worker.progress.connect(messages.append)
    worker.run()
    out = tmp_path / "out"
    csv = next(out.glob("*_Output_Batch_3d.csv"))
    meta = json.loads(next(out.glob("*_metadata.json")).read_text(encoding="utf-8"))
    import pandas as pd
    return pd.read_csv(csv), meta, messages


def test_worker_honours_gui_choice_and_records_applied_calibration(tmp_path):
    _write_imagej_pair(tmp_path, "fiji", 0.0145, 0.0625)
    detected = {"fiji.tif": {"pixel_size_nm": 14.5, "z_step_nm": 62.5}}
    df, meta, _ = _run_worker(tmp_path, "gui", detected)

    assert df["applied_pixel_size_nm"].unique().tolist() == [58.0]
    assert df["applied_z_step_nm"].unique().tolist() == [250.0]
    params = meta["parameters"]
    assert params["calibration_source"] == "gui"
    assert params["applied_calibration_by_file"]["fiji.tif"] == {
        "pixel_size_nm": 58.0, "z_step_nm": 250.0, "source": "gui",
    }
    assert params["pixel_size_nm"] == 58.0
    assert params["z_step_nm"] == 250.0


def test_worker_per_file_tiff_choice_records_tiff_source(tmp_path):
    _write_imagej_pair(tmp_path, "fiji", 0.0145, 0.0625)
    detected = {"fiji.tif": {"pixel_size_nm": 14.5, "z_step_nm": 62.5}}
    df, meta, _ = _run_worker(tmp_path, "tiff_per_file", detected)

    assert df["applied_pixel_size_nm"].iloc[0] == pytest.approx(14.5)
    applied = meta["parameters"]["applied_calibration_by_file"]["fiji.tif"]
    assert applied["source"] == "tiff"
    assert applied["pixel_size_nm"] == pytest.approx(14.5)


def test_worker_does_not_reread_tiff_behind_the_users_back(tmp_path):
    # Calibrated file, but no detected entry was passed: the worker must use the
    # GUI value (fallback), not silently read the TIFF itself.
    _write_imagej_pair(tmp_path, "fiji", 0.0145, 0.0625)
    df, meta, _ = _run_worker(tmp_path, "tiff_per_file", {})
    assert df["applied_pixel_size_nm"].iloc[0] == pytest.approx(58.0)
    assert (
        meta["parameters"]["applied_calibration_by_file"]["fiji.tif"]["source"]
        == "gui_fallback_no_tiff_metadata"
    )


def test_expansion_factor_accepts_two_decimals(window):
    """S7: ExF 4.25 must not be rounded to 4.2/4.3 (≈ 3.6 % volume error)."""
    window.exp_factor_spin.setValue(4.25)
    assert window.exp_factor_spin.value() == pytest.approx(4.25)
    window.exp_factor_spin.setValue(4.125)
    assert window.exp_factor_spin.value() == pytest.approx(4.125)
