# Detector offset from Leica .lif metadata: only the unambiguous photon-counting case sets
# a value; the source is recorded in the CSV and metadata.json.
import json
import os
import struct

import numpy as np
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from lif_metadata import detector_summary, offset_from_lif


def _detector(name, channel, mode, offset, gain="80", active="1"):
    return (f'<Detector Name="{name}" Type="SiPM" Channel="{channel}" IsActive="{active}" '
            f'Gain="{gain}" Offset="{offset}" AcquisitionMode="1" AcquisitionModeName="{mode}" '
            f'CanDoPhotonCounting="1"/>')


def _write_lif(path, detectors, images=2):
    """Minimal LIF: the same header layout as Leica files (0x70, size, 0x2A, UTF-16 XML)."""
    elements = "".join(f'<Element Name="img{i}"><Data><Image><Attachment>{"".join(detectors)}'
                       f"</Attachment></Image></Data></Element>" for i in range(images))
    xml = f"<LMSDataContainerHeader Version=\"2\">{elements}</LMSDataContainerHeader>"
    body = struct.pack("<bi", 0x2A, len(xml)) + xml.encode("utf-16-le")
    path.write_bytes(struct.pack("<ii", 0x70, len(body)) + body + b"\x00" * 64)
    return path


def test_photon_counting_detectors_give_zero_offset(tmp_path):
    lif = _write_lif(tmp_path / "a.lif", [
        _detector("HyD S 1", 1, "PhotonCounting", 0),
        _detector("HyD S 2", 2, "PhotonCounting", 0, gain="31.1"),
        _detector("HyD S 3", 3, "Standard", 7, active="0"),   # inactive: ignored
    ])
    summary = detector_summary(lif)
    assert summary["n_images"] == 2
    assert {d["name"] for d in summary["active_detectors"]} == {"HyD S 1", "HyD S 2"}
    value, source = offset_from_lif(lif)
    assert value == 0.0
    assert source["method"] == "lif_photon_counting"
    assert source["file"] == "a.lif" and len(source["sha256"]) == 64


def test_analog_mode_does_not_set_a_value(tmp_path):
    lif = _write_lif(tmp_path / "b.lif", [
        _detector("HyD S 1", 1, "PhotonCounting", 0),
        _detector("PMT 2", 2, "Standard", 2.5),
    ])
    value, source = offset_from_lif(lif)
    assert value is None
    assert "not photon counting" in source["message"]


def test_not_a_lif_file_is_rejected(tmp_path):
    bad = tmp_path / "x.lif"
    bad.write_bytes(b"not a lif file at all")
    with pytest.raises(ValueError, match="LIF"):
        detector_summary(bad)


def test_batch_records_offset_source(tmp_path):
    from Batch import process_condensates

    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    kwargs = dict(tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d",
                  min_voxels=1, auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0,
                  signal_channel=0, dapi_channel=0)
    assert set(process_condensates(**kwargs)["K_offset_source"]) == {"manual"}
    df = process_condensates(**kwargs, detector_offset_source="lif_photon_counting")
    assert set(df["K_offset_source"]) == {"lif_photon_counting"}
    assert set(df["K_offset_method"]) == {"explicit_setting"}


@pytest.fixture
def dialog(tmp_path):
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication
    from App import AdvancedSettingsDialog

    app = QApplication.instance() or QApplication([])
    dlg = AdvancedSettingsDialog()
    dlg.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)  # never the real registry
    dlg.detector_offset_spin.setValue(5.0)
    yield dlg
    dlg.close()
    app.processEvents()


def test_dialog_reads_lif_and_stores_source(dialog, tmp_path):
    from unittest.mock import patch
    from PySide6.QtWidgets import QFileDialog, QMessageBox

    lif = _write_lif(tmp_path / "exp.lif", [_detector("HyD S 1", 1, "PhotonCounting", 0)])
    with patch.object(QFileDialog, "getOpenFileName", return_value=(str(lif), "")), \
            patch.object(QMessageBox, "information") as info:
        dialog.read_offset_from_lif()
    assert info.called and "PhotonCounting" in info.call_args.args[2]
    assert dialog.detector_offset_spin.value() == 0.0
    dialog.accept()
    source = json.loads(dialog.settings.value("adv_detector_offset_source"))
    assert source["method"] == "lif_photon_counting" and source["file"] == "exp.lif"


def test_manual_change_after_lif_resets_source(dialog, tmp_path):
    from unittest.mock import patch
    from PySide6.QtWidgets import QFileDialog, QMessageBox

    lif = _write_lif(tmp_path / "exp.lif", [_detector("HyD S 1", 1, "PhotonCounting", 0)])
    with patch.object(QFileDialog, "getOpenFileName", return_value=(str(lif), "")), \
            patch.object(QMessageBox, "information"):
        dialog.read_offset_from_lif()
    dialog.detector_offset_spin.setValue(12.0)
    dialog.accept()
    assert json.loads(dialog.settings.value("adv_detector_offset_source"))["method"] == "manual"


def test_worker_writes_offset_source_to_metadata(tmp_path):
    from App import AnalysisWorker

    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    source = {"method": "lif_photon_counting", "file": "exp.lif", "sha256": "0" * 64}
    AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
        "detector_offset_adu": 0.0, "detector_offset_source": source,
    }).run()
    meta = json.loads(next((tmp_path / "out").glob("*_metadata.json")).read_text(encoding="utf-8"))
    assert meta["partitioning"]["offset_source"] == source
