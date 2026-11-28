# Small objects (POL II, 2026-11-21): objects too small for three radial layers were excluded and
# biased the primary objects towards large ones; FA_xy and the main-axis angle separate a Z elongation
# from shape; saturated (clipped) intensities are flagged because they make K_part unreliable.
import numpy as np
import pandas as pd
import tifffile

from Batch import process_condensates, saturation_level
from partitioning_plots import _summary_table_rows
from postprocessing import _prepare_reporting_frames, detect_red_flags, qc_policy_fingerprint
from rezim_a_anisotropy import shape_anisotropy

SHAPE = (12, 48, 48)
LARGE = (slice(2, 10), slice(6, 22), slice(6, 22))
SMALL = (slice(5, 7), slice(32, 35), slice(32, 35))


def test_fa_xy_and_angle_separate_z_elongation():
    z_rod = np.zeros((30, 10, 10), bool)
    z_rod[5:25, 4:6, 4:6] = True
    result = shape_anisotropy(z_rod)
    assert result["fractional_anisotropy"] > 0.9
    assert result["fa_xy"] < 0.05
    assert result["main_axis_angle_to_z_deg"] < 1
    xy_rod = np.zeros((6, 30, 10), bool)
    xy_rod[2:4, 5:25, 4:6] = True
    result = shape_anisotropy(xy_rod)
    assert result["fa_xy"] > 0.9
    assert abs(result["main_axis_angle_to_z_deg"] - 90) < 1


def test_saturation_level_uses_reached_bit_depth():
    image = np.zeros((4, 4), np.uint16)
    assert saturation_level(image) == 65535
    image[0, 0] = 4095
    assert saturation_level(image) == 4095
    image[0, 0] = 65535
    assert saturation_level(image) == 65535
    assert saturation_level(image.astype(np.float32)) is None


def _run(tmp_path, saturate=False, **kwargs):
    raw = np.full(SHAPE, 100, dtype=np.uint16)
    mask = np.zeros(SHAPE, dtype=np.uint8)
    for box in (LARGE, SMALL):
        mask[box] = 1
        raw[box] = 800
    raw[LARGE] = np.random.default_rng(0).integers(700, 900, raw[LARGE].shape)
    if saturate:
        raw[SMALL] = 65535
        raw[2, 6:22, 6:22] = 65535  # one plane of the large object: ~14 % of all object voxels
    roi = np.zeros(SHAPE, dtype=np.uint8)
    roi[:, 1:47, 1:47] = 1
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    return process_condensates(
        tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d", min_voxels=1,
        auto_roi=False, request_roi_func=lambda shape, is_3d: roi.copy(), roi_output_dir=tmp_path,
        pixel_size_nm=100.0, z_step_nm=100.0, signal_channel=0, dapi_channel=0,
        mode_a_enabled=True, mode_a_min_core_voxels=8, allow_z_edge=True, **kwargs,
    )


def _small(df):
    return df.loc[df["volume_px"] == df["volume_px"].min()].iloc[0]


def test_small_object_excluded_by_default(tmp_path):
    df = _run(tmp_path)
    small = _small(df)
    assert not bool(small["radial_profile_valid"])
    assert not bool(small["mode_a_primary_include"])
    _, primary, _, _, _ = _prepare_reporting_frames(df, 0.0, 1e9)
    assert len(primary) == 1


def test_keep_small_objects_keeps_size_k_and_object_fa(tmp_path):
    df = _run(tmp_path, keep_small_objects=True)
    small = _small(df)
    assert bool(small["small_objects_kept"]) and not bool(small["radial_profile_valid"])
    assert small["radial_profile_reason"]
    assert bool(small["mode_a_primary_include"]) and small["mode_a_qc_reason"] == ""
    assert np.isnan(small["A_shell"]) and np.isnan(small["Delta_A_core_shell"])
    assert np.isfinite(small["A_object"]) and np.isfinite(small["partition_coefficient"])
    assert "Unclassified (QC)" not in small["condensate_class"]
    _, primary, _, report, _ = _prepare_reporting_frames(df, 0.0, 1e9)
    assert len(primary) == 2
    # fa_complete still means a complete radial profile.
    assert report["fa_complete"].sum() == 1
    assert qc_policy_fingerprint({}, df, 0, 1)[0] != qc_policy_fingerprint({}, _run(tmp_path), 0, 1)[0]


def test_saturation_is_counted_and_flagged(tmp_path):
    df = _run(tmp_path, saturate=True)
    small = _small(df)
    assert small["saturated_voxels"] == small["volume_px"]
    assert df["saturated_voxels"].max() == 256
    assert set(df["saturation_level_adu"]) == {65535.0}
    flags = [flag for flag in detect_red_flags(df) if flag["flag"] == "saturated_signal"]
    assert len(flags) == 1 and "raw image" in flags[0]["message"]
    assert not [flag for flag in detect_red_flags(_run(tmp_path)) if flag["flag"] == "saturated_signal"]


def test_pdf_table_reports_fa_xy_angle_and_saturation():
    frame = pd.DataFrame({
        "volume_bio_um3": [0.1, 0.2], "volume_px": [100, 100], "saturated_voxels": [50, 0],
        "FA_xy": [0.2, 0.4], "main_axis_angle_to_z_deg": [10.0, 80.0], "partition_coefficient": [2.0, 3.0],
    })
    rows = dict(_summary_table_rows(frame, "s"))
    assert rows["Median FA in XY only (FA_xy)"] == "0.30"
    assert rows["Main axis within 30° of Z"].startswith("50.0%")
    assert rows["Saturated object voxels"].startswith("25.0%") and "unreliable" in rows["Saturated object voxels"]


def test_settings_dialog_stores_small_object_policy(tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication
    from App import AdvancedSettingsDialog

    app = QApplication.instance() or QApplication([])
    settings = QSettings(str(tmp_path / "s.ini"), QSettings.IniFormat)
    dialog = AdvancedSettingsDialog()
    dialog.settings = settings
    dialog.load_adv_settings()
    assert not dialog.keep_small_check.isChecked()
    dialog.keep_small_check.setChecked(True)
    dialog.accept()
    assert str(settings.value("qc_keep_small_objects")).lower() == "true"
