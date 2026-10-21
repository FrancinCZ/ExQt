# IF stacks that do not cover the whole nucleus in Z (feedback 2026-10-17): almost every object
# touches the first/last slice and was excluded as truncated. An explicit, recorded option keeps
# objects cut only by the Z end of the stack; objects at the XY border, alignment padding or the ROI
# boundary stay excluded. Default unchanged (exclude).
import numpy as np
import tifffile

from Batch import process_condensates
from postprocessing import _prepare_reporting_frames, qc_policy_fingerprint

SHAPE = (8, 40, 40)
INTERIOR = (slice(3, 5), slice(18, 22), slice(18, 22))
XY_EDGE = (slice(3, 5), slice(0, 4), slice(18, 22))
Z_EDGE = (slice(0, 2), slice(30, 34), slice(5, 9))
Z_TOP = (slice(6, 8), slice(5, 9), slice(30, 34))


def _run(tmp_path, **kwargs):
    raw = np.full(SHAPE, 100, dtype=np.uint16)
    mask = np.zeros(SHAPE, dtype=np.uint8)
    for box in (INTERIOR, XY_EDGE, Z_EDGE, Z_TOP):
        mask[box] = 1
        raw[box] = 800
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    return process_condensates(
        tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d", min_voxels=1,
        auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0, signal_channel=0, dapi_channel=0, **kwargs,
    )


def _primary_ids(df):
    _, primary, _, report, _ = _prepare_reporting_frames(df, 0.0, 1e9)
    return set(report.loc[report["primary_qc_valid"], "Z_px"].apply(np.floor).astype(int))


def test_z_edge_columns_are_recorded(tmp_path):
    df = _run(tmp_path)
    by_z = df.set_index(df["Z_px"].apply(np.floor).astype(int))
    assert bool(by_z.loc[0, "touches_z_edge_only"]) and bool(by_z.loc[6, "touches_z_edge_only"])
    xy = df[df["Y_px"] < 3].iloc[0]
    assert bool(xy["touches_image_edge"]) and not bool(xy["touches_z_edge_only"])
    assert not df["z_edge_allowed"].any()


def test_default_excludes_z_edge_objects(tmp_path):
    df = _run(tmp_path)
    assert len(df[df["touches_image_edge"].astype(bool)]) == 3
    assert _primary_ids(df) == {3}


def test_allow_z_edge_keeps_z_cut_objects_but_not_xy_edge(tmp_path):
    df = _run(tmp_path, allow_z_edge=True)
    assert df["z_edge_allowed"].all()
    assert _primary_ids(df) == {0, 3, 6}
    # The flag stays truthful: the object still touches the stack edge.
    assert df["touches_image_edge"].astype(bool).sum() == 3


def test_allow_z_edge_with_mode_a(tmp_path):
    df = _run(tmp_path, allow_z_edge=True, mode_a_enabled=True, mode_a_min_core_voxels=1)
    reasons = df.set_index(df["Z_px"].apply(np.floor).astype(int))["mode_a_qc_reason"].fillna("")
    assert "touches_image_edge" not in reasons.loc[0]
    xy_reason = df.loc[df["Y_px"] < 3, "mode_a_qc_reason"].iloc[0]
    assert "touches_image_edge" in xy_reason


def test_policy_changes_fingerprint(tmp_path):
    df = _run(tmp_path)
    allowed = _run(tmp_path, allow_z_edge=True)
    assert qc_policy_fingerprint({}, df, 0, 1)[0] != qc_policy_fingerprint({}, allowed, 0, 1)[0]


def test_settings_dialog_stores_z_edge_policy(tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication
    from App import AdvancedSettingsDialog

    app = QApplication.instance() or QApplication([])
    settings = QSettings(str(tmp_path / "s.ini"), QSettings.IniFormat)
    dialog = AdvancedSettingsDialog()
    dialog.settings = settings  # never the real registry
    dialog.load_adv_settings()
    assert not dialog.allow_z_edge_check.isChecked()
    dialog.allow_z_edge_check.setChecked(True)
    dialog.accept()
    assert str(settings.value("qc_allow_z_edge")).lower() == "true"
