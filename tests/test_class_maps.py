# Class images (feedback 2026-11-26): which object was counted as Globular / Elongated.
# Batch writes <image>_Objects.tif, reporting colours it into <image>_Classes.tif/.png.
# Also: per-layer voxel counts and the "Save next to input" output folder.
import os

import numpy as np
import pandas as pd
import tifffile

from Batch import object_map_path, process_condensates
from class_maps import CLASS_CODES, export_class_map, export_class_maps, object_class_codes
from partitioning_plots import _summary_table_rows
from postprocessing import REPORT_EXCLUDED_COLUMNS, REPORT_PRIMARY_COLUMNS, _prepare_reporting_frames

SHAPE = (12, 40, 40)
BIG = (slice(2, 10), slice(5, 15), slice(5, 15))
ROD = (slice(1, 11), slice(25, 28), slice(25, 28))
SMALL = (slice(5, 7), slice(30, 32), slice(5, 7))


def _run(tmp_path, **kwargs):
    raw = np.full(SHAPE, 100, dtype=np.uint16)
    mask = np.zeros(SHAPE, dtype=np.uint8)
    for box in (BIG, ROD, SMALL):
        mask[box] = 1
        raw[box] = 900
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    return process_condensates(
        tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d", min_voxels=1,
        auto_roi=True, pixel_size_nm=100.0, z_step_nm=100.0, signal_channel=0, dapi_channel=0,
        mode_a_enabled=True, mode_a_min_core_voxels=1, object_map_dir=tmp_path, **kwargs,
    )


def test_object_map_holds_measured_object_ids(tmp_path):
    df = _run(tmp_path)
    object_map = tifffile.imread(object_map_path("img.tif", tmp_path))
    assert object_map.shape == SHAPE
    assert set(np.unique(object_map)) - {0} == set(df["object_id"].astype(int))
    for row in df.itertuples():
        assert int(object_map[int(round(row.Z_px)), int(round(row.Y_px)), int(round(row.X_px))]) == row.object_id


def test_no_object_map_without_folder(tmp_path):
    raw = np.full(SHAPE, 100, dtype=np.uint16)
    mask = np.zeros(SHAPE, dtype=np.uint8)
    mask[BIG] = 1
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    process_condensates(tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d",
                        min_voxels=1, auto_roi=True, pixel_size_nm=100.0, z_step_nm=100.0,
                        signal_channel=0, dapi_channel=0)
    assert not object_map_path("img.tif", tmp_path).exists()


def test_layer_voxel_counts_add_up(tmp_path):
    df = _run(tmp_path)
    total = df["shell_voxels"] + df["middle_voxels"] + df["core_voxels"]
    assert (total == df["volume_px"]).all()
    assert (df["core_voxels"] == df["mode_a_core_voxels"]).all()
    assert "shell_voxels" in REPORT_PRIMARY_COLUMNS and "core_voxels" in REPORT_EXCLUDED_COLUMNS


def test_class_codes_follow_primary_and_class():
    report = pd.DataFrame({
        "object_id": [1, 2, 3, 4],
        "primary_qc_valid": [True, True, True, False],
        "condensate_class": ["Globular (Low FA) / Core-Enriched", "Elongated (High FA) / Profile-Unavailable",
                             "Punctate / Small Cluster", "Globular (Low FA) / Core-Enriched"],
    })
    assert object_class_codes(report) == {
        1: CLASS_CODES["Globular (Low FA)"], 2: CLASS_CODES["Elongated (High FA)"],
        3: CLASS_CODES["Primary, no shape class"], 4: CLASS_CODES["Not primary (QC / size range)"],
    }


def test_class_map_projection_prefers_primary(tmp_path):
    object_map = np.zeros((4, 6, 6), dtype=np.uint16)
    object_map[0, 1, 1] = 1  # excluded, below the elongated object in Z
    object_map[2, 1, 1] = 2  # elongated, primary
    report = pd.DataFrame({
        "filename": "a.tif", "object_id": [1, 2], "primary_qc_valid": [False, True],
        "condensate_class": ["Globular (Low FA) / Core-Enriched", "Elongated (High FA) / Core-Enriched"],
        "X_px": [1.0, 1.0], "Y_px": [1.0, 1.0], "Z_px": [0.0, 2.0],
    })
    tif_path, png_path = export_class_map(object_map, report, "a.tif", tmp_path, z_step_nm=300, pixel_size_nm=100)
    codes = tifffile.imread(tif_path)
    assert codes[0, 1, 1] == CLASS_CODES["Not primary (QC / size range)"]
    assert codes[2, 1, 1] == CLASS_CODES["Elongated (High FA)"]
    assert png_path.stat().st_size > 0
    from class_maps import _project
    assert _project(codes, 0)[1, 1] == CLASS_CODES["Elongated (High FA)"]


def test_run_writes_class_images(tmp_path):
    df = _run(tmp_path, keep_small_objects=True)
    df["filename"] = "img.tif"
    report = _prepare_reporting_frames(df, 0.0, 1e9)[3]
    written = export_class_maps(report, tmp_path)
    assert {p.name for p in written} == {"img_Classes.tif", "img_Classes.png"}
    codes = tifffile.imread(tmp_path / "img_Classes.tif")
    assert codes.shape == SHAPE and codes.max() > 0


def test_summary_counts_objects_without_profile():
    plot_df = pd.DataFrame({
        "volume_bio_um3": [0.1, 0.2, 0.3], "partition_coefficient": [2.0, 3.0, 4.0],
        "radial_profile_valid": ["True", "False", False],
        "condensate_class": ["Globular (Low FA) / Core-Enriched"] * 3,
    })
    rows = dict(_summary_table_rows(plot_df, "s"))
    assert rows["Too small for radial layers (no profile)"] == "2 (66.7%)"


def test_save_next_to_input_follows_input_folder(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    import App

    app = QApplication.instance() or QApplication([])
    window = App.ExQt.__new__(App.ExQt)
    from PySide6.QtWidgets import QCheckBox, QLineEdit
    window.folder_input = QLineEdit()
    window.output_path_edit = QLineEdit()
    window.same_folder_checkbox = QCheckBox()
    window.same_folder_checkbox.setChecked(True)
    window._sync_output_to_input("C:/data/run2")
    assert window.output_path_edit.text() == "C:/data/run2"
    window.same_folder_checkbox.setChecked(False)
    window._sync_output_to_input("C:/data/run3")
    assert window.output_path_edit.text() == "C:/data/run2"
