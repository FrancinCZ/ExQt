# Cosmetic items of review 2026-09-23: single source of defaults, report files named
# per run, a clear error for 2D multi-channel input, and a napari preview that only
# shows measured objects at the right Z position.
import inspect

import numpy as np
import pandas as pd
import pytest
import tifffile

from Batch import process_condensates


def test_defaults_come_from_one_place():
    import Batch
    import postprocessing
    import size_preview
    from defaults import DEFAULT_SETTINGS

    batch = inspect.signature(Batch.process_condensates).parameters
    assert batch["min_voxels"].default == DEFAULT_SETTINGS["raw_min_voxels"]
    assert batch["mode_a_min_core_voxels"].default == DEFAULT_SETTINGS["mode_a_min_core_voxels"]
    assert inspect.signature(size_preview.collect_size_preview).parameters["min_voxels"].default == DEFAULT_SETTINGS["raw_min_voxels"]
    plots = inspect.signature(postprocessing.generate_plots).parameters
    assert plots["min_size"].default == DEFAULT_SETTINGS["plot_min_size"]
    assert plots["max_size"].default == DEFAULT_SETTINGS["plot_max_size"]
    assert Batch.MODE_A_Z_SPLIT_REVIEW_FRACTION == DEFAULT_SETTINGS["mode_a_z_split_review_fraction"]


def test_no_hard_coded_default_literals_left():
    from pathlib import Path

    project = Path(__file__).resolve().parents[1]
    app = (project / "App.py").read_text(encoding="utf-8")
    assert '"plot_min_size", 0.0001' not in app and '"plot_max_size", 2.0' not in app
    assert '"min_voxels", 5' not in app and '"mode_a_min_core_voxels", 20' not in app
    assert "350.0" not in (project / "size_preview.py").read_text(encoding="utf-8")


def _objects_frame(mode):
    return pd.DataFrame({
        "filename": ["a.tif"], "cell_id": [1], "object_id": [1], "mode": [mode],
        "volume_bio_um3" if mode == "3d" else "area_bio_um2": [0.5], "volume_px" if mode == "3d" else "area_px": [100],
    })


def test_report_files_of_3d_and_2d_runs_do_not_overwrite_each_other(tmp_path):
    from postprocessing import generate_excel_stats

    results = {}
    for mode in ("3d", "2d"):
        csv = tmp_path / f"exp_Output_Batch_{mode}.csv"
        _objects_frame(mode).to_csv(csv, index=False)
        results[mode] = generate_excel_stats(str(csv), min_size=0.0, max_size=10.0)
    for key in ("excel", "primary_csv", "excluded_csv"):
        assert results["3d"][key] != results["2d"][key]
        assert "_3d_" in str(results["3d"][key]) and "_2d_" in str(results["2d"][key])


def test_merge_ignores_per_run_report_files(tmp_path):
    from postprocessing import find_merge_source_csvs, generate_excel_stats

    csv = tmp_path / "exp_Output_Batch_3d.csv"
    _objects_frame("3d").to_csv(csv, index=False)
    generate_excel_stats(str(csv), min_size=0.0, max_size=10.0)
    assert find_merge_source_csvs(tmp_path) == [csv]


def test_2d_multichannel_tiff_has_a_clear_error(tmp_path):
    tifffile.imwrite(tmp_path / "img.tif", np.ones((2, 20, 20), dtype=np.uint16), ome=True, metadata={"axes": "CYX"})
    tifffile.imwrite(tmp_path / "img_Mask.tif", np.ones((20, 20), dtype=np.uint8))
    with pytest.raises(ValueError, match="2D multi-channel"):
        process_condensates(
            tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="2d",
            min_voxels=1, auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0,
            signal_channel=1, dapi_channel=0,
        )


def test_napari_points_only_for_measured_objects_with_z_scale(tmp_path):
    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:16, 10:16] = 1                        # measured object
    mask[3, 25, 25] = 1                                # 1 voxel, below min_voxels
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    layers = []
    df = process_condensates(
        tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d",
        min_voxels=5, auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0,
        signal_channel=0, dapi_channel=0, send_layer_func=layers.append,
    )
    points = [layer for layer in layers if layer.get("type") == "points"]
    assert len(df) == 1 and len(points) == 1
    assert len(points[0]["data"]) == 1
    assert points[0]["kwargs"]["scale"] == (3.0, 1, 1)
