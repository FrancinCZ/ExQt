# S1 (review 2026-09-23): Batch must follow the aligner's OVERALL status, REVIEW
# stacks are analysed but kept out of primary, and excluded files are recorded.
import json
import os

import numpy as np
import pandas as pd
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from Batch import process_condensates, validate_alignment_qc


def _pair(folder, name="img"):
    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    raw[2:4, 10:14, 10:14] = 700
    tifffile.imwrite(folder / f"{name}.tif", raw)
    tifffile.imwrite(folder / f"{name}_Mask.tif", mask)


def _drift(folder, name, steps, overall=None):
    rows = [{"status": s, "reason": "insufficient_texture" if s == "FAIL" else "",
             "step_shift_y_px": 0.0, "step_shift_x_px": 0.0, "step_magnitude_px": 0.0} for s in steps]
    table = pd.DataFrame(rows)
    if overall is not None:
        table["overall_status"] = overall
    table.to_csv(folder / f"{name}_drift.csv", index=False)


def _process(folder, name="img"):
    return process_condensates(
        tif_path=folder / f"{name}.tif", mask_path=folder / f"{name}_Mask.tif", mode="3d",
        min_voxels=1, auto_roi=True, pixel_size_nm=100.0, z_step_nm=300.0,
        signal_channel=0, dapi_channel=0,
    )


def test_aligner_writes_overall_status(tmp_path):
    from scipy import ndimage as ndi
    from stack_aligner import align_tiff_stack

    y, x = np.indices((80, 92), dtype=float)
    base = (3200 * np.exp(-((y - 20) ** 2 + (x - 28) ** 2) / 85)
            + 2700 * np.exp(-((y - 56) ** 2 + (x - 64) ** 2) / 115))
    planes = [np.zeros_like(base)]  # empty first slice → insufficient_texture at the edge
    planes += [ndi.shift(base + 400, (i * 0.5, -i * 0.5), order=1, prefilter=False) for i in range(4)]
    planes += [np.zeros_like(base)]  # empty last slice
    source = tmp_path / "in" / "s.tif"
    source.parent.mkdir()
    tifffile.imwrite(source, np.stack(planes).astype(np.uint16), ome=True, metadata={"axes": "ZYX"})

    result = align_tiff_stack(source, tmp_path / "out")
    table = pd.read_csv(result.drift_csv)
    assert "FAIL" in set(table["status"])               # edge steps
    assert set(table["overall_status"]) == {result.alignment.status} == {"PASS"}


def test_overall_pass_with_edge_fail_steps_is_analysed(tmp_path):
    _pair(tmp_path)
    _drift(tmp_path, "img", ["FAIL", "PASS", "PASS", "PASS", "FAIL"], overall="PASS")
    df = _process(tmp_path)
    assert not df.empty
    assert set(df["alignment_status"]) == {"PASS"}


def test_review_is_analysed_but_not_primary(tmp_path):
    from postprocessing import _prepare_reporting_frames

    _pair(tmp_path)
    _drift(tmp_path, "img", ["PASS", "FAIL", "PASS", "PASS", "PASS"], overall="REVIEW")
    df = _process(tmp_path)
    assert set(df["alignment_status"]) == {"REVIEW"}
    _, primary, _, report_df, _ = _prepare_reporting_frames(df, 0.0, 1e6)
    assert primary.empty
    assert not report_df["primary_qc_valid"].any()


def test_not_aligned_input_is_marked(tmp_path):
    _pair(tmp_path)
    assert validate_alignment_qc(tmp_path / "img.tif") is None
    assert set(_process(tmp_path)["alignment_status"]) == {"not_aligned"}


def test_overall_fail_is_rejected(tmp_path):
    _pair(tmp_path)
    _drift(tmp_path, "img", ["PASS", "FAIL", "FAIL", "FAIL", "PASS"], overall="FAIL")
    with pytest.raises(ValueError, match="FAIL"):
        _process(tmp_path)


def test_legacy_drift_csv_without_overall_status(tmp_path):
    _pair(tmp_path)
    _drift(tmp_path, "img", ["PASS", "FAIL", "PASS", "PASS", "PASS"])
    with pytest.raises(ValueError, match="re-align"):
        validate_alignment_qc(tmp_path / "img.tif")
    _drift(tmp_path, "img", ["PASS", "REVIEW", "PASS", "PASS", "PASS"])
    assert validate_alignment_qc(tmp_path / "img.tif") == "REVIEW"


def test_worker_records_excluded_and_review_files(tmp_path):
    from App import AnalysisWorker

    for name in ("a_pass", "b_review", "c_fail"):
        _pair(tmp_path, name)
    _drift(tmp_path, "a_pass", ["PASS"] * 5, overall="PASS")
    _drift(tmp_path, "b_review", ["PASS", "REVIEW", "PASS", "PASS", "PASS"], overall="REVIEW")
    _drift(tmp_path, "c_fail", ["FAIL"] * 5, overall="FAIL")
    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    tifffile.imwrite(tmp_path / "d_nomask.tif", raw)

    worker = AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
    })
    worker.run()
    out = tmp_path / "out"
    df = pd.read_csv(next(out.glob("*_Output_Batch_3d.csv")))
    meta = json.loads(next(out.glob("*_metadata.json")).read_text(encoding="utf-8"))

    assert set(df["filename"]) == {"a_pass.tif", "b_review.tif"}
    excluded = {item["name"]: item["reason"] for item in meta["files"]["excluded_files"]}
    assert set(excluded) == {"c_fail.tif", "d_nomask.tif"}
    assert "FAIL" in excluded["c_fail.tif"]
    assert "mask" in excluded["d_nomask.tif"].lower()
    assert meta["files"]["alignment_review_files"] == ["b_review.tif"]
    assert meta["files"]["analysed_files"] == ["a_pass.tif", "b_review.tif"]
