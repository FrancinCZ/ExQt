# S4 (review 2026-09-23, VALIDATION_PROTOCOL §2): every run must record which code,
# which input files and which settings produced it; runs with a different meaning
# of K or of the classes must not share a QC fingerprint.
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

PROJECT = Path(__file__).resolve().parents[1]


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _pair(folder, name="img"):
    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    raw[2:4, 10:14, 10:14] = 700
    tifffile.imwrite(folder / f"{name}.tif", raw)
    tifffile.imwrite(folder / f"{name}_Mask.tif", mask)


def test_code_provenance_records_git_hashes_and_versions():
    from provenance import code_provenance

    info = code_provenance(PROJECT)
    assert info["source_sha256"]["Batch.py"] == _sha256(PROJECT / "Batch.py")
    assert info["source_sha256"]["rezim_a_metrics.py"] == _sha256(PROJECT / "rezim_a_metrics.py")
    assert info["packages"]["numpy"] == np.__version__
    assert info["python"]
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT, capture_output=True, text=True)
    if head.returncode == 0:
        assert info["git_commit"] == head.stdout.strip()
        assert isinstance(info["git_dirty"], bool)


def test_code_provenance_without_git(monkeypatch):
    import provenance

    def no_git(*args, **kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenance.subprocess, "run", no_git)
    info = provenance.code_provenance(PROJECT)
    assert info["git_commit"] is None and info["git_dirty"] is None
    assert info["source_sha256"]["Batch.py"] == _sha256(PROJECT / "Batch.py")


def test_input_record_hashes_raw_mask_and_alignment_sidecars(tmp_path):
    from provenance import input_file_record

    _pair(tmp_path)
    (tmp_path / "img_drift.csv").write_text("status\nPASS\n", encoding="utf-8")
    record = input_file_record(tmp_path / "img.tif", tmp_path / "img_Mask.tif")
    assert record["name"] == "img.tif"
    assert record["sha256"] == _sha256(tmp_path / "img.tif")
    assert record["size_bytes"] == (tmp_path / "img.tif").stat().st_size
    assert record["mask_sha256"] == _sha256(tmp_path / "img_Mask.tif")
    assert record["sidecars"]["img_drift.csv"] == _sha256(tmp_path / "img_drift.csv")


def test_worker_metadata_contains_provenance_inputs_and_auto_roi(tmp_path):
    from App import AnalysisWorker

    _pair(tmp_path)
    worker = AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
    })
    worker.run()
    meta = json.loads(next((tmp_path / "out").glob("*_metadata.json")).read_text(encoding="utf-8"))
    assert meta["parameters"]["auto_roi"] is True
    assert meta["provenance"]["source_sha256"]["Batch.py"] == _sha256(PROJECT / "Batch.py")
    inputs = meta["files"]["input_files"]
    assert [item["name"] for item in inputs] == ["img.tif"]
    assert inputs[0]["sha256"] == _sha256(tmp_path / "img.tif")


def _meta(auto_roi=False, offset_method="explicit_setting", offset=0.0, alpha=0.01):
    return {
        "parameters": {"mode": "3d", "expansion_factor": 4.0, "min_voxels": 5, "auto_roi": auto_roi},
        "channels": {"signal_channel": 1, "dapi_channel": 0},
        "partitioning": {"K_offset_method": offset_method, "detector_offset_adu": offset},
        "classification": {"gradient_test": "welch_t_two_sided", "gradient_alpha": alpha},
        "provenance": {"git_commit": "abc"},
    }


def test_fingerprint_separates_runs_with_different_k_or_class_meaning():
    from postprocessing import qc_policy_fingerprint

    df = pd.DataFrame({"volume_bio_um3": [0.1], "mode": ["3d"]})
    base, entries = qc_policy_fingerprint(_meta(), df, 0.0, 1.0)
    assert qc_policy_fingerprint(_meta(auto_roi=True), df, 0.0, 1.0)[0] != base
    assert qc_policy_fingerprint(_meta(offset_method="percentile_legacy"), df, 0.0, 1.0)[0] != base
    assert qc_policy_fingerprint(_meta(alpha=0.05), df, 0.0, 1.0)[0] != base
    # Offset value and code revision are shown but may legitimately differ between runs.
    assert qc_policy_fingerprint(_meta(offset=100.0), df, 0.0, 1.0)[0] == base
    shown = {(section, name) for section, name, *_ in entries}
    assert {("partitioning", "detector_offset_adu"), ("provenance", "git_commit")} <= shown


def test_fingerprint_falls_back_to_columns_for_old_metadata():
    from postprocessing import qc_policy_fingerprint

    manual = pd.DataFrame({"volume_bio_um3": [0.1], "roi_source": ["manual"], "K_offset_method": ["explicit_setting"]})
    auto = manual.assign(roi_source="auto_fov")
    assert qc_policy_fingerprint({}, manual)[0] != qc_policy_fingerprint({}, auto)[0]
