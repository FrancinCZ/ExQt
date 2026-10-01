# The ROI of a run is saved with that run's outputs, so a new run never overwrites the ROI
# of an earlier one and the input folder stays untouched.
import json
import os

import numpy as np
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _pair(folder, name="img"):
    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    raw[2:4, 10:14, 10:14] = 700
    tifffile.imwrite(folder / f"{name}.tif", raw)
    tifffile.imwrite(folder / f"{name}_Mask.tif", mask)


def _roi(shape, is_3d):
    roi = np.zeros(shape, dtype=np.int64)
    roi[3, 5:25, 5:25] = 1
    return roi


def test_process_condensates_saves_roi_in_given_folder(tmp_path):
    from Batch import process_condensates

    _pair(tmp_path)
    out = tmp_path / "run1"
    process_condensates(
        tif_path=tmp_path / "img.tif", mask_path=tmp_path / "img_Mask.tif", mode="3d",
        min_voxels=1, auto_roi=False, request_roi_func=_roi, pixel_size_nm=100.0, z_step_nm=300.0,
        signal_channel=0, dapi_channel=0, roi_output_dir=out,
    )
    saved = tifffile.imread(out / "img_ROI.tif")
    assert saved.shape == (6, 30, 30) and saved[3, 10, 10] == 1
    assert not (tmp_path / "img_ROI.tif").exists()


def _run_worker(folder, out):
    from App import AnalysisWorker

    worker = AnalysisWorker({
        "input_folder": str(folder), "output_folder": str(out), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": False, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
    })

    def answer(info):
        worker.user_roi_data = _roi(info["shape"], info["is_3d"])
        worker.roi_event.set()

    worker.request_roi_signal.connect(answer)
    worker.run()
    return json.loads(next(out.glob("*_metadata.json")).read_text(encoding="utf-8"))


def test_two_runs_keep_their_own_roi_and_record_it(tmp_path):
    from provenance import file_sha256

    _pair(tmp_path)
    meta1 = _run_worker(tmp_path, tmp_path / "Result")
    meta2 = _run_worker(tmp_path, tmp_path / "ResultNew")

    assert (tmp_path / "Result" / "img_ROI.tif").is_file()
    assert (tmp_path / "ResultNew" / "img_ROI.tif").is_file()
    assert not (tmp_path / "img_ROI.tif").exists()             # input folder untouched
    record = meta1["files"]["input_files"][0]
    assert record["roi_file"] == "img_ROI.tif"
    assert record["roi_sha256"] == file_sha256(tmp_path / "Result" / "img_ROI.tif")
    assert meta2["files"]["input_files"][0]["roi_sha256"] == file_sha256(tmp_path / "ResultNew" / "img_ROI.tif")


def test_leftover_roi_files_are_not_treated_as_source_images(tmp_path):
    from App import _source_tiff_files

    _pair(tmp_path)
    tifffile.imwrite(tmp_path / "img_ROI.tif", np.ones((6, 30, 30), dtype=np.uint8))
    assert [p.name for p in _source_tiff_files(tmp_path)] == ["img.tif"]
    meta = _run_worker(tmp_path, tmp_path / "out")
    assert meta["files"]["excluded_files"] == []
