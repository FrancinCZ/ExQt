# A REVIEW alignment only excludes a stack when a weak step lies inside the ROI's Z range
# (2026-11-07: weak steps Z125-Z132 above the nucleus excluded a whole SON stack).
import json

import numpy as np
import pandas as pd

from Batch import alignment_problem_steps, alignment_status_for_roi
from postprocessing import _prepare_reporting_frames


def _drift(folder, bad_steps, n_slices=10):
    rows = []
    for z in range(1, n_slices):
        status, reason = ("REVIEW", "weak_post_correlation") if z in bad_steps else ("PASS", "")
        rows.append({"reference_z": z - 1, "moving_z": z, "step_shift_y_px": 0.0, "step_shift_x_px": 0.0,
                     "step_magnitude_px": 0.0, "status": status, "reason": reason, "overall_status": "REVIEW"})
    pd.DataFrame(rows).to_csv(folder / "sample_drift.csv", index=False)
    return folder / "sample.tif"


def _roi(z_min, z_max, n_slices=10):
    roi = np.zeros((n_slices, 8, 8), dtype=bool)
    roi[z_min:z_max + 1, 2:6, 2:6] = True
    return roi


def test_problem_steps_are_read_from_drift_csv(tmp_path):
    tif = _drift(tmp_path, {8, 9})
    assert [z for z, _ in alignment_problem_steps(tif)] == [8, 9]


def test_steps_outside_roi_do_not_exclude(tmp_path):
    tif = _drift(tmp_path, {8, 9})
    assert alignment_status_for_roi(tif, _roi(1, 6)) == "REVIEW_OUTSIDE_ROI"


def test_step_inside_roi_keeps_review(tmp_path):
    tif = _drift(tmp_path, {4})
    assert alignment_status_for_roi(tif, _roi(1, 6)) == "REVIEW"


def test_step_at_roi_border_only_touches_one_roi_slice(tmp_path):
    # Step Z6->Z7 moves slice 7, which lies outside a ROI ending at Z6.
    tif = _drift(tmp_path, {7})
    assert alignment_status_for_roi(tif, _roi(1, 6)) == "REVIEW_OUTSIDE_ROI"
    # Step Z1->Z2 lies between two ROI slices.
    tif = _drift(tmp_path, {2})
    assert alignment_status_for_roi(tif, _roi(1, 6)) == "REVIEW"


def test_whole_fov_roi_keeps_review(tmp_path):
    tif = _drift(tmp_path, {8})
    assert alignment_status_for_roi(tif, np.ones((10, 8, 8), dtype=bool)) == "REVIEW"
    assert alignment_status_for_roi(tif, np.zeros((10, 8, 8), dtype=bool)) == "REVIEW"


def _objects(status):
    return pd.DataFrame({
        "filename": "a.tif", "cell_id": 1, "object_id": [1, 2], "mode": "3d", "is_3d": True,
        "volume_px": 600, "volume_bio_um3": [0.5, 0.6], "partition_coefficient": [3.0, 4.0],
        "touches_image_edge": False, "touches_roi_edge": False, "alignment_status": status,
    })


def test_review_outside_roi_stays_primary():
    _, primary, _, _, _ = _prepare_reporting_frames(_objects("REVIEW_OUTSIDE_ROI"), 0.1, 2.0)
    assert len(primary) == 2
    _, primary, excluded, _, _ = _prepare_reporting_frames(_objects("REVIEW"), 0.1, 2.0)
    assert len(primary) == 0 and len(excluded) == 2


def test_drift_csv_without_step_positions_stays_strict(tmp_path):
    pd.DataFrame({"status": ["PASS", "REVIEW"], "step_shift_y_px": 0.0, "step_shift_x_px": 0.0,
                  "step_magnitude_px": 0.0}).to_csv(tmp_path / "sample_drift.csv", index=False)
    assert alignment_status_for_roi(tmp_path / "sample.tif", _roi(1, 6)) == "REVIEW"
