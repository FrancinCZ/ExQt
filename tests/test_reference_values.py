# S6 (review 2026-09-23): reference values used by plots and classification live in one
# place, carry their source/status, and plot titles do not present interpretation as a test.
import inspect
import json
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from reference_values import REFERENCE_VALUES, reference_value


def _fa_prolate(ratio):
    eig = np.array([ratio ** 2, 1.0, 1.0])
    return float(np.sqrt(1.5 * np.sum((eig - eig.mean()) ** 2) / np.sum(eig ** 2)))


def test_every_reference_value_has_value_meaning_and_status():
    assert set(REFERENCE_VALUES) == {"fa_globular_boundary", "fa_psf_elongation_limit", "k_low_contrast"}
    for entry in REFERENCE_VALUES.values():
        assert np.isfinite(entry["value"])
        assert entry["meaning"] and entry["status"]
    assert all("UNVERIFIED" in REFERENCE_VALUES[k]["status"] for k in REFERENCE_VALUES)


def test_psf_limit_description_matches_its_geometry():
    # The stated axis-ratio equivalents must be true (FA of a homogeneous prolate ellipsoid).
    assert _fa_prolate(2.5) == pytest.approx(reference_value("fa_psf_elongation_limit"), abs=0.005)
    assert _fa_prolate(1.84) == pytest.approx(reference_value("fa_globular_boundary"), abs=0.005)


def test_classification_uses_the_central_boundary(monkeypatch):
    import reference_values
    from rezim_a_metrics import classify_particle

    kwargs = dict(delta_intensity_core_shell=10.0, voxel_count=5000, gradient_significant=False)
    assert classify_particle(0.64, **kwargs).startswith("Globular")
    assert classify_particle(0.66, **kwargs).startswith("Elongated")
    monkeypatch.setitem(reference_values.REFERENCE_VALUES["fa_globular_boundary"], "value", 0.70)
    assert classify_particle(0.66, **kwargs).startswith("Globular")


def test_plots_use_central_values_and_neutral_titles():
    import partitioning_plots

    source = inspect.getsource(partitioning_plots)
    for loaded in ("Phase Coexistence", "Lever Rule", "Surface Relaxation", "0.82", "0.65", "1.5 /", "< 1.5"):
        assert loaded not in source, loaded
    assert "reference_value(" in source


def test_reference_values_are_recorded_in_metadata_and_fingerprint(tmp_path):
    import pandas as pd
    import tifffile
    from App import AnalysisWorker
    from postprocessing import qc_policy_fingerprint

    raw = np.full((6, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros((6, 30, 30), dtype=np.uint8)
    mask[2:4, 10:14, 10:14] = 1
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
    }).run()
    meta = json.loads(next((tmp_path / "out").glob("*_metadata.json")).read_text(encoding="utf-8"))
    assert meta["reference_values"]["fa_globular_boundary"]["value"] == 0.65

    df = pd.DataFrame({"volume_bio_um3": [0.1]})
    changed = json.loads(json.dumps(meta))
    changed["reference_values"]["fa_globular_boundary"]["value"] = 0.7
    assert qc_policy_fingerprint(meta, df, 0, 1)[0] != qc_policy_fingerprint(changed, df, 0, 1)[0]
