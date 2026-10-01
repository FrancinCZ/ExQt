# K5 (review 2026-09-23): the geometric null ellipsoid must share the object's
# principal axes (eigenvectors), not be axis-aligned with the smallest axis on Z.
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from rezim_a_metrics import compute_core_shell_metrics

SAMPLING = (62.5, 14.5, 14.5)  # Z, Y, X in nm, as in the review


def _rotated_ellipsoid(semi_axes_nm, rotation, sampling=SAMPLING):
    """Homogeneous ellipsoid, arbitrary orientation, sampled on the anisotropic grid."""
    a = np.asarray(semi_axes_nm, dtype=float)
    s = np.asarray(sampling, dtype=float)
    half = np.ceil(a.max() / s).astype(int) + 2
    grids = np.meshgrid(*[np.arange(-h, h + 1) * step for h, step in zip(half, s)], indexing="ij")
    points = np.stack([g.ravel() for g in grids], axis=1)      # physical Z, Y, X
    body = points @ rotation.as_matrix()                         # coordinates in the ellipsoid frame
    inside = np.sum((body / a) ** 2, axis=1) <= 1.0
    return inside.reshape(grids[0].shape)


def _metrics(mask):
    return compute_core_shell_metrics(mask, sampling=SAMPLING, intensity_image=None, min_core_voxels=20)


def test_random_homogeneous_ellipsoids_have_no_excess():
    """Acceptance: 200 rotated homogeneous ellipsoids, true excess = 0."""
    rng = np.random.default_rng(0)
    excess, null_invalid = [], 0
    for i in range(200):
        axes = rng.uniform(100.0, 400.0, size=3)
        mask = _rotated_ellipsoid(axes, Rotation.random(random_state=i))
        m = _metrics(mask)
        if not m["null_valid"]:
            null_invalid += 1
        if np.isfinite(m["delta_FA_excess"]):
            excess.append(m["delta_FA_excess"])
    excess = np.asarray(excess)
    q1, q3 = np.percentile(excess, [25, 75])
    assert np.median(np.abs(excess)) < 0.005
    assert q3 - q1 < 0.01
    assert null_invalid / 200 < 0.02


def test_null_ellipsoid_matches_object_shape_for_elongated_object():
    # Long axis in XY, tilted: the old axis-aligned null put the shortest axis on Z.
    mask = _rotated_ellipsoid((120.0, 200.0, 380.0), Rotation.from_euler("zyx", [30, 20, 10], degrees=True))
    m = _metrics(mask)
    assert m["null_valid"]
    assert m["null_FA_object"] == pytest.approx(m["A_object"], abs=0.01)
    assert abs(m["delta_FA_excess"]) < 0.02


def test_principal_std_columns_are_ascending_and_named_by_rank():
    mask = _rotated_ellipsoid((120.0, 200.0, 380.0), Rotation.random(random_state=3))
    m = _metrics(mask)
    stds = [m["principal_std_1_nm"], m["principal_std_2_nm"], m["principal_std_3_nm"]]
    assert stds == sorted(stds)
    # Homogeneous ellipsoid: std = a / sqrt(5).
    assert stds == pytest.approx(np.array([120.0, 200.0, 380.0]) / np.sqrt(5), rel=0.05)
    for old in ("principal_std_z_nm", "principal_std_y_nm", "principal_std_x_nm"):
        assert old not in m


def test_null_invalid_reason_is_recorded():
    from null_model import null_model_delta_fa

    bad = null_model_delta_fa([np.nan, 1.0, 1.0], SAMPLING, principal_axes=np.eye(3))
    assert not bad["null_valid"] and bad["null_invalid_reason"] == "non_finite_principal_std"

    tiny = null_model_delta_fa([20.0, 25.0, 30.0], SAMPLING, principal_axes=np.eye(3))
    assert not tiny["null_valid"] and tiny["null_invalid_reason"] == "too_few_voxels"

    ok = _metrics(_rotated_ellipsoid((150.0, 250.0, 350.0), Rotation.random(random_state=5)))
    assert ok["null_valid"] and ok["null_invalid_reason"] == ""


def test_unexpected_null_model_errors_are_not_swallowed(monkeypatch):
    import null_model

    def boom(*args, **kwargs):
        raise RuntimeError("bug in null model")

    monkeypatch.setattr(null_model, "null_model_delta_fa", boom)
    with pytest.raises(RuntimeError, match="bug in null model"):
        _metrics(_rotated_ellipsoid((150.0, 250.0, 350.0), Rotation.random(random_state=5)))
