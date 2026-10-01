# Regression locks for results that review 2026-09-23 verified once but no test covered
# (VALIDATION_PROTOCOL §5). The documentation refers to these tests.
import numpy as np
import pytest

from Batch import _assess_z_split_topology
from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell


def _ellipsoid(semi_axes_nm, sampling):
    s = np.asarray(sampling, dtype=float)
    a = np.asarray(semi_axes_nm, dtype=float)
    half = np.ceil(a / s).astype(int) + 2
    z, y, x = np.meshgrid(*[np.arange(-h, h + 1) * step for h, step in zip(half, s)], indexing="ij")
    return (z / a[0]) ** 2 + (y / a[1]) ** 2 + (x / a[2]) ** 2 <= 1.0


ANISO = (60.0, 15.0, 15.0)  # Z step 4x coarser than XY, as in ExM


def _fa(semi_axes):
    """Analytic FA of a homogeneous ellipsoid (covariance eigenvalues = a_i^2 / 5)."""
    eig = np.asarray(semi_axes, dtype=float) ** 2 / 5.0
    return float(np.sqrt(1.5 * np.sum((eig - eig.mean()) ** 2) / np.sum(eig ** 2)))


@pytest.mark.parametrize("semi_axes", [(600, 600, 600), (1200, 600, 600), (1800, 600, 300), (400, 900, 600)])
def test_fa_matches_analytic_value_on_anisotropic_grid(semi_axes):
    result = shape_anisotropy(_ellipsoid(semi_axes, ANISO), sampling=ANISO)
    assert result["fractional_anisotropy"] == pytest.approx(_fa(semi_axes), abs=0.01)
    assert result["principal_std"] == pytest.approx(np.sort(semi_axes) / np.sqrt(5), rel=0.02)


def test_fa_without_sampling_would_be_a_false_anisotropy():
    """rules.md §2: the same sphere measured in voxel units looks elongated (FA ≈ 0.66)."""
    sphere = _ellipsoid((600, 600, 600), ANISO)
    assert shape_anisotropy(sphere, sampling=ANISO)["fractional_anisotropy"] < 0.01
    assert shape_anisotropy(sphere)["fractional_anisotropy"] > 0.6


def test_sphere_layer_fractions_match_radial_thirds():
    """Theory for EDT thirds of a sphere: shell 1-(2/3)^3, middle (2/3)^3-(1/3)^3, core (1/3)^3."""
    sphere = _ellipsoid((600, 600, 600), (15.0, 15.0, 15.0))
    layers = split_core_middle_shell(sphere, sampling=(15.0, 15.0, 15.0))
    total = sphere.sum()
    fractions = [layers[name].sum() / total for name in ("shell", "middle", "core")]
    assert fractions == pytest.approx([19 / 27, 7 / 27, 1 / 27], abs=0.02)
    assert (layers["shell"] | layers["middle"] | layers["core"]).sum() == total
    assert not (layers["shell"] & layers["core"]).any()


def test_z_topology_single_blob_passes_and_split_object_does_not():
    blob = np.zeros((6, 20, 20), dtype=bool)
    blob[1:5, 5:15, 5:15] = True
    assert _assess_z_split_topology(blob)["status"] == "pass"

    split = np.zeros((6, 30, 30), dtype=bool)
    split[1:5, 5:12, 5:12] = True                     # two substantial parts on every slice,
    split[1:5, 18:25, 18:25] = True                   # joined only on one slice
    split[2, 11:19, 11:19] = True
    assert _assess_z_split_topology(split)["status"] in {"review", "fail"}
