"""
null_model.py — Geometric null model for radial FA layering.

For each measured 3D object, generates a synthetic homogeneous ellipsoid with
the same principal semi-axes (derived from the covariance eigenvalues) and the
same voxel sampling, then computes FA_shell / FA_middle / FA_core using the
identical EDT-thirds algorithm as the real analysis.

This isolates the purely geometric component of the radial FA gradient so that
the biological excess can be quantified as:

    delta_FA_excess = delta_FA_measured - delta_FA_null

If delta_FA_excess ≈ 0 for most objects, the observed radial gradient is
entirely explained by ellipsoidal geometry.  If delta_FA_excess > 0
consistently, the biological structure (e.g. SON fibrillar core) produces
additional anisotropy beyond what geometry alone predicts.
"""

import numpy as np
from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell


def _make_ellipsoid(semi_axes_voxels, sampling):
    """Build a binary 3D ellipsoid mask with the given semi-axes (in voxels)."""
    radii = np.asarray(semi_axes_voxels, dtype=float)
    half = np.ceil(radii).astype(int) + 2          # pad by 2 voxels
    z, y, x = np.mgrid[-half[0]:half[0]+1,
                        -half[1]:half[1]+1,
                        -half[2]:half[2]+1]
    inside = (z / radii[0])**2 + (y / radii[1])**2 + (x / radii[2])**2 <= 1.0
    return inside


def null_model_delta_fa(principal_std_nm, sampling, min_core_voxels=20):
    """
    Given the principal standard deviations (in nm, Z/Y/X order) of a real
    object's voxel distribution, construct a homogeneous ellipsoid with the
    same shape and compute its radial FA layer values.

    Parameters
    ----------
    principal_std_nm : array-like, shape (3,)
        Standard deviations along the three principal axes of the real object
        in physical units (nm), as returned by shape_anisotropy()['principal_std'].
    sampling : tuple of 3 floats
        Physical voxel size (Z, Y, X) in nm.
    min_core_voxels : int
        Minimum voxels required for valid FA.

    Returns
    -------
    dict with keys:
        null_FA_shell, null_FA_middle, null_FA_core, null_FA_object,
        null_delta_FA_core_shell, null_valid
    """
    principal_std_nm = np.asarray(principal_std_nm, dtype=float)
    sampling = np.asarray(sampling, dtype=float)

    if not np.all(np.isfinite(principal_std_nm)) or np.any(principal_std_nm <= 0):
        return _nan_result()

    # Convert physical standard deviations to semi-axes in voxel units.
    # For a uniform ellipsoid with semi-axis a, the standard deviation of the
    # uniform distribution along that axis is a / sqrt(5).
    # So: a = std * sqrt(5).
    semi_axes_nm = principal_std_nm * np.sqrt(5.0)
    semi_axes_voxels = semi_axes_nm / sampling

    # Ensure the ellipsoid is large enough to produce meaningful layers.
    if np.any(semi_axes_voxels < 2):
        return _nan_result()

    ellipsoid = _make_ellipsoid(semi_axes_voxels, sampling)

    if ellipsoid.sum() < min_core_voxels * 3:
        return _nan_result()

    try:
        layers = split_core_middle_shell(ellipsoid, min_core_voxels=min_core_voxels, sampling=tuple(sampling))
    except ValueError:
        return _nan_result()

    fa_values = {}
    for name, mask in [("object", ellipsoid),
                       ("shell", layers["shell"]),
                       ("middle", layers["middle"]),
                       ("core", layers["core"])]:
        if np.any(mask):
            result = shape_anisotropy(mask, sampling=tuple(sampling), min_voxels=min_core_voxels)
            fa_values[name] = result["fractional_anisotropy"]
        else:
            fa_values[name] = np.nan

    delta = fa_values["core"] - fa_values["shell"]
    valid = all(np.isfinite(fa_values[k]) for k in ("shell", "middle", "core"))

    return {
        "null_FA_object": round(fa_values["object"], 6) if np.isfinite(fa_values["object"]) else np.nan,
        "null_FA_shell": round(fa_values["shell"], 6) if np.isfinite(fa_values["shell"]) else np.nan,
        "null_FA_middle": round(fa_values["middle"], 6) if np.isfinite(fa_values["middle"]) else np.nan,
        "null_FA_core": round(fa_values["core"], 6) if np.isfinite(fa_values["core"]) else np.nan,
        "null_delta_FA_core_shell": round(delta, 6) if np.isfinite(delta) else np.nan,
        "null_valid": valid,
    }


def _nan_result():
    return {
        "null_FA_object": np.nan,
        "null_FA_shell": np.nan,
        "null_FA_middle": np.nan,
        "null_FA_core": np.nan,
        "null_delta_FA_core_shell": np.nan,
        "null_valid": False,
    }
