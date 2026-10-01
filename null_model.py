import numpy as np
from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell


def _make_oriented_ellipsoid(semi_axes_nm, principal_axes, sampling):
    """Binary ellipsoid with the given physical semi-axes along the given principal directions,
    sampled on the same anisotropic (Z, Y, X) grid as the measured object."""
    semi_axes_nm = np.asarray(semi_axes_nm, dtype=float)
    rotation = np.asarray(principal_axes, dtype=float)      
    sampling = np.asarray(sampling, dtype=float)
    #Exact half-extent of the rotated ellipsoid along each grid axis, plus 2 voxels padding.
    extent_nm = np.sqrt(np.sum((rotation * semi_axes_nm) ** 2, axis=1))
    half = np.ceil(extent_nm / sampling).astype(int) + 2
    grids = np.meshgrid(
        *[np.arange(-h, h + 1) * step for h, step in zip(half, sampling)], indexing="ij"
    )
    points = np.stack([g.ravel() for g in grids], axis=1)
    body = points @ rotation                                 # coordinates in the object's frame
    inside = np.sum((body / semi_axes_nm) ** 2, axis=1) <= 1.0
    return inside.reshape(grids[0].shape)


def null_model_delta_fa(principal_std_nm, sampling, min_core_voxels=20, *, principal_axes):
    """Expected radial ΔFA of a homogeneous ellipsoid with the object's principal stds AND
    principal directions (eigenvectors), discretised on the object's anisotropic grid."""
    principal_std_nm = np.asarray(principal_std_nm, dtype=float)
    principal_axes = np.asarray(principal_axes, dtype=float)
    sampling = np.asarray(sampling, dtype=float)

    if not np.all(np.isfinite(principal_std_nm)) or np.any(principal_std_nm <= 0):
        return _nan_result("non_finite_principal_std")
    if principal_axes.shape != (3, 3) or not np.all(np.isfinite(principal_axes)):
        return _nan_result("non_finite_principal_axes")

    #A homogeneous solid ellipsoid has std = semi-axis / sqrt(5) along each principal axis.
    semi_axes_nm = principal_std_nm * np.sqrt(5.0)
    ellipsoid = _make_oriented_ellipsoid(semi_axes_nm, principal_axes, sampling)

    if ellipsoid.sum() < min_core_voxels * 3:
        return _nan_result("too_few_voxels")

    try:
        layers = split_core_middle_shell(ellipsoid, min_core_voxels=min_core_voxels, sampling=tuple(sampling))
    except ValueError:
        return _nan_result("layer_split_failed")

    fa_values = {}
    layer_valid = {}
    for name, mask in [("object", ellipsoid),
                        ("shell", layers["shell"]),
                        ("middle", layers["middle"]),
                        ("core", layers["core"])]:
        if np.any(mask):
            result = shape_anisotropy(mask, sampling=tuple(sampling), min_voxels=min_core_voxels)
            fa_values[name] = result["fractional_anisotropy"]
            layer_valid[name] = bool(result["anisotropy_valid"])
        else:
            fa_values[name] = np.nan
            layer_valid[name] = False

    delta = fa_values["core"] - fa_values["shell"]
    valid = all(layer_valid[k] for k in ("shell", "middle", "core"))

    return {
        "null_FA_object": round(fa_values["object"], 6) if np.isfinite(fa_values["object"]) else np.nan,
        "null_FA_shell": round(fa_values["shell"], 6) if np.isfinite(fa_values["shell"]) else np.nan,
        "null_FA_middle": round(fa_values["middle"], 6) if np.isfinite(fa_values["middle"]) else np.nan,
        "null_FA_core": round(fa_values["core"], 6) if np.isfinite(fa_values["core"]) else np.nan,
        "null_delta_FA_core_shell": round(delta, 6) if np.isfinite(delta) else np.nan,
        "null_valid": valid,
        "null_invalid_reason": "" if valid else "invalid_layer_fa",
    }


def _nan_result(reason):
    return {
        "null_FA_object": np.nan,
        "null_FA_shell": np.nan,
        "null_FA_middle": np.nan,
        "null_FA_core": np.nan,
        "null_delta_FA_core_shell": np.nan,
        "null_valid": False,
        "null_invalid_reason": reason,
    }
