import numpy as np

from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell


MODE_A_LAYER_SCHEME = "edt_over_max_thirds_v1"


def compute_core_shell_metrics(
    object_mask,
    *,
    sampling,
    intensity_image=None,
    min_core_voxels=20,
    primary_include=True,
    primary_exclusion_reason="",
):

    #Radial FA Profiling always receives physical sampling in Z/Y/X order from Batch.
    if len(sampling) != 3 or any(float(value) <= 0 for value in sampling):
        raise ValueError("sampling must contain three positive values ordered (Z, Y, X)")

    sampling = tuple(float(value) for value in sampling)
    layers = split_core_middle_shell(
        object_mask,
        min_core_voxels=min_core_voxels,
        sampling=sampling,
    )

    #Keep layer metrics and their validity flags together so Batch can copy the same schema into every object row.
    results = {}
    empty_layers = []
    object_principal_std = None
    for name, layer_mask in (
        ("object", object_mask),
        ("shell", layers["shell"]),
        ("middle", layers["middle"]),
        ("core", layers["core"]),
    ):
        if np.any(layer_mask):
            result = shape_anisotropy(
                layer_mask,
                sampling=sampling,
                min_voxels=min_core_voxels,
            )
            layer_mean_int = (
                float(np.mean(intensity_image[layer_mask]))
                if intensity_image is not None
                else np.nan
            )
            # Capture the object-level principal standard deviations for the null model.
            if name == "object":
                object_principal_std = result.get("principal_std")
        else:
            empty_layers.append(name)
            result = {
                "fractional_anisotropy": np.nan,
                "voxel_count": 0,
                "anisotropy_valid": False,
            }
            layer_mean_int = np.nan
        results[f"A_{name}"] = result["fractional_anisotropy"]
        results[f"{name}_voxels"] = result["voxel_count"]
        results[f"A_{name}_valid"] = result["anisotropy_valid"]
        results[f"mean_intensity_{name}"] = layer_mean_int

    #The delta is a geometric comparison and remains NaN when a layer is empty, preserving the QC signal instead of inventing a value.
    results["delta_A_middle_shell"] = results["A_middle"] - results["A_shell"]
    results["delta_A_core_middle"] = results["A_core"] - results["A_middle"]
    results["delta_A_core_shell"] = results["A_core"] - results["A_shell"]

    #Intensity differences across layers (Core vs Middle vs Shell)
    results["delta_intensity_middle_shell"] = results["mean_intensity_middle"] - results["mean_intensity_shell"]
    results["delta_intensity_core_middle"] = results["mean_intensity_core"] - results["mean_intensity_middle"]
    results["delta_intensity_core_shell"] = results["mean_intensity_core"] - results["mean_intensity_shell"]

    results["core_valid"] = bool(layers["qc"]["core_valid"])
    results["layer_qc"] = layers["qc"]
    results["layers"] = layers
    results["mode_a_empty_layers"] = ";".join(empty_layers)

    # Geometric null model: expected radial ΔFA for a homogeneous ellipsoid
    # with the same principal axes as this object.
    if object_principal_std is not None and len(object_principal_std) == 3:
        results["principal_std_z_nm"] = round(float(object_principal_std[0]), 4)
        results["principal_std_y_nm"] = round(float(object_principal_std[1]), 4)
        results["principal_std_x_nm"] = round(float(object_principal_std[2]), 4)
        try:
            from null_model import null_model_delta_fa
            null = null_model_delta_fa(object_principal_std, sampling, min_core_voxels)
            for key, value in null.items():
                results[key] = value
            # Biological excess: measured minus geometric expectation
            if np.isfinite(results["delta_A_core_shell"]) and np.isfinite(null["null_delta_FA_core_shell"]):
                results["delta_FA_excess"] = round(
                    results["delta_A_core_shell"] - null["null_delta_FA_core_shell"], 6
                )
            else:
                results["delta_FA_excess"] = np.nan
        except Exception:
            results["null_FA_object"] = np.nan
            results["null_FA_shell"] = np.nan
            results["null_FA_middle"] = np.nan
            results["null_FA_core"] = np.nan
            results["null_delta_FA_core_shell"] = np.nan
            results["null_valid"] = False
            results["delta_FA_excess"] = np.nan
    else:
        results["principal_std_z_nm"] = np.nan
        results["principal_std_y_nm"] = np.nan
        results["principal_std_x_nm"] = np.nan
        results["null_FA_object"] = np.nan
        results["null_FA_shell"] = np.nan
        results["null_FA_middle"] = np.nan
        results["null_FA_core"] = np.nan
        results["null_delta_FA_core_shell"] = np.nan
        results["null_valid"] = False
        results["delta_FA_excess"] = np.nan



    results["mode_a_sampling_z_nm"] = sampling[0]
    results["mode_a_sampling_y_nm"] = sampling[1]
    results["mode_a_sampling_x_nm"] = sampling[2]
    results["mode_a_sampling_order"] = "Z,Y,X"
    #FA is computed from voxel-coordinate PCA (shape), not from intensity gradients.
    results["mode_a_fa_type"] = "geometric_shape"
    results["mode_a_min_core_voxels"] = int(min_core_voxels)
    results["mode_a_primary_include"] = bool(primary_include)
    results["mode_a_qc_reason"] = str(primary_exclusion_reason)

    # Automatic Morphological Particle Classification
    particle_class = classify_particle(
        a_object=results["A_object"],
        delta_intensity_core_shell=results["delta_intensity_core_shell"],
        voxel_count=results["object_voxels"],
        min_core_voxels=min_core_voxels,
    )
    results["particle_class"] = particle_class
    results["condensate_class"] = particle_class  # backward compatibility alias
    return results


def classify_particle(
    a_object: float,
    delta_intensity_core_shell: float,
    voxel_count: int,
    min_core_voxels: int = 20,
) -> str:
    """Classify a 3D segmented particle into objective morphological categories based on shape and radial profile."""
    if not np.isfinite(a_object):
        return "Unclassified"
    if voxel_count < min_core_voxels * 3:
        return "Punctate / Small Cluster"
    if a_object < 0.65:
        if not np.isfinite(delta_intensity_core_shell):
            return "Globular (Low FA) / Profile-Unavailable"
        if abs(delta_intensity_core_shell) < 1e-3:
            return "Globular (Low FA) / No-Gradient"
        if delta_intensity_core_shell < 0:
            return "Globular (Low FA) / Shell-Enriched"
        return "Globular (Low FA) / Core-Enriched"
    else:
        if not np.isfinite(delta_intensity_core_shell):
            return "Elongated (High FA) / Profile-Unavailable"
        if abs(delta_intensity_core_shell) < 1e-3:
            return "Elongated (High FA) / No-Gradient"
        if delta_intensity_core_shell < 0:
            return "Elongated (High FA) / Shell-Enriched"
        return "Elongated (High FA) / Core-Enriched"


# Backward compatibility alias
classify_condensate = classify_particle