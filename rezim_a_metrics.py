import numpy as np
from scipy import stats

import null_model
from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell
from reference_values import reference_value


MODE_A_LAYER_SCHEME = "edt_over_max_thirds_v1"
#Two-sided Welch t-test of core vs shell voxels, alpha = 0.01.
GRADIENT_TEST = "welch_t_two_sided"
GRADIENT_ALPHA = 0.01
UNCLASSIFIED_QC = "Unclassified (QC)"


#Welch t-test of core vs shell voxel intensities; returns (t, p) or (nan, nan) when not testable.
def _core_shell_gradient_test(core_values, shell_values):
    core = np.asarray(core_values, dtype=float)
    shell = np.asarray(shell_values, dtype=float)
    if core.size < 2 or shell.size < 2:
        return np.nan, np.nan
    if np.var(core) == 0 and np.var(shell) == 0:
        # Noise-free input: any difference is exact, no difference is no gradient.
        difference = float(np.mean(core) - np.mean(shell))
        return (float(np.sign(difference)) * np.inf, 0.0) if difference != 0 else (0.0, 1.0)
    result = stats.ttest_ind(core, shell, equal_var=False)
    return float(result.statistic), float(result.pvalue)


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
    object_principal_axes = None
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
                object_principal_axes = result.get("principal_axes")
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

    results["delta_intensity_middle_shell"] = results["mean_intensity_middle"] - results["mean_intensity_shell"]
    results["delta_intensity_core_middle"] = results["mean_intensity_core"] - results["mean_intensity_middle"]
    results["delta_intensity_core_shell"] = results["mean_intensity_core"] - results["mean_intensity_shell"]

    #The class needs a noise-aware test: the raw difference alone is random in sign for a homogeneous object.
    if intensity_image is not None and np.any(layers["core"]) and np.any(layers["shell"]):
        t_value, p_value = _core_shell_gradient_test(
            intensity_image[layers["core"]], intensity_image[layers["shell"]]
        )
    else:
        t_value, p_value = np.nan, np.nan
    results["gradient_t_core_shell"] = t_value
    results["gradient_p_core_shell"] = p_value
    results["gradient_significant"] = bool(np.isfinite(p_value) and p_value < GRADIENT_ALPHA)
    results["gradient_test"] = GRADIENT_TEST
    results["gradient_alpha"] = GRADIENT_ALPHA

    results["core_valid"] = bool(layers["qc"]["core_valid"])
    results["layer_qc"] = layers["qc"]
    results["layers"] = layers
    results["mode_a_empty_layers"] = ";".join(empty_layers)

    #Geometric null model: expected radial ΔFA for a homogeneous ellipsoid with the same
    #principal stds AND principal directions as this object, on the same anisotropic grid.
    #principal_std_1..3 are ascending (eigh order); they are NOT the Z/Y/X axes.
    if object_principal_std is not None and len(object_principal_std) == 3:
        for rank, value in enumerate(object_principal_std, start=1):
            results[f"principal_std_{rank}_nm"] = round(float(value), 4)
        # Expected failures come back as null_valid=False with null_invalid_reason;
        # anything else is a bug and must not be hidden.
        null = null_model.null_model_delta_fa(
            object_principal_std, sampling, min_core_voxels, principal_axes=object_principal_axes
        )
        results.update(null)
        # Excess: measured minus geometric expectation (mask shape only).
        if np.isfinite(results["delta_A_core_shell"]) and np.isfinite(null["null_delta_FA_core_shell"]):
            results["delta_FA_excess"] = round(
                results["delta_A_core_shell"] - null["null_delta_FA_core_shell"], 6
            )
        else:
            results["delta_FA_excess"] = np.nan
    else:
        for rank in (1, 2, 3):
            results[f"principal_std_{rank}_nm"] = np.nan
        results.update(null_model._nan_result("object_anisotropy_unavailable"))
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

    particle_class = classify_particle(
        a_object=results["A_object"],
        delta_intensity_core_shell=results["delta_intensity_core_shell"],
        voxel_count=results["object_voxels"],
        min_core_voxels=min_core_voxels,
        gradient_significant=results["gradient_significant"] if np.isfinite(p_value) else None,
        qc_passed=bool(primary_include),
    )
    results["particle_class"] = particle_class
    results["condensate_class"] = particle_class  # backward compatibility alias
    return results


def classify_particle(
    a_object: float,
    delta_intensity_core_shell: float,
    voxel_count: int,
    min_core_voxels: int = 20,
    gradient_significant: bool | None = None,
    qc_passed: bool = True,
) -> str:
    """Classify a 3D segmented particle into objective morphological categories based on shape and radial profile.

    Core-/Shell-Enriched requires a significant core vs shell difference (see GRADIENT_TEST,
    GRADIENT_ALPHA); None = not testable. Core-Enriched is also the expected result of PSF
    blurring of a homogeneous object, so it is not evidence of biological enrichment on its own.
    """
    if not qc_passed:
        return UNCLASSIFIED_QC
    if not np.isfinite(a_object):
        return "Unclassified"
    if voxel_count < min_core_voxels * 3:
        return "Punctate / Small Cluster"
    # Boundary is a working value [UNVERIFIED], defined once in reference_values.py.
    shape = "Globular (Low FA)" if a_object < reference_value("fa_globular_boundary") else "Elongated (High FA)"
    if not np.isfinite(delta_intensity_core_shell) or gradient_significant is None:
        return f"{shape} / Profile-Unavailable"
    if not gradient_significant:
        return f"{shape} / No Significant Gradient"
    if delta_intensity_core_shell < 0:
        return f"{shape} / Shell-Enriched"
    return f"{shape} / Core-Enriched"


# Backward compatibility alias
classify_condensate = classify_particle