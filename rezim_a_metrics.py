
import numpy as np

from rezim_a_anisotropy import shape_anisotropy
from rezim_a_core_shell import split_core_middle_shell


MODE_A_LAYER_SCHEME = "edt_over_max_thirds_v1"


def compute_core_shell_metrics(
    object_mask,
    *,
    sampling,
    min_core_voxels=20,
    primary_include=True,
    primary_exclusion_reason="",
):

    #Rezim A always receives physical sampling in Z/Y/X order from Batch.
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
        else:
            empty_layers.append(name)
            result = {
                "fractional_anisotropy": np.nan,
                "voxel_count": 0,
                "anisotropy_valid": False,
            }
        results[f"A_{name}"] = result["fractional_anisotropy"]
        results[f"{name}_voxels"] = result["voxel_count"]
        results[f"A_{name}_valid"] = result["anisotropy_valid"]

    #The delta is a geometric comparison and remains NaN when a layer is empty, preserving the QC signal instead of inventing a value.
    results["delta_A_middle_shell"] = results["A_middle"] - results["A_shell"]
    results["delta_A_core_middle"] = results["A_core"] - results["A_middle"]
    results["delta_A_core_shell"] = results["A_core"] - results["A_shell"]
    results["core_valid"] = bool(layers["qc"]["core_valid"])
    results["layer_qc"] = layers["qc"]
    results["layers"] = layers
    results["mode_a_empty_layers"] = ";".join(empty_layers)


    results["mode_a_sampling_z_nm"] = sampling[0]
    results["mode_a_sampling_y_nm"] = sampling[1]
    results["mode_a_sampling_x_nm"] = sampling[2]
    results["mode_a_sampling_order"] = "Z,Y,X"
    results["mode_a_min_core_voxels"] = int(min_core_voxels)
    results["mode_a_primary_include"] = bool(primary_include)
    results["mode_a_qc_reason"] = str(primary_exclusion_reason)
    return results
