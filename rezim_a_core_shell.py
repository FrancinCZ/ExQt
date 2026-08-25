import numpy as np
from scipy.ndimage import distance_transform_edt


def split_core_middle_shell(
    object_mask,
    min_core_voxels=20,
    sampling=None,
    shell_end=1 / 3,
    core_start=2 / 3,
):

    mask = np.asarray(object_mask, dtype=bool)

    if mask.ndim not in (2, 3):
        raise ValueError("object_mask must be 2D or 3D.")

    if not np.any(mask):
        raise ValueError("object_mask contains no foreground object.")

    if not 0 < shell_end < core_start < 1:
        raise ValueError("Layer boundaries must satisfy 0 < shell_end < core_start < 1.")

    #Distance from the boundary gives a size-normalized radial coordinate that works with anisotropic voxel spacing when sampling is provided.
    distance_from_edge = distance_transform_edt(mask, sampling=sampling)
    max_distance = float(distance_from_edge.max())

    if max_distance == 0:
        raise ValueError("Object is too small to form radial layers.")


    normalized_distance = distance_from_edge / max_distance

    #Fixed thirds make the layer scheme reproducible across objects. Batch records the exact versioned EDT/max definition with every output row.
    shell = mask & (normalized_distance <= shell_end)
    middle = mask & (normalized_distance > shell_end) & (normalized_distance <= core_start)
    core = mask & (normalized_distance > core_start)

    core_voxels = int(core.sum())
    complete_coverage = int(shell.sum() + middle.sum() + core.sum()) == int(mask.sum())

    return {
        "core": core,
        "middle": middle,
        "shell": shell,
        "normalized_distance": normalized_distance,
        "qc": {
            "object_voxels": int(mask.sum()),
            "core_voxels": core_voxels,
            "max_distance": max_distance,
            "core_valid": core_voxels >= min_core_voxels,
            "shell_end": shell_end,
            "core_start": core_start,
            "complete_coverage": complete_coverage,
        },
    }


if __name__ == "__main__":
    print("This module provides split_core_middle_shell().")
    print("Next step: test it on one real region.image from a Labkit mask.")
