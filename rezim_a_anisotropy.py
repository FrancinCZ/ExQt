
import numpy as np


# Calculate geometric fractional anisotropy for one 2D or 3D mask.
def shape_anisotropy(mask, sampling=None, min_voxels=20):
    binary_mask = np.asarray(mask, dtype=bool)
    if binary_mask.ndim not in (2, 3):
        raise ValueError("Anisotropy mask must be 2D or 3D.")

    #Convert foreground voxels into physical coordinates so anisotropy uses the same Z/Y/X calibration supplied by Batch.process_condensates.
    coordinates = np.argwhere(binary_mask)
    voxel_count = int(coordinates.shape[0])
    if voxel_count == 0:
        raise ValueError("Mask contains no foreground voxels.")

    if sampling is None:
        scale = np.ones(binary_mask.ndim, dtype=float)
    else:
        scale = np.asarray(sampling, dtype=float)
        if scale.shape != (binary_mask.ndim,) or np.any(scale <= 0):
            raise ValueError("sampling must contain one positive value per mask axis.")

    physical_coordinates = coordinates * scale
    centered = physical_coordinates - physical_coordinates.mean(axis=0)

    #Eigenvalues of the coordinate covariance describe spread along the object's principal axes.
    covariance = centered.T @ centered / voxel_count
    eigenvalues = np.linalg.eigvalsh(covariance)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    eigenvalue_sum_of_squares = float(np.sum(eigenvalues**2))

    if eigenvalue_sum_of_squares == 0.0:
        fractional_anisotropy = np.nan
        principal_std_nm = np.full(binary_mask.ndim, np.nan)
    else:
        eigenvalue_mean = float(np.mean(eigenvalues))
        fractional_anisotropy = float(
            np.sqrt(
                binary_mask.ndim
                / (binary_mask.ndim - 1)
                * np.sum((eigenvalues - eigenvalue_mean) ** 2)
                / eigenvalue_sum_of_squares
            )
        )
        principal_std_nm = np.sqrt(eigenvalues)

    return {
        "fractional_anisotropy": fractional_anisotropy,
        "principal_std": principal_std_nm,
        "voxel_count": voxel_count,
        "anisotropy_valid": (
            voxel_count >= min_voxels and np.isfinite(fractional_anisotropy)
        ),
    }
