import numpy as np


#Calculate geometric fractional anisotropy for one 2D or 3D mask.
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
    # eigh returns ascending eigenvalues; column k of eigenvectors belongs to eigenvalue k.
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    eigenvalue_sum_of_squares = float(np.sum(eigenvalues**2))

    if eigenvalue_sum_of_squares == 0.0:
        fractional_anisotropy = np.nan
        principal_std_nm = np.full(binary_mask.ndim, np.nan)
        principal_axes = np.full((binary_mask.ndim, binary_mask.ndim), np.nan)
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
        principal_axes = eigenvectors

    #Shape in the XY plane only and the direction of the longest axis: an elongation along Z
    #(PSF, axial scaling, anisotropic expansion) does not change FA_xy, and shows up as a small
    #angle between the main axis and Z (0 deg = along Z, 90 deg = in the XY plane).
    fa_xy = np.nan
    angle_to_z_deg = np.nan
    if binary_mask.ndim == 3 and eigenvalue_sum_of_squares > 0.0:
        xy_eigenvalues = np.clip(np.linalg.eigvalsh(covariance[1:, 1:]), 0.0, None)
        xy_sum_of_squares = float(np.sum(xy_eigenvalues**2))
        if xy_sum_of_squares > 0.0:
            fa_xy = float(np.sqrt(2.0 * np.sum((xy_eigenvalues - xy_eigenvalues.mean()) ** 2) / xy_sum_of_squares))
        main_axis = eigenvectors[:, -1]
        angle_to_z_deg = float(np.degrees(np.arccos(np.clip(abs(main_axis[0]), 0.0, 1.0))))

    return {
        "fractional_anisotropy": fractional_anisotropy,
        "fa_xy": fa_xy,
        "main_axis_angle_to_z_deg": angle_to_z_deg,
        "principal_std": principal_std_nm,
        #Columns are unit principal directions in physical (Z, Y, X) space, ascending by std.
        "principal_axes": principal_axes,
        "voxel_count": voxel_count,
        "anisotropy_valid": (
            voxel_count >= min_voxels and np.isfinite(fractional_anisotropy)
        ),
    }