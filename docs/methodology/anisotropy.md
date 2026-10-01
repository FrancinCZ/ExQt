# 3D Fractional Anisotropy (FA)

Fractional Anisotropy (FA) is a number between **0** and **1** that describes how elongated the 3D mask of an object is:

| FA | Axis ratio of a prolate ellipsoid (a:1:1) |
| :--- | :---: |
| 0.0 | 1.0 (sphere) |
| 0.47 | 1.5 |
| 0.71 | 2.0 |
| 0.82 | 2.5 |
| → 1 | → ∞ |

## What Affects FA

FA describes the segmented mask, not the cause of its shape. Besides the structure itself it is affected by:

- **optics** – the confocal PSF is elongated along Z, so small objects appear elongated;
- **segmentation** – where the mask boundary lies;
- **expansion and axial scale** – non-uniform expansion or an incorrect Z-step change the shape;
- **discretisation** – objects with few voxels give noisy values.

## How It Is Calculated

ExQt takes the physical coordinates of all mask voxels (voxel index × sampling in nm, order Z, Y, X), computes their covariance matrix and its eigenvalues $\lambda_1 \le \lambda_2 \le \lambda_3$ (variances along the object's own principal axes, not along X/Y/Z):

$$\text{FA} = \sqrt{\frac{3}{2}} \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$

If all three eigenvalues are equal, FA = 0. The physical sampling is essential: the same sphere measured in voxel units on a grid with a 4× coarser Z step gives FA ≈ 0.66.

## Implementation

- Code: `rezim_a_anisotropy.shape_anisotropy(mask, sampling, min_voxels)`; also returns `principal_std` (√λ, ascending) and `principal_axes` (eigenvectors).
- A layer FA is valid only with at least `min_voxels` voxels (`A_*_valid`).
