# 3D Fractional Anisotropy (FA)

**How stretched out is a condensate?**

Fractional Anisotropy (FA) is a single number between **0** and **1** that describes the 3D shape of an object:

| FA value | Prolate ellipsoid axis ratio (a:1:1) | What the shape looks like |
| :--- | :---: | :--- |
| **0.0** | 1.0 | Sphere |
| **0.47** | 1.5 | Moderately elongated |
| **0.71** | 2.0 | Clearly elongated |
| **0.82** | 2.5 | Strongly elongated |
| **→ 1** | → ∞ | Rod / filament |

(Values for a homogeneous ellipsoid. FA is computed from the covariance eigenvalues, i.e. variances, so it is not proportional to the axis ratio and approaches 1 only slowly.)

## What FA Does and Does Not Tell You

FA describes the geometry of the **segmented mask**. It does not by itself explain why an object has that shape. The measured value combines:

- the shape of the structure itself;
- **optics** – a confocal PSF is elongated along Z, so small objects look elongated even when round (see `fa_psf_elongation_limit` in `reference_values.py`, a working value not yet validated for the ExM gel);
- **segmentation** – threshold and mask quality set where the boundary lies;
- **expansion** – local gel inhomogeneity can distort shape;
- **discretisation** – objects with few voxels have noisy eigenvalues.

FA describes the shape of the mask only; it does not identify why an object has that shape.

## How It Is Calculated

ExQt takes the physical coordinates of all mask voxels (voxel index × sampling in nm, order Z, Y, X), computes their covariance matrix and its eigenvalues $\lambda_1 \le \lambda_2 \le \lambda_3$ (variances along the object's own principal axes, not along X/Y/Z):

$$\text{FA} = \sqrt{\frac{3}{2}} \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$

If all three eigenvalues are equal, FA = 0. The physical sampling is essential: the same sphere measured in voxel units on a grid with a 4× coarser Z step gives FA ≈ 0.66.

## Implementation

- Code: `rezim_a_anisotropy.shape_anisotropy(mask, sampling, min_voxels)`; also returns `principal_std` (√λ, ascending) and `principal_axes` (eigenvectors).
- A layer FA is valid only with at least `min_voxels` voxels (`A_*_valid`).
