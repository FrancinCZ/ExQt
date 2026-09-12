# 3D Fractional Anisotropy (FA)

**How stretched out is a condensate?**

Fractional Anisotropy (FA) is a single number between **0** and **1** that describes the 3D shape of an object:

| FA Value | What the shape looks like | Everyday analogy |
| :--- | :--- | :--- |
| **0.0** | Completely spherical | A round ball or bubble |
| **0.3 – 0.5** | Moderately elongated | An egg or rugby ball |
| **> 0.7** | Highly stretched | A rod, tube, or filament |

---

## Why Does Condensate Shape Change?

Inside the nucleus, condensates are not floating in empty space — they are surrounded by an elastic network of chromatin fibers:

- **Small condensates:** Surface tension easily overcomes the surrounding chromatin mesh, pulling small droplets into round spheres (low FA), much like tiny raindrops in air.
- **Larger condensates:** As condensates grow larger, their rounding force becomes weaker relative to their size. The surrounding chromatin network resists expansion, squeezing the condensate into irregular or elongated shapes (higher FA).

---

## How It Is Calculated

ExQt computes the 3D spread of all voxels along the three principal axes ($X, Y, Z$) using a covariance matrix:

$$\text{FA} = \sqrt{\frac{3}{2}} \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$

Where $\lambda_1, \lambda_2, \lambda_3$ represent the spread along the longest, middle, and shortest axis, and $\bar{\lambda}$ is their average. If all three axes are equal, the numerator becomes zero and $\text{FA} = 0$.

---

## Production Code: 3D FA Engine (`rezim_a_anisotropy.py`)

Below is the production implementation from [`rezim_a_anisotropy.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/rezim_a_anisotropy.py):

```python
#EXCERPT FROM rezim_a_anisotropy.py: shape_anisotropy()

def shape_anisotropy(mask: np.ndarray, sampling: tuple | None = None, min_voxels: int = 20) -> dict:
    binary_mask = np.asarray(mask, dtype=bool)
    coordinates = np.argwhere(binary_mask)
    voxel_count = int(coordinates.shape[0])
    if voxel_count == 0:
        raise ValueError("Mask contains no foreground voxels.")

    #1. Scale voxel indices by physical pixel dimensions
    scale = np.ones(binary_mask.ndim, dtype=float) if sampling is None else np.asarray(sampling, dtype=float)
    physical_coordinates = coordinates * scale
    centered = physical_coordinates - physical_coordinates.mean(axis=0)

    #2. Covariance matrix & eigenvalue decomposition
    covariance = centered.T @ centered / voxel_count
    eigenvalues = np.linalg.eigvalsh(covariance)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    eigenvalue_sum_of_squares = float(np.sum(eigenvalues ** 2))

    if eigenvalue_sum_of_squares == 0.0:
        fractional_anisotropy = np.nan
        principal_std_nm = np.full(binary_mask.ndim, np.nan)
    else:
        eigenvalue_mean = float(np.mean(eigenvalues))
        fractional_anisotropy = float(
            np.sqrt(
                binary_mask.ndim / (binary_mask.ndim - 1)
                * np.sum((eigenvalues - eigenvalue_mean) ** 2)
                / eigenvalue_sum_of_squares
            )
        )
        principal_std_nm = np.sqrt(eigenvalues)

    return {
        "fractional_anisotropy": fractional_anisotropy,
        "principal_std": principal_std_nm,
        "voxel_count": voxel_count,
        "anisotropy_valid": (voxel_count >= min_voxels and np.isfinite(fractional_anisotropy)),
    }
```