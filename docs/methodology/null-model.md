# Geometric Null Model: Separating Shape from Biology

When you peel layers off an elongated object, a curious geometric effect happens: **the inside naturally looks more stretched out than the outside**, even if the object is completely uniform.

ExQt includes a **Geometric Null Model** (`null_model.py`) to prevent mistaking this natural geometric effect for real biological structure.

---

## The Problem: The "Cucumber Effect"

Imagine an elongated cucumber or a rugby ball. If you peel layers off the surface toward the center:
1. The innermost core that remains is naturally skinnier and more elongated than the outer shell.
2. Even in a completely solid, uniform object with no internal structure, the core will have a higher shape anisotropy than the surface ($\text{FA}_{\text{core}} > \text{FA}_{\text{shell}}$).
3. If you didn't account for this, you might falsely claim that every elongated condensate has a special "internal scaffold".

---

## The Solution: A Synthetic Reference

To solve this, ExQt creates a computer-generated "twin" for every condensate:

1. **Measure the real condensate** — get its 3D length, width, and height.
2. **Build a synthetic smooth ellipsoid** — with the exact same dimensions, but completely uniform inside.
3. **Slice both into layers** — using the same 3-layer method (Core, Middle, Shell).
4. **Compare the difference:**

$$\Delta\text{FA}_{\text{excess}} = \Delta\text{FA}_{\text{measured}} - \Delta\text{FA}_{\text{null}}$$

### What Does the Result Mean?

| Result | Plain-language meaning | Biological interpretation |
| :--- | :--- | :--- |
| **Positive** ($\Delta\text{FA}_{\text{excess}} > 0$) | The core is more stretched than geometry predicts. | **Real biological structure:** e.g., an internal protein scaffold or fibrillar assembly. |
| **Near Zero** ($\Delta\text{FA}_{\text{excess}} \approx 0$) | The core elongation matches the fake ellipsoid. | **Pure geometry:** The elongated core is just a natural consequence of the outer shape. |
| **Negative** ($\Delta\text{FA}_{\text{excess}} < 0$) | The core is rounder than the outer shape. | **Fluid core:** A round liquid droplet enclosed inside an elongated interface. |

---

## Production Code: Synthetic Ellipsoid Null Generator (`null_model.py`)

Below is the production implementation from [`null_model.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/null_model.py):

```python
#EXCERPT FROM null_model.py

def _make_ellipsoid(semi_axes_voxels: np.ndarray) -> np.ndarray:
    """Constructs a binary 3D ellipsoid mask with given semi-axes (in voxels)."""
    radii = np.asarray(semi_axes_voxels, dtype=float)
    half = np.ceil(radii).astype(int) + 2  #2-voxel padding
    z, y, x = np.mgrid[
        -half[0]:half[0] + 1,
        -half[1]:half[1] + 1,
        -half[2]:half[2] + 1
    ]
    inside = (z / radii[0]) ** 2 + (y / radii[1]) ** 2 + (x / radii[2]) ** 2 <= 1.0
    return inside


def null_model_delta_fa(principal_std_nm: np.ndarray, sampling: tuple, min_core_voxels: int = 20) -> dict:
    """
    Constructs a homogeneous ellipsoid matching the real object's principal standard deviations
    and computes its expected geometric radial FA gradient.
    """
    principal_std_nm = np.asarray(principal_std_nm, dtype=float)
    sampling = np.asarray(sampling, dtype=float)

    if not np.all(np.isfinite(principal_std_nm)) or np.any(principal_std_nm <= 0):
        return {"null_valid": False}

    #For a uniform solid ellipsoid with semi-axis R, standard deviation sigma = R / sqrt(5) therefore, semi-axis in voxels = sqrt(5) * sigma_nm / voxel_sampling_nm
    semi_axes_voxels = np.sqrt(5.0) * principal_std_nm / sampling
    if np.any(semi_axes_voxels < 1.0):
        return {"null_valid": False}

    #1. Synthesize binary ellipsoid
    synthetic_mask = _make_ellipsoid(semi_axes_voxels)

    #2. Partition using identical EDT-thirds algorithm
    layers = split_core_middle_shell(synthetic_mask, sampling=sampling, min_core_voxels=min_core_voxels)
    if not layers["qc"]["core_valid"]:
        return {"null_valid": False}

    #3. Compute FA across synthetic layers
    fa_shell = shape_anisotropy(layers["shell"], sampling=sampling)["fractional_anisotropy"]
    fa_middle = shape_anisotropy(layers["middle"], sampling=sampling)["fractional_anisotropy"]
    fa_core = shape_anisotropy(layers["core"], sampling=sampling)["fractional_anisotropy"]

    return {
        "null_FA_shell": fa_shell,
        "null_FA_middle": fa_middle,
        "null_FA_core": fa_core,
        "null_delta_FA_core_shell": fa_core - fa_shell,
        "null_valid": True,
    }
```