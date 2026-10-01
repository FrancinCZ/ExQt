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

1. **Measure the real condensate** — principal standard deviations $\sigma_1 \le \sigma_2 \le \sigma_3$ **and** the principal directions (eigenvectors of the voxel-coordinate covariance, in physical Z/Y/X nm).
2. **Build a synthetic homogeneous ellipsoid** — semi-axes $a_i = \sqrt{5}\,\sigma_i$, oriented along the **same eigenvectors**, sampled on the **same anisotropic grid** (`sampling = (Z, Y, X)` nm), so its discretisation matches the object's.
3. **Slice both into layers** — using the same 3-layer method (Core, Middle, Shell).
4. **Compare the difference:**

$$\Delta\text{FA}_{\text{excess}} = \Delta\text{FA}_{\text{measured}} - \Delta\text{FA}_{\text{null}}$$

### What Does the Result Mean?

The layers and their FA are computed from the **mask geometry only** (`mode_a_fa_type = geometric_shape`); intensity is not used. $\Delta\text{FA}_{\text{excess}}$ therefore compares the shape of the mask's inner and outer layers with those of a homogeneous ellipsoid that has the same principal axes. It says how far the mask shape departs from an ellipsoid, not how the protein is distributed inside.

| Result | Meaning |
| :--- | :--- |
| **Near zero** | The layer shapes are what an ellipsoid of the same size and orientation gives. On synthetic ellipsoids the spread is about ±0.004 (IQR), so smaller values cannot be distinguished from zero. |
| **Positive** | The inner layers are more elongated than the ellipsoid predicts – the outer outline is less elongated or less ellipsoidal than the core region (e.g. irregular or lobed masks). |
| **Negative** | The inner layers are rounder than the ellipsoid predicts – e.g. an elongated outline around a rounder central region. |

Non-ellipsoidal outlines (bent, pear- or bean-shaped masks), segmentation of the boundary and PSF elongation along Z all change the excess, so a non-zero value needs a check of the masks before any biological reading.

## Implementation and Accuracy

Implemented in `null_model.null_model_delta_fa(principal_std_nm, sampling, min_core_voxels, principal_axes=...)`, called from `rezim_a_metrics.compute_core_shell_metrics`. Output columns: `null_FA_object/shell/middle/core`, `null_delta_FA_core_shell`, `null_valid`, `null_invalid_reason` (`non_finite_principal_std`, `non_finite_principal_axes`, `too_few_voxels`, `layer_split_failed`, `invalid_layer_fa`, `object_anisotropy_unavailable`), and `principal_std_1_nm ≤ principal_std_2_nm ≤ principal_std_3_nm` (ranked, **not** Z/Y/X). Unexpected errors are not caught.

**Accuracy.** On 200 homogeneous, randomly rotated synthetic ellipsoids (semi-axes 100–400 nm, sampling 62.5 × 14.5 × 14.5 nm, true excess = 0) the model gives a median |ΔFA_excess| of 0.0035 (IQR −0.003 to +0.004), and the null model is valid for all of them. Excess values within about ±0.005 therefore cannot be distinguished from zero.
