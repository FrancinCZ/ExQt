# Geometric Null Model

Layers defined by the distance from the surface have different shapes even in a homogeneous object: in an elongated ellipsoid the inner layers are more elongated than the outer ones, so $\Delta\text{FA}_{\text{core-shell}}$ is positive without any internal structure. The null model gives this geometric expectation for each object, so the measured value can be compared with it.

---

## How It Is Calculated

1. From the object's voxel coordinates (physical Z/Y/X nm): principal standard deviations $\sigma_1 \le \sigma_2 \le \sigma_3$ and the principal directions (eigenvectors).
2. A homogeneous ellipsoid with semi-axes $a_i = \sqrt{5}\,\sigma_i$, oriented along the same eigenvectors and sampled on the same anisotropic grid as the object.
3. The ellipsoid is split into layers with the same method as the object, and its layer FA are computed.

$$\Delta\text{FA}_{\text{excess}} = \Delta\text{FA}_{\text{measured}} - \Delta\text{FA}_{\text{null}}$$

---

## What the Result Describes

Layers and their FA are computed from the **mask geometry only**; intensity is not used. $\Delta\text{FA}_{\text{excess}}$ therefore describes how the layer shapes of the mask differ from those of an ellipsoid with the same principal axes.

| Result | Description |
| :--- | :--- |
| **Near zero** | Layer shapes as expected for the ellipsoid. |
| **Positive** | Inner layers more elongated than for the ellipsoid. |
| **Negative** | Inner layers rounder than for the ellipsoid. |

The value depends on the outline of the mask, on segmentation and on the optical blur along Z. ExQt does not determine the cause of a non-zero value.

On 200 homogeneous, randomly rotated synthetic ellipsoids (semi-axes 100–400 nm, sampling 62.5 × 14.5 × 14.5 nm, true excess = 0) the median |$\Delta\text{FA}_{\text{excess}}$| is 0.0035 (IQR −0.003 to +0.004); values within about ±0.005 cannot be distinguished from zero.

---

## Implementation

`null_model.null_model_delta_fa`, called from `rezim_a_metrics.compute_core_shell_metrics`. Output columns: `null_FA_object/shell/middle/core`, `null_delta_FA_core_shell`, `delta_FA_excess`, `null_valid`, `null_invalid_reason`, and `principal_std_1_nm ≤ principal_std_2_nm ≤ principal_std_3_nm` (ranked, not Z/Y/X).
