# Partition Coefficient (K_part)

**How concentrated is the protein inside the condensate?**

The partition coefficient ($K_{\text{part}}$) tells you how many times more concentrated (brighter) the target protein is inside the droplet compared to the surrounding nucleoplasm.

---

## How It Is Calculated

$$K_{\text{part}} = \frac{I_{\text{condensate}} - \text{offset}}{I_{\text{nucleoplasm}} - \text{offset}}$$

Where:
- $I_{\text{condensate}}$: Average brightness inside the condensate.
- $I_{\text{nucleoplasm}}$: Average brightness of the diffuse nucleoplasm (inside the nucleus, but outside any condensates).
- $\text{offset}$: Dark camera baseline noise.

---

## Safety Checks

To prevent bugs and skewed statistics, ExQt applies three automatic safety checks:

1. **Dark offset inside the nucleus:** Camera baseline is measured strictly inside the nucleus, avoiding zero-padded image borders.
2. **Division-by-zero protection:** If the nucleoplasm background is almost as dark as the camera offset ($I_{\text{nucleoplasm}} - \text{offset} \le 1.0$), division would explode to infinity. ExQt sets $K_{\text{part}} = \text{NaN}$ instead.
3. **True condensation check:** A real condensate must be more concentrated than its surroundings ($K_{\text{part}} \ge 1.0$). If it is darker than the background, it is flagged as invalid (`K_valid = False`).

---

## Production Code: Partitioning Engine (`Batch.py`)

Below is the production implementation from [`Batch.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py):

```python
#EXCERPT FROM Batch.py: process_batch_pair()

#1. Estimate dark camera offset from lowest 0.5% percentile INSIDE nuclear ROI
roi_pixels = img_intensity[roi_mask] if np.any(roi_mask) else img_intensity.ravel()
camera_offset = float(np.percentile(roi_pixels, 0.5)) if roi_pixels.size > 0 else 0.0

#2. Extract nucleoplasm intensity excluding all segmented condensates
nuc_pixels = img_intensity[roi_mask & (~(labeled_mask > 0))]
nucleoplasm_mean = float(np.mean(nuc_pixels)) if nuc_pixels.size > 0 else float("nan")

#3. Denominator-safe evaluation per segmented region
denom = nucleoplasm_mean - camera_offset
if np.isnan(nucleoplasm_mean) or denom <= 1.0:
    k_val = float("nan")
    k_valid = False
    k_invalid_reason = "Nucleoplasm intensity <= camera offset + 1.0"
else:
    raw_k = (region.mean_intensity - camera_offset) / denom
    if raw_k < 1.0:
        k_val = float("nan")
        k_valid = False
        k_invalid_reason = "K < 1.0 (signal below nucleoplasm)"
    else:
        k_val = float(raw_k)
        k_valid = True
        k_invalid_reason = ""

row["camera_offset"] = camera_offset
row["nucleoplasm_mean_intensity"] = nucleoplasm_mean
row["partition_coefficient"] = k_val
row["K_valid"] = k_valid
row["K_invalid_reason"] = k_invalid_reason
```