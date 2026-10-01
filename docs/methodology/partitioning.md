# Partition Coefficient (K_part)

**How concentrated is the protein inside the condensate?**

The partition coefficient ($K_{\text{part}}$) tells you how many times more concentrated (brighter) the target protein is inside the droplet compared to the surrounding nucleoplasm.

---

## How It Is Calculated

$$K_{\text{part}} = \frac{I_{\text{condensate}} - \text{offset}}{I_{\text{nucleoplasm}} - \text{offset}}$$

Where:
- $I_{\text{condensate}}$: Average brightness inside the condensate.
- $I_{\text{nucleoplasm}}$: Average brightness of the diffuse nucleoplasm (inside the nucleus, but outside any condensates).
- $\text{offset}$: Detector/camera offset in ADU, entered in **Settings → Detector offset (ADU)** (default `0` = no subtraction).

### Which offset to use

| Detector mode | Offset | How to set it |
| :--- | :--- | :--- |
| Photon counting (e.g. Leica HyD S/X/R in counting mode) | 0 ADU – photons are counted without an electronic offset | Keep 0, or press **From .lif…** to read the detector settings and record the .lif file as the source. |
| Analog (PMT, analog HyD, camera) | Unknown from the file – Leica stores "Offset" there as an instrument setting, not in ADU | Acquire a dark frame (laser off, same detector settings) and enter its median. |

The source of the value is saved per object in `K_offset_source` (`manual` or `lif_photon_counting`) and in `*_metadata.json → partitioning.offset_source` (file name, SHA-256 and detector settings).

---

## Rules Implemented in `Batch.process_condensates`

1. **Explicit offset, never estimated from the image.** A low percentile of ROI pixels is not a camera offset: it follows the darkest thing inside the ROI and the noise level, so K would depend on ROI shape and noise, and a tight ROI inside the nucleus could inflate K several-fold. The applied value is stored per object in `K_offset_adu` with `K_offset_method = "explicit_setting"`, and in `*_metadata.json → partitioning`.
2. **Nucleoplasm per nucleus.** $I_{\text{nucleoplasm}}$ is the mean of the nucleus ROI (`cell_id`) minus all segmented objects, restricted to the Z-range in which that nucleus has objects.
3. **Alignment padding excluded.** Aligned stacks come with `<stem>_alignment_valid.tif` (written by `Tools → Align Z-stacks`). Voxels marked 0 (zero padding, shifted-in borders) never enter the nucleoplasm mean. Stacks aligned with an older ExQt have no such mask; a warning is printed and they should be re-aligned.
4. **No K for Auto-ROI.** Without a nucleus outline the "nucleoplasm" would include gel and background, so `partition_coefficient = NaN`, `K_valid = False`, `K_invalid_reason = "auto_roi_no_nucleoplasm"`.
5. **Non-positive denominator.** If $I_{\text{nucleoplasm}} - \text{offset} \le 0$, K is NaN with `K_invalid_reason = "negative_denominator"`; with no nucleoplasm voxels, `"no_nucleoplasm"`. K < 1 is **not** filtered out.

The numerator is clipped at 0 (`max(I_condensate - offset, 0)`).
