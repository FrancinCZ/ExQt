# Partition Coefficient (K_part)

$K_{\text{part}}$ is the ratio of the mean offset-corrected intensity inside an object to that of the surrounding nucleoplasm.

---

## How It Is Calculated

$$K_{\text{part}} = \frac{\max(I_{\text{object}} - \text{offset},\ 0)}{I_{\text{nucleoplasm}} - \text{offset}}$$

- $I_{\text{object}}$: mean intensity inside the object mask.
- $I_{\text{nucleoplasm}}$: mean intensity of the nucleus ROI (`cell_id`) without any segmented objects, in the Z-range where that nucleus has objects; alignment padding is excluded.
- $\text{offset}$: detector offset in ADU from **Settings → Detector offset (ADU)** (default `0` = no subtraction). It is never estimated from the image.

### Which offset to use

| Detector mode | Offset | How to set it |
| :--- | :--- | :--- |
| Photon counting | 0 ADU | Keep 0, or press **From .lif…** to read the detector settings from a Leica .lif file and record it as the source. |
| Analog (PMT, camera) | Not available from the file | Acquire a dark frame (laser off, same detector settings) and enter its median. |

The value and its source are saved per object (`K_offset_adu`, `K_offset_source`) and in `*_metadata.json → partitioning`.

---

## When K_part Is Not Reported

- **Auto-ROI:** there is no nucleus outline, so no nucleoplasm is defined (`K_invalid_reason = "auto_roi_no_nucleoplasm"`).
- **No nucleoplasm voxels** (`"no_nucleoplasm"`) or $I_{\text{nucleoplasm}} - \text{offset} \le 0$ (`"negative_denominator"`).

Values below 1 are reported, not filtered out.

---

## Things That Change K_part

- Pixels outside the nucleus inside the ROI lower $I_{\text{nucleoplasm}}$ and raise $K_{\text{part}}$; draw the ROI along the nuclear boundary.
- Blurring mixes the intensity of small objects with their surroundings.
- Intensity is used as measured; a ratio of intensities is not necessarily a ratio of concentrations.
