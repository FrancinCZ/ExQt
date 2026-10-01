# Radial Profiling (Mode A)

Mode A divides each 3D object into three concentric layers – **Core**, **Middle** and **Shell** – and reports the FA and the mean intensity of each layer.

---

## The Three Layers

Every voxel gets its distance to the object surface from an anisotropic Euclidean distance transform (`sampling = (Z, Y, X)` nm), normalised by the maximum distance $d_{\max}$:

| Layer | Normalised distance $d / d_{\max}$ | Share of a sphere's volume |
| :--- | :--- | :---: |
| **Shell** | $0 < d \le 1/3$ | ≈ 70 % |
| **Middle** | $1/3 < d \le 2/3$ | ≈ 26 % |
| **Core** | $d > 2/3$ | ≈ 4 % |

Because the core holds only about 4 % of the volume, an object needs roughly 500 voxels or more for the core to reach `min_core_voxels` (default 20); the primary set therefore contains mainly larger objects.

---

## Intensity Difference Between Layers

ExQt reports the mean intensity of each layer and $\Delta I_{\text{core-shell}} = I_{\text{core}} - I_{\text{shell}}$. The class (`condensate_class`) uses a two-sided Welch t-test of core vs shell voxel intensities at α = 0.01:

- significant and $\Delta I > 0$ → *Core-Enriched*; significant and $\Delta I < 0$ → *Shell-Enriched*;
- not significant → *No Significant Gradient*; objects that failed QC → *Unclassified (QC)*.

A significant result means only that the central and outer voxels of the mask differ in mean intensity; it does not identify the cause. Optical blurring alone makes the centre of a homogeneous object brighter than its rim, so *Core-Enriched* also occurs without internal structure. The test assumes independent noise in each voxel.

---

## Implementation

- Layers: `rezim_a_core_shell.split_core_middle_shell` (scheme `edt_over_max_thirds_v1`).
- Metrics and class: `rezim_a_metrics.compute_core_shell_metrics`, `rezim_a_metrics.classify_particle`; the FA class boundary is a working value in `reference_values.py`.
