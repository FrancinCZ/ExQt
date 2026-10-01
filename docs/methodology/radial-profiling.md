# Radial Profiling (Mode A: Core–Shell Architecture)

Mode A (Radial FA Profiling) divides each 3D object into three concentric layers – **Core**, **Middle** and **Shell** – and reports the shape (FA) and mean intensity of each layer, so the inner and outer parts of an object can be compared.

---

## The Three Layers

Every voxel gets its distance to the object surface from an anisotropic Euclidean distance transform (`sampling = (Z, Y, X)` nm), normalised by the maximum distance $d_{\max}$:

| Layer | Normalised distance $d / d_{\max}$ | Share of a sphere's volume |
| :--- | :--- | :---: |
| **Shell** | $0 < d \le 1/3$ | ≈ 70 % |
| **Middle** | $1/3 < d \le 2/3$ | ≈ 26 % |
| **Core** | $d > 2/3$ | ≈ 4 % |

Because the core is only about 4 % of the volume, an object needs roughly 500 voxels or more for the core to reach `min_core_voxels` (default 20). This is why the primary set over-represents large objects (see *Known Limitation* in the QC pipeline page).

## Intensity Difference Between Layers

For each layer ExQt reports the mean intensity ($I_{\text{core}}$, $I_{\text{middle}}$, $I_{\text{shell}}$) and the difference $\Delta I_{\text{core-shell}} = I_{\text{core}} - I_{\text{shell}}$.

The class (`condensate_class`) uses a two-sided Welch t-test of core vs shell voxel intensities at α = 0.01 (`gradient_p_core_shell`, `gradient_significant`):

- significant and $\Delta I > 0$ → *Core-Enriched*; significant and $\Delta I < 0$ → *Shell-Enriched*;
- not significant → *No Significant Gradient*; objects that failed QC → *Unclassified (QC)*.

How to read it: a significant difference only says that the centre and the rim of the mask differ in brightness. Blurring by the PSF makes the centre of even a perfectly homogeneous object brighter than its rim, so *Core-Enriched* is expected without any internal structure. Antibody penetration, self-quenching or segmentation that includes blurred edge voxels can also shift the difference. The test assumes voxel noise is independent (true for photon noise, not after deconvolution or smoothing).

## Implementation

- Layers: `rezim_a_core_shell.split_core_middle_shell(mask, min_core_voxels, sampling)` (scheme `edt_over_max_thirds_v1`).
- Metrics and class: `rezim_a_metrics.compute_core_shell_metrics`, `rezim_a_metrics.classify_particle`; the FA class boundary comes from `reference_values.py`.
