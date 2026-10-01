# Data Dictionary: Exported CSV Schema

This reference document defines all 30+ columns exported in ExQt's primary measurement tables (`_Output_Batch_3d.csv` and `<run>_Primary_Condensates.csv`).

---

## 1. Object Identification & Centroid Coordinates

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `filename` | - | `string` | Basename of the processed 3D TIFF acquisition file. |
| `cell_id` | - | `integer` | Identifier of the specific nucleus / cell within the field of view. |
| `roi_source` | - | `string` | How the analysed region was defined: `auto_fov` (Auto-ROI, whole field of view), `manual` (ROI drawn by the user), `empty_fallback_fov` (nothing drawn, user explicitly confirmed the whole field of view). |
| `object_id` | - | `integer` | Unique connected-component label assigned to the condensate. |
| `Z_px` | pixels | `float` | Axial centroid coordinate within the raw optical stack. |
| `Y_px` | pixels | `float` | Vertical lateral centroid coordinate. |
| `X_px` | pixels | `float` | Horizontal lateral centroid coordinate. |

---

## 2. Volumetric & Geometric Scaling

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `volume_px` | voxels | `integer` | Total number of foreground voxels ($N_{\text{vox}}$) comprising the 3D mask. |
| `volume_gel_um3` | $\mu\text{m}^3$ | `float` | Physical volume in hydrogel: $N_{\text{vox}} \times \frac{\Delta x \Delta y \Delta z}{10^9}$. |
| `volume_bio_um3` | $\mu\text{m}^3$ | `float` | Pre-expansion biological volume: $\text{Vol}_{\text{gel}} / \text{ExF}^3$. |
| `equivalent_diameter_um`| $\mu\text{m}$ | `float` | Diameter of an equivalent sphere: $(6 \cdot \text{Vol}_{\text{bio}} / \pi)^{1/3}$. |
| `applied_pixel_size_nm` | nm | `float` | Unscaled lateral pixel size $\Delta x$ actually applied to this file (Settings or, if chosen, the file's TIFF calibration; source in `*_metadata.json` → `applied_calibration_by_file`). |
| `applied_z_step_nm` | nm | `float` | Unscaled axial slice spacing $\Delta z$ actually applied to this file (same source rules as `applied_pixel_size_nm`). |

---

## 3. Photometric Intensity & Partitioning

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mean_intensity` | ADU | `float` | Mean fluorescence intensity inside the 3D condensate mask ($I_{\text{cond}}$). |
| `max_intensity` | ADU | `float` | Maximum peak pixel intensity recorded within the condensate. |
| `integrated_density` | ADU | `float` | Sum of all voxel intensities: $\sum_{v \in \text{cond}} I(v)$. |
| `nucleoplasm_mean_intensity` | ADU | `float` | Mean intensity of the diffuse nucleoplasm ($I_{\text{nuc}}$) within this nucleus ROI (objects and alignment padding excluded, object Z-range only). NaN in Auto-ROI. |
| `K_offset_adu` | ADU | `float` | Detector offset subtracted for $K_{\text{part}}$; explicit setting from Settings (0 = none). |
| `K_offset_method` | - | `string` | How the offset was obtained; currently always `explicit_setting`. |
| `K_offset_source` | - | `string` | Where the explicit value came from: `manual` or `lif_photon_counting` (read from a Leica .lif file, details in `*_metadata.json`). |
| `partition_coefficient` | - | `float` | Concentration ratio: $\max(I_{\text{cond}} - \text{offset}, 0) / (I_{\text{nuc}} - \text{offset})$. NaN when `K_valid` is `False`. |
| `K_valid` | - | `boolean` | `True` if a manual ROI nucleoplasm exists and $I_{\text{nuc}} - \text{offset} > 0$. K < 1 is not filtered. |
| `K_invalid_reason` | - | `string` | `auto_roi_no_nucleoplasm`, `no_nucleoplasm`, or `negative_denominator`. |

---

## 4. Mode A: Radial Profiling & Concentration Gradients

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mean_intensity_core` | ADU | `float` | Mean intensity within the innermost 33% distance layer ($d_{\text{norm}} > 2/3$). |
| `mean_intensity_middle` | ADU | `float` | Mean intensity within the intermediate layer ($1/3 < d_{\text{norm}} \le 2/3$). |
| `mean_intensity_shell` | ADU | `float` | Mean intensity within the outer boundary layer ($0 < d_{\text{norm}} \le 1/3$). |
| `Delta_intensity_core_shell` | ADU | `float` | Radial concentration difference: $I_{\text{core}} - I_{\text{shell}}$. |
| `condensate_class` | - | `string` | `<shape> / <gradient>`: shape `Globular (Low FA)` (A_object < `fa_globular_boundary`, 0.65, working value in `reference_values.py`) or `Elongated (High FA)`; gradient `Core-Enriched` / `Shell-Enriched` only if `gradient_significant`, else `No Significant Gradient`, or `Profile-Unavailable`. Also `Punctate / Small Cluster` (< 3 × min core voxels), `Unclassified` (no FA) and `Unclassified (QC)` (object failed Mode A QC). **Core-Enriched is expected from PSF blurring of a homogeneous object** – not evidence of biological enrichment without a null/PSF comparison. |
| `gradient_t_core_shell` | - | `float` | Welch t statistic, core vs shell voxel intensities. |
| `gradient_p_core_shell` | - | `float` | Two-sided Welch p-value (NaN if a layer has < 2 voxels). |
| `gradient_significant` | - | `boolean` | `gradient_p_core_shell < gradient_alpha`. |
| `gradient_test` / `gradient_alpha` | - | `string` / `float` | `welch_t_two_sided`, 0.01. |

---

## 5. 3D Fractional Anisotropy & Geometric Null Model

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `A_object` | - | `float` | 3D Fractional Anisotropy (elongation) of the entire condensate mask. |
| `A_core` | - | `float` | 3D Fractional Anisotropy of the core layer. |
| `A_middle` | - | `float` | 3D Fractional Anisotropy of the intermediate shell. |
| `A_shell` | - | `float` | 3D Fractional Anisotropy of the outermost boundary shell. |
| `Delta_A_core_shell` | - | `float` | Measured radial FA difference: $\text{FA}_{\text{core}} - \text{FA}_{\text{shell}}$. |
| `principal_std_1_nm`, `principal_std_2_nm`, `principal_std_3_nm` | nm | `float` | Principal standard deviations of the object's voxel coordinates (physical, after ExF), **ascending**; ranked, not Z/Y/X. |
| `null_delta_FA_core_shell` | - | `float` | Expected geometric gradient of a homogeneous ellipsoid with the object's principal stds **and eigenvectors**, sampled on the same anisotropic grid. |
| `null_valid` | - | `boolean` | `True` if the null shell/middle/core FA are all valid. |
| `null_invalid_reason` | - | `string` | Why the null model is unavailable (empty if valid): `non_finite_principal_std`, `non_finite_principal_axes`, `too_few_voxels`, `layer_split_failed`, `invalid_layer_fa`, `object_anisotropy_unavailable`. |
| `delta_FA_excess` | - | `float` | $\Delta\text{FA}_{\text{measured}} - \Delta\text{FA}_{\text{null}}$: how far the mask's layer shapes depart from a same-shape homogeneous ellipsoid (mask geometry only, not protein distribution). |

---

## 6. Quality Control & Topological Metrics

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mode_a_core_voxels` | voxels | `integer` | Count of voxels residing in the core layer (must be $\ge$ `min_core_voxels`, default 20, for valid FA). |
| `mode_a_z_topology_status` | - | `string` | Topological section continuity: `'pass'`, `'review'`, or `'fail'`. |
| `mode_a_z_split_slice_fraction`| - | `float` | Fraction of occupied $Z$-slices containing $\ge 2$ disconnected components. |
| `touches_image_edge` | - | `boolean` | `True` if the object reaches the image/stack border or alignment padding (valid-voxel mask). All modes; excluded from primary. |
| `touches_roi_edge` | - | `boolean` | `True` if the ROI boundary cuts the object. All modes; excluded from primary. |
| `mode_a_object_touches_edge` | - | `boolean` | Mode A copy of `touches_image_edge`. |
| `mode_a_object_touches_roi_edge` | - | `boolean` | Mode A copy of `touches_roi_edge`. |
| `alignment_status` | - | `string` | Aligner overall status of the input stack: `PASS`, `REVIEW` (analysed, never primary) or `not_aligned` (no `<stem>_drift.csv`). `FAIL` stacks are not analysed. |
| `primary_qc_valid` | - | `boolean` | `True` if object passes all criteria to join the primary benchmark cohort. |
| `mode_a_qc_reason` | - | `string` | Semicolon-delimited list of rejection reasons for excluded objects. |
