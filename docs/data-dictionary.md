# Data Dictionary: Exported CSV Schema

This reference document defines all 30+ columns exported in ExQt's primary measurement tables (`_Output_Batch_3d.csv` and `Result_Primary_Condensates.csv`).

---

## 1. Object Identification & Centroid Coordinates

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `filename` | - | `string` | Basename of the processed 3D TIFF acquisition file. |
| `cell_id` | - | `integer` | Identifier of the specific nucleus / cell within the field of view. |
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
| `applied_pixel_size_nm` | nm | `float` | Unscaled lateral pixel size $\Delta x$ supplied during batch configuration. |
| `applied_z_step_nm` | nm | `float` | Unscaled axial slice spacing $\Delta z$ supplied during batch configuration. |

---

## 3. Photometric Intensity & Partitioning

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mean_intensity` | ADU | `float` | Mean fluorescence intensity inside the 3D condensate mask ($I_{\text{cond}}$). |
| `max_intensity` | ADU | `float` | Maximum peak pixel intensity recorded within the condensate. |
| `integrated_density` | ADU | `float` | Sum of all voxel intensities: $\sum_{v \in \text{cond}} I(v)$. |
| `nucleoplasm_mean_intensity` | ADU | `float` | Mean intensity of the diffuse nucleoplasm ($I_{\text{nuc}}$) within the nuclear ROI. |
| `camera_offset` | ADU | `float` | Dark camera electronic offset derived from the 0.5th percentile of ROI pixels. |
| `partition_coefficient` | - | `float` | Concentration ratio: $(I_{\text{cond}} - \text{offset}) / (I_{\text{nuc}} - \text{offset})$. |
| `K_valid` | - | `boolean` | `True` if denominator $\ge 1.0$ and $K_{\text{part}} \ge 1.0$; `False` otherwise. |
| `K_invalid_reason` | - | `string` | Diagnostic text describing denominator collapse or non-condensation. |

---

## 4. Mode A: Radial Profiling & Concentration Gradients

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mean_intensity_core` | ADU | `float` | Mean intensity within the innermost 33% distance layer ($d_{\text{norm}} > 2/3$). |
| `mean_intensity_middle` | ADU | `float` | Mean intensity within the intermediate layer ($1/3 < d_{\text{norm}} \le 2/3$). |
| `mean_intensity_shell` | ADU | `float` | Mean intensity within the outer boundary layer ($0 < d_{\text{norm}} \le 1/3$). |
| `Delta_intensity_core_shell` | ADU | `float` | Radial concentration difference: $I_{\text{core}} - I_{\text{shell}}$. |
| `condensate_class` | - | `string` | Categorical classification: `Core-Enriched`, `Shell-Enriched`, or `Uniform`. |

---

## 5. 3D Fractional Anisotropy & Geometric Null Model

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `A_object` | - | `float` | 3D Fractional Anisotropy (elongation) of the entire condensate mask. |
| `A_core` | - | `float` | 3D Fractional Anisotropy of the central core shell. |
| `A_middle` | - | `float` | 3D Fractional Anisotropy of the intermediate shell. |
| `A_shell` | - | `float` | 3D Fractional Anisotropy of the outermost boundary shell. |
| `Delta_A_core_shell` | - | `float` | Measured radial FA difference: $\text{FA}_{\text{core}} - \text{FA}_{\text{shell}}$. |
| `null_delta_FA_core_shell` | - | `float` | Expected geometric gradient from a synthetic homogeneous matching ellipsoid. |
| `delta_FA_excess` | - | `float` | Biological excess anisotropy: $\Delta\text{FA}_{\text{measured}} - \Delta\text{FA}_{\text{null}}$. |

---

## 6. Quality Control & Topological Metrics

| Column Name | Physical Unit | Data Type | Mathematical Definition & Description |
| :--- | :---: | :---: | :--- |
| `mode_a_core_voxels` | voxels | `integer` | Count of voxels residing in the core layer (must be $\ge 20$ for valid FA). |
| `mode_a_z_topology_status` | - | `string` | Topological section continuity: `'pass'`, `'review'`, or `'fail'`. |
| `mode_a_z_split_slice_fraction`| - | `float` | Fraction of occupied $Z$-slices containing $\ge 2$ disconnected components. |
| `mode_a_object_touches_edge` | - | `boolean` | `True` if object bounding box intersects the 3D image border. |
| `mode_a_object_touches_roi_edge` | - | `boolean` | `True` if object voxels intersect the nuclear envelope mask boundary. |
| `primary_qc_valid` | - | `boolean` | `True` if object passes all criteria to join the primary benchmark cohort. |
| `mode_a_qc_reason` | - | `string` | Semicolon-delimited list of rejection reasons for excluded objects. |