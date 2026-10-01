# Quality Control Exclusion Rules

Every object rejected from the primary benchmark population is assigned an explicit, machine-readable exclusion code in `_QC_Excluded.csv`. This prevents "black box" filtering and provides complete traceability for academic audit.

---

## Catalog of Rejection Criteria

### 1. `touches_image_edge` (Field-of-View Truncation)
- **Criterion:** The object bounding box reaches $Z = 0$, $Z = Z_{\max}$, $Y = 0$, $Y = Y_{\max}$, $X = 0$, or $X = X_{\max}$, **or** (aligned stacks) the object touches alignment padding marked 0 in `<stem>_alignment_valid.tif`, i.e. the original field-of-view border moved inside the array by `expand_canvas`.
- **Scope:** computed for every object in every mode (column `touches_image_edge`); Mode A keeps the copy `mode_a_object_touches_edge`.
- **Physical Rationale:** An object bisected by the camera field of view has an artificially reduced volume and distorted Fractional Anisotropy.

### 2. `touches_roi_edge` (Nuclear Envelope Crossing)
- **Criterion:** Object voxels are adjacent to voxels outside the ROI, i.e. the ROI cut the object (column `touches_roi_edge`, all modes; never set when the ROI is the whole field of view).
- **Physical Rationale:** Nuclear speckles reside inside the nucleoplasm. Condensates overlapping the nuclear envelope risk cytoplasmic signal leakage and nucleoplasmic baseline contamination.

### 3. `z_topology_split` / `z_topology_review` (Axial Fragmentation)
- **Criterion:** When stepping through optical $Z$-slices, the object splits into multiple discrete disconnected components.
- **Physical Rationale:** While a true biological condensate may have an irregular shape, frequent topological splitting across axial slices often indicates that two distinct neighboring structures were improperly merged into a single mask during thresholding.

### 4. `core_below_min_voxels` (Insufficient Core Volume)
- **Criterion:** The innermost layer (Core: $d_{\text{norm}} > 2/3$) contains fewer than `min_core_voxels` voxels (default 20).
- **Physical Rationale:** Calculating a $3\times3$ covariance matrix and eigenvalue decomposition on fewer than 20 voxels is statistically underpowered and yields noisy, unreliable Fractional Anisotropy estimates.

### 5. `empty_layer_{name}` (Zero-Voxel Layer)
- **Criterion:** A concentric shell contains zero foreground voxels.
- **Physical Rationale:** Occurs when an object is extremely flat or thin, preventing the distance transform from forming intermediate layers.

---

## Implementation

- Flags are computed in `Batch.process_condensates`; Z-topology in `Batch._assess_z_split_topology` (policy `substantial_component_fraction_v1`, thresholds recorded in `*_metadata.json → mode_a`).
- The primary set is decided in `postprocessing._prepare_reporting_frames` (see the QC pipeline page).
