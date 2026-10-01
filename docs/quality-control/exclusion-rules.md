# Quality Control Exclusion Rules

Every size-eligible object that is not primary has a machine-readable reason in `<run>_QC_Excluded.csv`.

---

## Rejection Criteria

### 1. `touches_image_edge` (Image Edge)
- **Criterion:** the object reaches the first or last position in Z, Y or X, or (aligned stacks) touches the alignment padding marked 0 in `<stem>_alignment_valid.tif`.
- **Reason:** the object may continue outside the image, so its volume and shape are incomplete.

### 2. `touches_roi_edge` (ROI Edge)
- **Criterion:** object voxels are adjacent to voxels outside the ROI (never set when the ROI is the whole field of view).
- **Reason:** the ROI cut the object, so its volume, shape and intensity are incomplete.

### 3. `z_topology_review` / `z_topology_fail` (Split Along Z)
- **Criterion:** on a substantial fraction of its Z-slices the object consists of several separate parts.
- **Reason:** the mask may join separate structures, so its shape metrics would not describe a single object.

### 4. `core_below_min_voxels` (Small Core)
- **Criterion:** the core layer has fewer than `min_core_voxels` voxels (default 20).
- **Reason:** FA from very few voxels is noisy.

### 5. `empty_layer_{name}` (Empty Layer)
- **Criterion:** a layer contains no voxels.
- **Reason:** occurs for very thin or flat objects; layer metrics cannot be computed.

---

## Implementation

- Flags: `Batch.process_condensates`; Z-topology: `Batch._assess_z_split_topology` (thresholds recorded in `*_metadata.json → mode_a`).
- Primary set: `postprocessing._prepare_reporting_frames` (see the QC pipeline page).
