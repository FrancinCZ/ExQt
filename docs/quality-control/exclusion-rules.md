# Quality Control Exclusion Rules

Every object rejected from the primary benchmark population is assigned an explicit, machine-readable exclusion code in `_QC_Excluded.csv`. This prevents "black box" filtering and provides complete traceability for academic audit.

---

## Catalog of Rejection Criteria

### 1. `touches_image_edge` (Field-of-View Truncation)
- **Criterion:** The object bounding box reaches $Z = 0$, $Z = Z_{\max}$, $Y = 0$, $Y = Y_{\max}$, $X = 0$, or $X = X_{\max}$.
- **Physical Rationale:** An object bisected by the camera field of view has an artificially reduced volume and distorted Fractional Anisotropy.

### 2. `touches_roi_edge` (Nuclear Envelope Crossing)
- **Criterion:** Foreground voxels intersect the boundary of the nuclear ROI mask.
- **Physical Rationale:** Nuclear speckles reside inside the nucleoplasm. Condensates overlapping the nuclear envelope risk cytoplasmic signal leakage and nucleoplasmic baseline contamination.

### 3. `z_topology_split` / `z_topology_review` (Axial Fragmentation)
- **Criterion:** When stepping through optical $Z$-slices, the object splits into multiple discrete disconnected components.
- **Physical Rationale:** While a true biological condensate may have an irregular shape, frequent topological splitting across axial slices often indicates that two distinct neighboring structures were improperly merged into a single mask during thresholding.

### 4. `core_below_min_voxels` (Insufficient Core Volume)
- **Criterion:** The innermost layer (Core: $d_{\text{norm}} > 2/3$) contains fewer than 20 voxels.
- **Physical Rationale:** Calculating a $3\times3$ covariance matrix and eigenvalue decomposition on fewer than 20 voxels is statistically underpowered and yields noisy, unreliable Fractional Anisotropy estimates.

### 5. `empty_layer_{name}` (Zero-Voxel Layer)
- **Criterion:** A concentric shell contains zero foreground voxels.
- **Physical Rationale:** Occurs when an object is extremely flat or thin, preventing the distance transform from forming intermediate layers.

---

## Production Code: Z-Split Topological Assessment (`Batch.py`)

Below is the production implementation from [`Batch.py:140–195`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py#L140-L195) demonstrating the optical section component analysis:

```python
#EXCERPT FROM Batch.py: _assess_z_split_topology()

def _assess_z_split_topology(
    mask: np.ndarray,
    min_component_voxels: int = 5,
    min_component_fraction: float = 0.05,
    pass_fraction: float = 0.15,
    review_fraction: float = 0.35,
) -> dict:
    """
    Evaluates topological continuity of a 3D object across optical Z-sections.
    Detects whether an object splits into multiple independent islands within individual slices.
    """
    occupied_slices = 0
    split_slices = 0
    raw_max_components = 0
    substantial_max_components = 0

    for z_slice in mask:
        foreground_voxels = int(np.count_nonzero(z_slice))
        if foreground_voxels == 0:
            continue
        occupied_slices += 1

        labeled_slice = label(z_slice)
        component_areas = np.bincount(labeled_slice.ravel())[1:]
        raw_max_components = max(raw_max_components, int(component_areas.size))

        # Filter out minor single-voxel perimeter spurs
        minimum_area = max(float(min_component_voxels), min_component_fraction * foreground_voxels)
        substantial_components = int(np.count_nonzero(component_areas >= minimum_area))
        substantial_max_components = max(substantial_max_components, substantial_components)
        
        if substantial_components >= 2:
            split_slices += 1

    split_fraction = split_slices / occupied_slices if occupied_slices else 0.0
    
    if split_fraction <= pass_fraction:
        status = "pass"
    elif split_fraction <= review_fraction:
        status = "review"
    else:
        status = "fail"

    return {
        "occupied_slice_count": occupied_slices,
        "split_slice_count": split_slices,
        "split_slice_fraction": float(split_fraction),
        "status": status,
    }
```