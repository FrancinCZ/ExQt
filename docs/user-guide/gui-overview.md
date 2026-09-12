# Graphical User Interface (GUI) & Optical Calibration

The ExQt user interface (`app.py`) provides an integrated control environment designed to eliminate manual parameter entry errors and enforce reproducible optical calibration across experimental batches.

---

## Interface Layout & Functional Panels

The main window is organized into four functional configuration sections:

1. **Data Selection & File Pairing:** Path to raw 3D `.tif` intensity stacks, directory containing binary masks (`*_Mask.tif`), and target output directory.
2. **Optical Calibration & Geometry:** Lateral $XY$ pixel size in nanometers, axial $Z$-step in nanometers, and nominal linear Expansion Factor ($\text{ExF}$).
3. **Analysis Modules & Filtering:** Toggle for Mode A (Radial FA Profiling and Geometric Null Model), nuclear Auto-ROI detection, and publication volume filter bounds ($V_{\min}, V_{\max}$).
4. **Execution & Interactive Verification:** Direct launch of the Napari 3D viewer or initiation of the automated batch processing pipeline.

*(Interface screenshot can be placed here as `images/gui_overview.png`)*

---

## Optical Calibration & Voxel Scaling

In biological microscopy, voxels are typically anisotropic (axial $Z$-spacing is larger than lateral $XY$-pixel dimensions). Furthermore, in **Expansion Microscopy (ExM)**, the physical size recorded by the microscope camera must be scaled down by the linear expansion factor ($\text{ExF}$) to recover the true pre-expansion biological scale.

### Mathematical Scaling Formulation

The effective biological sampling increments $(\Delta x_{\text{bio}}, \Delta y_{\text{bio}}, \Delta z_{\text{bio}})$ are defined as:

$$\Delta x_{\text{bio}} = \frac{\Delta x_{\text{camera}}}{\text{ExF}}, \quad \Delta y_{\text{bio}} = \frac{\Delta y_{\text{camera}}}{\text{ExF}}, \quad \Delta z_{\text{bio}} = \frac{\Delta z_{\text{camera}}}{\text{ExF}}$$

The unit volume of a single voxel in biological space ($V_{\text{voxel, bio}}$) and in physical hydrogel space ($V_{\text{voxel, gel}}$) is:

$$V_{\text{voxel, bio}} = \frac{(\Delta x_{\text{bio}})^2 \times \Delta z_{\text{bio}}}{10^9} \quad [\mu\text{m}^3]$$

$$V_{\text{voxel, gel}} = \frac{(\Delta x_{\text{camera}})^2 \times \Delta z_{\text{camera}}}{10^9} \quad [\mu\text{m}^3]$$

---

## Production Code: Calibration Consistency Check

To prevent users from inadvertently combining acquisitions taken with different objective lenses or zoom factors into a single pooled dataset, ExQt verifies calibration consistency across all TIFF files using [`calibration_policy.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/calibration_policy.py):

```python
#--- EXCERPT FROM calibration_policy.py ---

def summarize_calibrations(metadata_by_file: dict) -> dict:
    """
    Detects whether TIFF image calibrations are uniform across the batch.
    Extracts unique (pixel_size_nm, z_step_nm) tuples from all processed files.
    """
    pairs = {
        (
            round(float(values["pixel_size_nm"]), 6),
            round(float(values["z_step_nm"]), 6),
        )
        for values in metadata_by_file.values()
        if "pixel_size_nm" in values and "z_step_nm" in values
    }
    return {
        "calibration_count": len(pairs),
        "calibrations": sorted(pairs),
        "has_mismatch": len(pairs) > 1,  
    }
```

### Explanatory Breakdown:
- **`round(float(...), 6)`**: Guards against floating-point rounding jitter in TIFF metadata headers.
- **`has_mismatch`**: If more than one unique calibration tuple is detected, the GUI alerts the user before generating pooled figures, preventing fatal calibration skew.