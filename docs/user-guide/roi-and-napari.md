# Nuclear ROI & Napari 3D Preview

Accurate spatial delineation of the cell nucleus is the fundamental prerequisite for quantitative condensate biophysics. Biomolecular condensates must be evaluated exclusively within the nucleoplasm; extracellular fluorescence or cytoplasmic background must be excluded to prevent systematic bias in partition coefficient ($K_{\text{part}}$) calculations.

---

## The Dual Role of the Nuclear ROI

1. **Denying False Positives:** Delineates the nuclear boundary so that cytoplasmic granules, membrane debris, or nonspecific antibody precipitates outside the nucleus are rejected before feature extraction.
2. **Establishing the Nucleoplasmic Baseline:** The mean fluorescence intensity of the surrounding nucleoplasm ($I_{\text{nuc}}$) serves as the reference denominator in calculating the partition coefficient:
   $$K_{\text{part}} = \frac{I_{\text{condensate}} - \text{offset}}{I_{\text{nuc}} - \text{offset}}$$
   If extracellular dark pixels (values near 0) are accidentally included in the nucleoplasmic mask, $I_{\text{nuc}} - \text{offset}$ approaches zero, causing the measured partition coefficient to artificially diverge into astronomical numbers ($10^6$–$10^8$).

---

## Interactive Workflow: Napari Layer Binding & 3D Extrusion

The workflow proceeds in three steps:
1. **Layer Dispatch:** ExQt sends the multi-channel 3D stack and initial condensate labels to Napari, applying anisotropic spatial scaling `(z_step_nm / pixel_size_nm, 1, 1)`.
2. **Interactive Contouring:** The user outlines the nucleus (or multiple nuclei) on the representative focal plane using Napari's native drawing tools.
3. **Volumetric Extrusion:** ExQt automatically sweeps the 2D contour along the entire optical axis, generating a cylindrical 3D binary volume mask spanning all $Z$-slices.

*(Napari preview screenshot can be placed here as `images/napari_roi.png`)*

---

## Production Code: ROI Handling & Audit Persistence (`Batch.py`)

Below is the production implementation from [`Batch.py:475–506`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py#L475-L506) handling the interactive handoff, $Z$-extrusion, and disk persistence:

```python
# --- EXCERPT FROM Batch.py: process_batch_pair() ---

# 1. Anisotropic Napari Layer Binding
if send_layer_func:
    # In ExM, Z-step is typically larger than XY pixel size — pass scale so the preview is physically accurate
    _roi_scale = (z_step_nm / pixel_size_nm, 1, 1) if is_3d else None
    _sig_kw: dict = {"colormap": "gray", "blending": "additive"}
    _msk_kw: dict = {"opacity": 0.6, "blending": "additive"}
    if _roi_scale is not None:
        _sig_kw["scale"] = _roi_scale
        _msk_kw["scale"] = _roi_scale
    send_layer_func({
        "type": "image", "name": f"Signal ({tif_path.name})",
        "data": img_intensity, "kwargs": _sig_kw
    })
    send_layer_func({
        "type": "labels", "name": "Mask Condensates",
        "data": img_mask_process.astype(int),
        "kwargs": _msk_kw
    })

print("Waiting for user to draw ROI in Napari...")
extruded_mask = request_roi_func(img_intensity.shape, is_3d)

# 2. 3D Volumetric Extrusion
if is_3d and extruded_mask.max() > 0:
    mask_2d = extruded_mask.max(axis=0)
    extruded_mask = np.repeat(mask_2d[np.newaxis, :, :], extruded_mask.shape[0], axis=0)
roi_mask = extruded_mask > 0

# 3. Permanent Audit Persistence on Disk
roi_save_path = mask_path.with_name(f"{tif_path.stem}_ROI.tif")
try:
    tifffile.imwrite(str(roi_save_path), extruded_mask.astype(np.uint8), compression="zlib")
    print(f"      [ROI] Saved ROI mask -> {roi_save_path.name}")
except Exception as e:
    print(f"      [ROI] Warning: could not save ROI mask: {e}")
```

### Key Technical Safeguards:
1. **Aspect Ratio Preservation (`_roi_scale`)**: Passing `scale=(z_step_nm / pixel_size_nm, 1, 1)` ensures that Napari does not render thin Z-slices as flattened disks, allowing accurate volumetric assessment of whether a condensate touches the nuclear envelope.
2. **Audit Verification (`*_ROI.tif`)**: Saving the extruded ROI mask alongside the raw data allows external auditors, supervisors, or peer reviewers to independently inspect the exact nuclear boundary used in every individual calculation.