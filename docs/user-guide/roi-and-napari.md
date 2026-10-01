# Nuclear ROI & Napari 3D Preview

Accurate spatial delineation of the cell nucleus is the fundamental prerequisite for quantitative condensate biophysics. Biomolecular condensates must be evaluated exclusively within the nucleoplasm; extracellular fluorescence or cytoplasmic background must be excluded to prevent systematic bias in partition coefficient ($K_{\text{part}}$) calculations.

---

## The Dual Role of the Nuclear ROI

1. **Denying False Positives:** Delineates the nuclear boundary so that cytoplasmic granules, membrane debris, or nonspecific antibody precipitates outside the nucleus are rejected before feature extraction.
2. **Establishing the Nucleoplasmic Baseline:** The mean fluorescence intensity of the surrounding nucleoplasm ($I_{\text{nuc}}$) serves as the reference denominator in calculating the partition coefficient:
   $$K_{\text{part}} = \frac{I_{\text{condensate}} - \text{offset}}{I_{\text{nuc}} - \text{offset}}$$
   Pixels outside the nucleus (gel, background, alignment padding) lower $I_{\text{nuc}}$ and therefore raise $K_{\text{part}}$. On synthetic data an ROI reaching 10 px (of a 30 px nucleus radius) into the gel raised K from 5.0 to 8.7. Alignment padding is excluded automatically; the ROI outline is the user's responsibility. With Auto-ROI no nucleoplasm is defined and K is not reported.

---

## Interactive Workflow

1. **Layer dispatch:** ExQt shows the signal stack and the mask in Napari with scale `(z_step_nm / pixel_size_nm, 1, 1)`, so slices are displayed with their true spacing.
2. **Drawing:** draw one shape per nucleus on one or more Z-slices.
3. **Conversion to a 3D ROI** (`Batch._interpolate_or_extrude_roi`):
   - shapes on different slices that overlap in XY belong to the same nucleus and are interpolated between that nucleus's own slices (signed distance fields); slices outside its range stay empty;
   - a nucleus drawn on a single slice is extruded through the whole stack;
   - every nucleus gets its own `cell_id`; shapes that touch on the same slice are rejected with a message so they can be fixed.
4. **Nothing drawn:** ExQt asks before using the whole field of view (`roi_source = empty_fallback_fov`).

## Implementation

- `App.confirm_roi` converts the Napari shapes to labels and checks the grouping before the worker continues; `Batch._interpolate_or_extrude_roi` builds the 3D ROI.
- The ROI is saved in the run's output folder as `<stem>_ROI.tif` (one label per nucleus), and its SHA-256 is recorded in `*_metadata.json → files.input_files`, so every result can be traced to the exact ROI; `roi_source` in the CSV records how it was defined.
