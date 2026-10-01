# Nuclear ROI & Napari

The ROI defines which part of the image belongs to each nucleus:

1. **Objects outside the ROI are not analysed.**
2. **Nucleoplasm:** the ROI without the segmented objects gives $I_{\text{nucleoplasm}}$ for $K_{\text{part}}$. Pixels outside the nucleus inside the ROI lower this value and raise $K_{\text{part}}$, so draw the ROI along the nuclear boundary. Alignment padding is excluded automatically. With Auto-ROI no nucleoplasm is defined and $K_{\text{part}}$ is not reported.

---

## Drawing the ROI

1. ExQt shows the signal stack and the mask in Napari with the true slice spacing.
2. Draw one shape per nucleus on one or more Z-slices and click **Confirm ROI and Continue**.
3. ExQt converts the shapes to a 3D ROI (`Batch._interpolate_or_extrude_roi`):
   - shapes on different slices that overlap in XY belong to the same nucleus and are interpolated between its slices; slices outside that range are not included;
   - a nucleus drawn on a single slice is extended through the whole stack;
   - each nucleus gets its own `cell_id`; shapes that touch on the same slice are rejected with a message so they can be fixed.
4. If nothing is drawn, ExQt asks before using the whole field of view (`roi_source = empty_fallback_fov`).

---

## Implementation

- `App.confirm_roi` converts the Napari shapes to labels and checks them; `Batch._interpolate_or_extrude_roi` builds the 3D ROI.
- The ROI is saved in the run's output folder as `<stem>_ROI.tif` (one label per nucleus); its SHA-256 is recorded in `*_metadata.json → files.input_files`, and `roi_source` in the CSV records how it was defined.
