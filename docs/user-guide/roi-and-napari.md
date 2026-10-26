# Nuclear ROI & Napari

The ROI defines which part of the image belongs to each nucleus:

1. **Objects outside the ROI are not analysed.**
2. **Nucleoplasm:** the ROI without the segmented objects gives $I_{\text{nucleoplasm}}$ for $K_{\text{part}}$. Pixels outside the nucleus inside the ROI lower this value and raise $K_{\text{part}}$, so draw the ROI along the nuclear boundary. Alignment padding is excluded automatically. With Auto-ROI no nucleoplasm is defined and $K_{\text{part}}$ is not reported.

---

## Drawing the ROI

1. ExQt shows the signal stack and the mask in Napari with the true slice spacing.
2. The left panel shows the steps and the number of shapes drawn. Draw one shape per nucleus (ellipse **E**, polygon **P**, select **S**, delete **Backspace**) on one or more Z-slices. Several nuclei in one image (e.g. confocal/IF) are simply several shapes that do not touch.
3. **Preview ROI** shows the final per-nucleus ROI as a labels layer and writes on each shape the nucleus number (`cell_id`) it will get. It warns when a nucleus ROI spans many Z-slices without any mask object (for example a stack that starts below the cells or ends inside the nucleus). The nucleoplasm is already limited to the slices where that nucleus has objects, so $K_{\text{part}}$ is not affected; the warning points to a stack that may not cover the nucleus, where most objects touch the first/last slice and are excluded by QC (see the Z-edge option in *Exclusion rules*).
4. Click **Confirm ROI and Continue**. If the output folder already has `<stem>_ROI.tif` from a previous run, **Use saved ROI** reuses it instead of drawing (`roi_source = saved_roi`).
5. ExQt converts the shapes to a 3D ROI (`Batch._interpolate_or_extrude_roi`):
   - shapes on different slices that overlap in XY belong to the same nucleus and are interpolated between its slices; slices outside that range are not included;
   - a nucleus drawn on a single slice is extended through the whole stack;
   - each nucleus gets its own `cell_id`; shapes that touch on the same slice are rejected with a message so they can be fixed.
6. If nothing is drawn, ExQt asks before using the whole field of view (`roi_source = empty_fallback_fov`). Size and FA are measured, but $K_{\text{part}}$ is not computed, because there is no nucleus outline for the nucleoplasm.

---

## Implementation

- `App.preview_roi` / `App.confirm_roi` convert the Napari shapes to labels and check them (helpers in `roi_tools.py`); `Batch._interpolate_or_extrude_roi` builds the 3D ROI.
- The ROI is saved in the run's output folder as `<stem>_ROI.tif` (one label per nucleus); its SHA-256 is recorded in `*_metadata.json → files.input_files`, and `roi_source` in the CSV records how it was defined.
