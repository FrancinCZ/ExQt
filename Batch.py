import math
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import laplace
from skimage.measure import regionprops, label
from rezim_a_metrics import compute_core_shell_metrics


def _touches_image_edge(region, image_shape):
    #Return True when a region bbox reaches any image or stack boundary.
    bbox_min = region.bbox[:region.image.ndim]
    bbox_max = region.bbox[region.image.ndim:]
    return any(start == 0 or end == size for start, end, size in zip(bbox_min, bbox_max, image_shape))


def _max_components_per_z_slice(object_mask):
    #Detect split silhouettes without deciding biological object identity.

    if object_mask.ndim != 3:
        return 1
    return max(int(label(object_mask[z]).max()) for z in range(object_mask.shape[0]))

def _select_focus_slice(volume):
    #Choose the sharpest Z slice for single-slice processing
    # Laplacian variance favors slices with the strongest visible structure
    scores = [laplace(volume[z].astype(np.float32)).var() for z in range(volume.shape[0])]
    return int(np.argmax(scores))

def _resolve_channel_axis(img_raw, expected_axis=1, max_channels=6):
    #Validate or infer the channel axis used to extract signal images.
    if img_raw.ndim != 4:
        return expected_axis

    sizes = img_raw.shape
    if sizes[expected_axis] <= max_channels:
        return expected_axis

    candidates = [ax for ax, s in enumerate(sizes) if s <= max_channels]
    if len(candidates) == 1:
        found_axis = candidates[0]
        print(f"      WARNING: axis {expected_axis} has size {sizes[expected_axis]}, which doesn't "
            f"look like a channel axis. Found a plausible channel axis at position {found_axis} "
            f"(size {sizes[found_axis]}) instead - auto-adjusting.")
        return found_axis

    raise ValueError(f"Could not confidently identify the channel axis for TIF shape {sizes}.")

def get_metadata_from_tif(tif_path):
    #Read physical pixel spacing from ImageJ metadata and TIFF tags.

    with tifffile.TiffFile(tif_path) as tif:
        try:
            imagej_metadata = tif.imagej_metadata
            meta = {}
            if imagej_metadata:
                if 'spacing' in imagej_metadata:
                    meta['z_step'] = float(imagej_metadata['spacing']) * 1000 
            
            tags = tif.pages[0].tags
            if 'XResolution' in tags:
                res = tags['XResolution'].value
                if res[0] > 0 and res[1] > 0:
                    pixel_size_um = res[1] / res[0]
                    meta['pixel_size'] = pixel_size_um * 1000 
            
            return meta if meta else None
        except Exception:
            return None


def process_condensates(
    tif_path, mask_path, mode="3d", target_z_slice=None, expansion_factor=1.0,
    min_voxels=5, show_napari=True, pixel_size_nm=None, z_step_nm=None,
    signal_channel=1, dapi_channel=0, channel_axis=1, auto_roi=True,
    send_layer_func=None, request_roi_func=None, mode_a_enabled=False,
    mode_a_min_core_voxels=20, mode_a_exclude_split_slices=True
):
    #Process one raw/mask TIFF pair and return one row per valid object.

    tif_path = Path(tif_path)
    mask_path = Path(mask_path)
    print(f"\n[1/5] Loading Pair: {tif_path.name} & {mask_path.name}")

    if send_layer_func:
        send_layer_func({"type": "clear_layers"})

    #Load both inputs together so all later measurements refer to the same
    #field of view and can be displayed through the GUI callback.
    img_raw = tifffile.imread(tif_path)
    img_mask = tifffile.imread(mask_path)
    img_mask = np.squeeze(img_mask)
    img_raw = np.squeeze(img_raw)

    meta = get_metadata_from_tif(tif_path)
    if meta:
        if pixel_size_nm is None:
            pixel_size_nm = meta.get('pixel_size')
        if z_step_nm is None:
            z_step_nm = meta.get('z_step')
        print(f"      [Info] TIF metadata available: XY={meta.get('pixel_size')}nm, Z={meta.get('z_step')}nm")
    else:
        print(f"      [Info] Metadata not found, using default values from GUI.")

    if pixel_size_nm is None or z_step_nm is None:
        raise ValueError("Pixel size XY and Z-step must be provided by ExQt Settings or TIF metadata.")

    #Channel selection is isolated here - downstream processing works with a
    #single intensity volume regardless of the original TIFF layout.
    if img_raw.ndim == 4:
        ch_axis = _resolve_channel_axis(img_raw, expected_axis=channel_axis)
        img_dapi = np.take(img_raw, dapi_channel, axis=ch_axis)
        img_intensity = np.take(img_raw, signal_channel, axis=ch_axis)
    else:
        img_dapi = img_raw
        img_intensity = img_raw

    if img_mask.shape != img_intensity.shape:
        if img_mask.shape == img_intensity.shape[::-1]: 
            img_mask = np.transpose(img_mask)
            print("Mask was automatically transposed.")
        else:
            raise ValueError(f"CRITICAL ERROR: TIF shape {img_intensity.shape} does not match "
                            f"Mask shape {img_mask.shape}. Ensure your segmentation software outputs correct dimensions.")

    is_stack = img_intensity.ndim == 3
    print(f"      Mode: {mode}")

    #Convert the requested mode into one processing image/mask pair. Rezim A
    #is intentionally available only when this conversion preserves 3D data.
    if mode == "single_slice":
        if not is_stack:
            raise ValueError("single_slice mode needs a 3D (Z,Y,X) stack.")
        z_idx = (target_z_slice if target_z_slice is not None else _select_focus_slice(img_intensity))
        print(f"      Using Z-slice {z_idx} (auto-focus)")
        img_intensity = img_intensity[z_idx]
        img_dapi_process = img_dapi[z_idx]
        img_mask_process = img_mask[z_idx]
        is_3d = False

    elif mode == "2d":
        img_intensity = img_intensity if not is_stack else img_intensity.max(axis=0)
        img_dapi_process = img_dapi if not is_stack else img_dapi.max(axis=0)
        img_mask_process = img_mask if not is_stack else img_mask.max(axis=0)
        is_3d = False

    elif mode == "3d":
        if not is_stack:
            raise ValueError("3d mode needs a 3D (Z,Y,X) stack.")
        img_intensity = img_intensity
        img_dapi_process = img_dapi
        img_mask_process = img_mask
        is_3d = True

    print(f"\n[2/5] ROI extraction (Auto-ROI: {auto_roi})")
    #A manual ROI is requested through callbacks so the worker can pause
    #without touching Qt widgets from its background thread.
    if auto_roi or request_roi_func is None:
        roi_mask = np.ones_like(img_mask_process, dtype=bool)
        extruded_mask = np.ones_like(img_mask_process, dtype=int)
    else:
        if send_layer_func:
            send_layer_func({
                "type": "image", "name": f"Signal ({tif_path.name})",
                "data": img_intensity, "kwargs": {"colormap": "gray", "blending": "additive"}
            })
            send_layer_func({
                "type": "labels", "name": "Mask Condensates",
                "data": img_mask_process.astype(int),
                "kwargs": {"opacity": 0.6, "blending": "additive"}
            })
        print("Waiting for user to draw ROI in Napari...")
        extruded_mask = request_roi_func(img_intensity.shape, is_3d)
        if is_3d and extruded_mask.max() > 0:
            mask_2d = extruded_mask.max(axis=0) 
            extruded_mask = np.repeat(mask_2d[np.newaxis, :, :], extruded_mask.shape[0], axis=0)
        roi_mask = extruded_mask > 0

    print("\n[3/5] Applying provided Mask...")

    img_mask_process = img_mask_process * roi_mask

    #Binary masks need connected-component labeling pre-labeled masks keep
    #their object IDs so cell/object associations survive the import.
    if img_mask_process.max() == 1:
        labeled_mask = label(img_mask_process > 0)
    else:
        labeled_mask = img_mask_process.astype(int)

    print("[4/5] Extracting true signal metrics...")
    eff_pixel_size_nm = pixel_size_nm / expansion_factor
    eff_z_step_nm = z_step_nm / expansion_factor

    if mode_a_enabled and not is_3d:
        print("      [Rezim A] Skipped: core-shell FA is available only in 3D mode.")
        mode_a_enabled = False
    elif mode_a_enabled:
        print(
            "      [Rezim A] Enabled: core-shell FA with effective sampling "
            f"(Z,Y,X)=({eff_z_step_nm:.3f}, {eff_pixel_size_nm:.3f}, {eff_pixel_size_nm:.3f}) nm"
        )

    #Regionprops supplies geometry and intensity statistics used to build the
    #stable CSV row schema consumed by postprocessing.py.
    props = regionprops(labeled_mask, intensity_image=img_intensity)
    objects_data = []

    for region in props:
        mean_int = region.intensity_mean
        if region.area < min_voxels:
            continue
            
        centroid_coords = tuple(int(round(c)) for c in region.centroid)
        if len(region.centroid) >= 3:
            z_px, y_px, x_px = [round(float(v), 3) for v in region.centroid[:3]]
        elif len(region.centroid) == 2:
            z_px, y_px, x_px = 0.0, round(float(region.centroid[0]), 3), round(float(region.centroid[1]), 3)
        else:
            z_px = y_px = x_px = np.nan
        
        try:
            cell_id = extruded_mask[centroid_coords] 
        except IndexError:
            continue
            
        if cell_id == 0:
            continue

        row = {
            "filename": tif_path.name,
            "mode": mode,
            "is_3d": is_3d,
            "cell_id": cell_id,
            "object_id": region.label,
            "Z_px": z_px, "Y_px": y_px, "X_px": x_px,
            "mean_intensity": round(mean_int, 2),
            "max_intensity": round(region.intensity_max, 2),
            "integrated_density": round(region.area * mean_int, 2),
        }

        #Convert pixel counts into calibrated biological units while keeping
        #raw counts for auditability and downstream QC.
        if is_3d:
            voxel_volume_bio_um3 = ((eff_pixel_size_nm**2) * eff_z_step_nm) / 1e9
            row["volume_px"] = region.area
            row["volume_bio_um3"] = round(region.area * voxel_volume_bio_um3, 5)
            row["shape_metric_bio"] = row["volume_bio_um3"]
        else:
            pixel_area_bio_um2 = (eff_pixel_size_nm**2) / 1e6
            row["area_px"] = region.area
            row["area_bio_um2"] = round(region.area * pixel_area_bio_um2, 5)
            row["shape_metric_bio"] = row["area_bio_um2"]

        #Optional Rezim A metrics operate on each local region mask and are
        if mode_a_enabled:
            object_mask = region.image.astype(bool)
            touches_edge = _touches_image_edge(region, labeled_mask.shape)
            max_slice_components = _max_components_per_z_slice(object_mask)

            qc_reasons = []
            if touches_edge:
                qc_reasons.append("touches_image_edge")
            if mode_a_exclude_split_slices and max_slice_components > 1:
                qc_reasons.append("split_in_z_slice")

            metrics = compute_core_shell_metrics(
                object_mask,
                sampling=(eff_z_step_nm, eff_pixel_size_nm, eff_pixel_size_nm),
                min_core_voxels=mode_a_min_core_voxels,
                primary_include=not qc_reasons,
                primary_exclusion_reason=";".join(qc_reasons),
            )

            empty_layer_names = [
                name for name in metrics["mode_a_empty_layers"].split(";") if name
            ]
            qc_reasons.extend(f"empty_layer_{name}" for name in empty_layer_names)

            if not metrics["core_valid"]:
                qc_reasons.append("core_below_min_voxels")

            if qc_reasons:
                metrics["mode_a_primary_include"] = False
                metrics["mode_a_qc_reason"] = ";".join(qc_reasons)

            row.update({
                "A_object": metrics["A_object"],
                "A_shell": metrics["A_shell"],
                "A_middle": metrics["A_middle"],
                "A_core": metrics["A_core"],
                "Delta_A_core_shell": metrics["delta_A_core_shell"],
                "A_object_valid": metrics["A_object_valid"],
                "A_shell_valid": metrics["A_shell_valid"],
                "A_middle_valid": metrics["A_middle_valid"],
                "A_core_valid": metrics["A_core_valid"],
                "mode_a_core_voxels": metrics["core_voxels"],
                "mode_a_core_valid": metrics["core_valid"],
                "mode_a_empty_layers": metrics["mode_a_empty_layers"],
                "mode_a_layer_complete_coverage": metrics["layer_qc"]["complete_coverage"],
                "mode_a_object_touches_edge": touches_edge,
                "mode_a_max_components_per_z_slice": max_slice_components,
                "mode_a_primary_include": metrics["mode_a_primary_include"],
                "mode_a_qc_reason": metrics["mode_a_qc_reason"],
                "mode_a_sampling_z_nm": metrics["mode_a_sampling_z_nm"],
                "mode_a_sampling_y_nm": metrics["mode_a_sampling_y_nm"],
                "mode_a_sampling_x_nm": metrics["mode_a_sampling_x_nm"],
                "mode_a_sampling_order": metrics["mode_a_sampling_order"],
                "mode_a_min_core_voxels": metrics["mode_a_min_core_voxels"],
                "mode_a_layer_scheme": "baseline_thirds",
            })

        objects_data.append(row)

    print(f"      Found {len(objects_data)} valid condensates inside ROI.")

    #Preview callbacks are optional, allowing the same numerical function 
    print("[5/5] Generating Preview")
    if send_layer_func:
        send_layer_func({"type": "image", "name": f"Raw Signal ({tif_path.name})", "data": img_intensity, "kwargs": {"colormap": "gray", "blending": "additive"}})
        send_layer_func({"type": "labels", "name": "ROI Boundaries", "data": roi_mask.astype(int), "kwargs": {"opacity": 0.2}})
        send_layer_func({"type": "labels", "name": "Segmentation Mask", "data": labeled_mask, "kwargs": {"opacity": 0.6, "blending": "additive"}})

        coords = [region.centroid for region in props]
        if objects_data and len(coords) > 0:
            if is_3d:
                sizes = [max((3 * region.area / (4 * math.pi))**(1/3) * 2.0, 3.0) for region in props]
            else:
                sizes = [max(math.sqrt(region.area / math.pi) * 2.0, 3.0) for region in props]
                
            send_layer_func({"type": "points", "name": "Detected Condensates", "data": coords, "kwargs": {"size": sizes, "symbol": "disc", "face_color": "yellow" if is_3d else "cyan", "out_of_slice_display": False}})

    return pd.DataFrame(objects_data)
