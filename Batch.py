import math
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import laplace
from skimage.measure import regionprops, label

def _select_focus_slice(volume):
    scores = [laplace(volume[z].astype(np.float32)).var() for z in range(volume.shape[0])]
    return int(np.argmax(scores))

def _resolve_channel_axis(img_raw, expected_axis=1, max_channels=6):
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
    with tifffile.TiffFile(tif_path) as tif:
        try:
            description = tif.pages[0].description
            return None
        except Exception:
            return None


def process_condensates(
    tif_path, mask_path, mode="3d", target_z_slice=None, expansion_factor=1.0,
    min_voxels=5, show_napari=True, pixel_size_nm=None, z_step_nm=None,
    signal_channel=1, dapi_channel=0, channel_axis=1, auto_roi=True,
    send_layer_func=None, request_roi_func=None
):
    tif_path = Path(tif_path)
    mask_path = Path(mask_path)
    print(f"\n[1/5] Loading Pair: {tif_path.name} & {mask_path.name}")

    if send_layer_func:
        send_layer_func({"type": "clear_layers"})

    img_raw = tifffile.imread(tif_path)
    img_mask = tifffile.imread(mask_path)
    img_mask = np.squeeze(img_mask)
    img_raw = np.squeeze(img_raw)

    meta = get_metadata_from_tif(tif_path)
    if meta:
        pixel_size_nm = meta.get('pixel_size', pixel_size_nm)
        z_step_nm = meta.get('z_step', z_step_nm)
        print(f"      [Auto] Metadata načtena: XY={pixel_size_nm}nm, Z={z_step_nm}nm")
    else:
        print(f"      [Info] Metadata nenalezena, používám výchozí hodnoty z GUI.")

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

    if img_mask_process.max() == 1:
        labeled_mask = label(img_mask_process > 0)
    else:
        labeled_mask = img_mask_process.astype(int)

    print("[4/5] Extracting true signal metrics...")
    eff_pixel_size_nm = pixel_size_nm / expansion_factor
    eff_z_step_nm = z_step_nm / expansion_factor

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

        objects_data.append(row)

    print(f"      Found {len(objects_data)} valid condensates inside ROI.")

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
