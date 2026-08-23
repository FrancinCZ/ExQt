import math
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation, laplace
from skimage.measure import regionprops, label
from rezim_a_metrics import MODE_A_LAYER_SCHEME, compute_core_shell_metrics


MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS = 20
MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION = 0.10
MODE_A_Z_SPLIT_PASS_FRACTION = 0.10
MODE_A_Z_SPLIT_REVIEW_FRACTION = 0.30


def _positive_finite(value, name):
    #Return a finite positive float or fail before calibrated measurements
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a positive number; received {value!r}.") from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and greater than zero; received {value!r}.")
    return result


def _physical_value_to_nm(value, unit):
    #Convert common OME/ImageJ physical-length units to nanometres
    if value is None or unit is None:
        return None
    normalized = str(unit).strip().lower().replace("μ", "u").replace("µ", "u")
    factors = {
        "nm": 1.0,
        "nanometer": 1.0,
        "nanometers": 1.0,
        "um": 1_000.0,
        "micron": 1_000.0,
        "microns": 1_000.0,
        "micrometer": 1_000.0,
        "micrometers": 1_000.0,
        "mm": 1_000_000.0,
        "millimeter": 1_000_000.0,
        "millimeters": 1_000_000.0,
        "cm": 10_000_000.0,
        "m": 1_000_000_000.0,
        "inch": 25_400_000.0,
        "in": 25_400_000.0,
    }
    factor = factors.get(normalized)
    if factor is None:
        return None
    converted = float(value) * factor
    return converted if np.isfinite(converted) and converted > 0 else None


def _prepare_labeled_mask(mask):
    #Validate a binary/instance mask and never silently merge 0/255 objects.
    array = np.asarray(mask)
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"Mask must be numeric; received dtype {array.dtype}.")
    if not np.all(np.isfinite(array)):
        raise ValueError("Mask contains NaN or infinite values.")
    if np.any(array < 0):
        raise ValueError("Mask contains negative labels.")
    rounded = np.rint(array)
    if not np.allclose(array, rounded):
        raise ValueError(
            "Mask contains non-integer values. Export a binary mask or an integer instance-label mask, "
            "not probabilities or an interpolated image."
        )

    integer_mask = rounded.astype(np.int64, copy=False)
    positive_labels = np.unique(integer_mask[integer_mask > 0])
    if positive_labels.size <= 1:
        #Covers conventional 0/1 and 0/255 binary masks.
        return label(integer_mask > 0)

    #Multiple positive values are interpreted as instance IDs. Reused IDs on disconnected objects would make regionprops merge unrelated structures.
    for object_id in positive_labels:
        if int(label(integer_mask == object_id).max()) > 1:
            raise ValueError(
                f"Mask label {int(object_id)} occurs in multiple disconnected components. "
                "Export a binary mask for automatic connected-component labeling, or a true "
                "instance-label mask with one unique ID per 3D object."
            )
    return integer_mask


def _touches_image_edge(region, image_shape):
    #Return True when a region bbox reaches any image or stack boundary.
    bbox_min = region.bbox[:region.image.ndim]
    bbox_max = region.bbox[region.image.ndim:]
    return any(start == 0 or end == size for start, end, size in zip(bbox_min, bbox_max, image_shape))


def _touches_roi_edge(region, roi_mask):
    """Return True when an object is directly adjacent to excluded ROI space."""
    roi = np.asarray(roi_mask, dtype=bool)
    if roi.ndim != region.image.ndim:
        raise ValueError("ROI and labeled object must have the same dimensionality.")
    if np.all(roi):
        return False

    ndim = region.image.ndim
    bbox_min = tuple(int(value) for value in region.bbox[:ndim])
    bbox_max = tuple(int(value) for value in region.bbox[ndim:])
    crop_min = tuple(max(0, start - 1) for start in bbox_min)
    crop_max = tuple(min(size, end + 1) for end, size in zip(bbox_max, roi.shape))
    crop_slices = tuple(slice(start, end) for start, end in zip(crop_min, crop_max))

    local_object = np.zeros(tuple(end - start for start, end in zip(crop_min, crop_max)), dtype=bool)
    object_slices = tuple(
        slice(start - crop_start, end - crop_start)
        for start, end, crop_start in zip(bbox_min, bbox_max, crop_min)
    )
    local_object[object_slices] = region.image.astype(bool)
    adjacent = binary_dilation(local_object) & ~local_object
    return bool(np.any(adjacent & ~roi[crop_slices]))


def _max_components_per_z_slice(object_mask):
    #Detect split silhouettes without deciding biological object identity.

    if object_mask.ndim != 3:
        return 1
    return max(int(label(object_mask[z]).max()) for z in range(object_mask.shape[0]))


def _assess_z_split_topology(
    object_mask,
    min_component_voxels=20,
    min_component_fraction=0.10,
    pass_fraction=0.10,
    review_fraction=0.30,
):

    mask = np.asarray(object_mask, dtype=bool)
    if mask.ndim != 3:
        raise ValueError("Z-split topology assessment requires a 3D object mask.")

    min_component_voxels = int(min_component_voxels)
    min_component_fraction = float(min_component_fraction)
    pass_fraction = float(pass_fraction)
    review_fraction = float(review_fraction)
    if min_component_voxels < 1:
        raise ValueError("min_component_voxels must be at least 1.")
    if not 0.0 <= min_component_fraction <= 1.0:
        raise ValueError("min_component_fraction must be between 0 and 1.")
    if not 0.0 <= pass_fraction <= review_fraction <= 1.0:
        raise ValueError("Z-split fractions must satisfy 0 <= pass <= review <= 1.")

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

        minimum_area = max(
            float(min_component_voxels),
            min_component_fraction * foreground_voxels,
        )
        substantial_components = int(np.count_nonzero(component_areas >= minimum_area))
        substantial_max_components = max(
            substantial_max_components,
            substantial_components,
        )
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
        "raw_max_components_per_slice": raw_max_components,
        "substantial_max_components_per_slice": substantial_max_components,
        "occupied_slice_count": occupied_slices,
        "split_slice_count": split_slices,
        "split_slice_fraction": float(split_fraction),
        "status": status,
    }

def _select_focus_slice(volume):
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
    #Read calibrated XY/Z spacing without guessing an unspecified TIFF unit
    with tifffile.TiffFile(tif_path) as tif:
        try:
            meta = {}
            sources = {}

            #OME PhysicalSize fields are the most explicit source because each value carries its own unit.
            if tif.ome_metadata:
                root = ET.fromstring(tif.ome_metadata)
                pixels = next(
                    (item for item in root.iter() if item.tag.rsplit("}", 1)[-1] == "Pixels"),
                    None,
                )
                if pixels is not None:
                    pixel_size = _physical_value_to_nm(
                        pixels.attrib.get("PhysicalSizeX"),
                        pixels.attrib.get("PhysicalSizeXUnit", "um"),
                    )
                    z_step = _physical_value_to_nm(
                        pixels.attrib.get("PhysicalSizeZ"),
                        pixels.attrib.get("PhysicalSizeZUnit", "um"),
                    )
                    if pixel_size is not None:
                        meta["pixel_size"] = pixel_size
                        sources["pixel_size"] = "OME PhysicalSizeX"
                    if z_step is not None:
                        meta["z_step"] = z_step
                        sources["z_step"] = "OME PhysicalSizeZ"

            imagej_metadata = tif.imagej_metadata or {}
            imagej_unit = imagej_metadata.get("unit")
            if imagej_metadata:
                if "z_step" not in meta and "spacing" in imagej_metadata:
                    z_step = _physical_value_to_nm(imagej_metadata["spacing"], imagej_unit)
                    if z_step is not None:
                        meta["z_step"] = z_step
                        sources["z_step"] = f"ImageJ spacing ({imagej_unit})"
            
            tags = tif.pages[0].tags
            if "pixel_size" not in meta and "XResolution" in tags and "ResolutionUnit" in tags:
                numerator, denominator = tags["XResolution"].value
                resolution_unit = tags["ResolutionUnit"].value
                unit_code = getattr(resolution_unit, "value", resolution_unit)
                try:
                    unit_code = int(unit_code)
                except (TypeError, ValueError):
                    unit_code = None
                unit_name = {2: "inch", 3: "cm"}.get(unit_code)
                if unit_name is None and _physical_value_to_nm(1.0, imagej_unit) is not None:
                    unit_name = imagej_unit
                if numerator > 0 and denominator > 0 and unit_name is not None:
                    pixel_size = _physical_value_to_nm(denominator / numerator, unit_name)
                    if pixel_size is not None:
                        meta["pixel_size"] = pixel_size
                        sources["pixel_size"] = f"TIFF XResolution ({unit_name})"

            if not {"pixel_size", "z_step"}.intersection(meta):
                return None
            meta["sources"] = sources
            meta["complete"] = "pixel_size" in meta and "z_step" in meta
            return meta
        except (ET.ParseError, TypeError, ValueError, KeyError, IndexError):
            return None


def process_condensates(
    tif_path, mask_path, mode="3d", target_z_slice=None, expansion_factor=1.0,
    min_voxels=5, show_napari=True, pixel_size_nm=None, z_step_nm=None,
    signal_channel=1, dapi_channel=0, channel_axis=1, auto_roi=True,
    send_layer_func=None, request_roi_func=None, mode_a_enabled=False,
    mode_a_min_core_voxels=20, mode_a_exclude_split_slices=True,
    mode_a_z_split_min_component_voxels=MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS,
    mode_a_z_split_min_component_fraction=MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION,
    mode_a_z_split_pass_fraction=MODE_A_Z_SPLIT_PASS_FRACTION,
    mode_a_z_split_review_fraction=MODE_A_Z_SPLIT_REVIEW_FRACTION,
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

    pixel_size_nm = _positive_finite(pixel_size_nm, "Pixel size XY")
    z_step_nm = _positive_finite(z_step_nm, "Z-step")
    expansion_factor = _positive_finite(expansion_factor, "Expansion factor")

    #Channel selection is isolated here - downstream processing works with a
    #single intensity volume regardless of the original TIFF layout.
    if img_raw.ndim == 4:
        ch_axis = _resolve_channel_axis(img_raw, expected_axis=channel_axis)
        channel_count = img_raw.shape[ch_axis]
        if not 0 <= int(signal_channel) < channel_count:
            raise ValueError(
                f"Signal channel {signal_channel} is outside channel axis {ch_axis} "
                f"with {channel_count} channels."
            )
        if not 0 <= int(dapi_channel) < channel_count:
            raise ValueError(
                f"DAPI channel {dapi_channel} is outside channel axis {ch_axis} "
                f"with {channel_count} channels."
            )
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

    #Convert the requested mode into one processing image/mask pair. 
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
    #A manual ROI is requested through callbacks so the worker can pause and wait for user input. The ROI is applied to the mask before regionprops.
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

    #Binary 0/1 and 0/255 masks are labeled consistently.
    labeled_mask = _prepare_labeled_mask(img_mask_process)

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

    #Regionprops supplies geometry and intensity statistics used to build the stable CSV row schema consumed by postprocessing.py.
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

        #Convert pixel counts into calibrated biological units while keeping raw counts for auditability and downstream QC.
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
            touches_roi_edge = _touches_roi_edge(region, roi_mask)
            topology = _assess_z_split_topology(
                object_mask,
                min_component_voxels=mode_a_z_split_min_component_voxels,
                min_component_fraction=mode_a_z_split_min_component_fraction,
                pass_fraction=mode_a_z_split_pass_fraction,
                review_fraction=mode_a_z_split_review_fraction,
            )
            max_slice_components = topology["raw_max_components_per_slice"]

            qc_reasons = []
            if touches_edge:
                qc_reasons.append("touches_image_edge")
            if touches_roi_edge:
                qc_reasons.append("touches_roi_edge")
            if mode_a_exclude_split_slices and topology["status"] != "pass":
                qc_reasons.append(f"z_topology_{topology['status']}")

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

            for layer_name in ("object", "shell", "middle"):
                if (
                    layer_name not in empty_layer_names
                    and not metrics[f"A_{layer_name}_valid"]
                ):
                    qc_reasons.append(f"invalid_anisotropy_{layer_name}")

            if not metrics["layer_qc"]["complete_coverage"]:
                qc_reasons.append("layer_coverage_error")

            if qc_reasons:
                metrics["mode_a_primary_include"] = False
                metrics["mode_a_qc_reason"] = ";".join(qc_reasons)

            row.update({
                "A_object": metrics["A_object"],
                "A_shell": metrics["A_shell"],
                "A_middle": metrics["A_middle"],
                "A_core": metrics["A_core"],
                "Delta_A_middle_shell": metrics["delta_A_middle_shell"],
                "Delta_A_core_middle": metrics["delta_A_core_middle"],
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
                "mode_a_object_touches_roi_edge": touches_roi_edge,
                "mode_a_max_components_per_z_slice": max_slice_components,
                "mode_a_z_topology_status": topology["status"],
                "mode_a_z_occupied_slices": topology["occupied_slice_count"],
                "mode_a_z_split_slices": topology["split_slice_count"],
                "mode_a_z_split_slice_fraction": topology["split_slice_fraction"],
                "mode_a_z_max_substantial_components_per_slice": topology["substantial_max_components_per_slice"],
                "mode_a_z_split_min_component_voxels": int(mode_a_z_split_min_component_voxels),
                "mode_a_z_split_min_component_fraction": float(mode_a_z_split_min_component_fraction),
                "mode_a_z_split_pass_fraction": float(mode_a_z_split_pass_fraction),
                "mode_a_z_split_review_fraction": float(mode_a_z_split_review_fraction),
                "mode_a_primary_include": metrics["mode_a_primary_include"],
                "mode_a_qc_reason": metrics["mode_a_qc_reason"],
                "mode_a_sampling_z_nm": metrics["mode_a_sampling_z_nm"],
                "mode_a_sampling_y_nm": metrics["mode_a_sampling_y_nm"],
                "mode_a_sampling_x_nm": metrics["mode_a_sampling_x_nm"],
                "mode_a_sampling_order": metrics["mode_a_sampling_order"],
                "mode_a_min_core_voxels": metrics["mode_a_min_core_voxels"],
                "mode_a_layer_scheme": MODE_A_LAYER_SCHEME,
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
