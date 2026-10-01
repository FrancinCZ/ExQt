import math
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation, laplace, distance_transform_edt
from skimage.measure import regionprops, label
from rezim_a_metrics import MODE_A_LAYER_SCHEME, UNCLASSIFIED_QC, compute_core_shell_metrics
from defaults import DEFAULT_SETTINGS


MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS = DEFAULT_SETTINGS["mode_a_z_split_min_component_voxels"]
MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION = DEFAULT_SETTINGS["mode_a_z_split_min_component_fraction"]
MODE_A_Z_SPLIT_PASS_FRACTION = DEFAULT_SETTINGS["mode_a_z_split_pass_fraction"]
MODE_A_Z_SPLIT_REVIEW_FRACTION = DEFAULT_SETTINGS["mode_a_z_split_review_fraction"]


def _positive_finite(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a positive number; received {value!r}.") from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and greater than zero; received {value!r}.")
    return result


def _non_negative_finite(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a number >= 0; received {value!r}.") from error
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be finite and >= 0; received {value!r}.")
    return result


def _physical_value_to_nm(value, unit):
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


def source_tiff_files(folder):
    """Source TIFFs of a folder: masks, ROIs, alignment sidecars and "final" exports are skipped."""
    return sorted(
        path for path in Path(folder).iterdir()
        if path.is_file()
        and path.suffix.lower() in {".tif", ".tiff"}
        and not path.stem.lower().endswith(("_mask", "_roi"))
        and "final" not in path.stem.lower()
        and "_alignment_" not in path.stem.lower()
    )


def matching_mask_path(source_path):
    """The conventional <stem>_Mask.tif/.tiff next to a source TIFF (the .tif path if neither exists)."""
    for suffix in (".tif", ".tiff"):
        candidate = source_path.with_name(f"{source_path.stem}_Mask{suffix}")
        if candidate.is_file():
            return candidate
    return source_path.with_name(f"{source_path.stem}_Mask.tif")


#Where the ROI of one image is saved: the run's output folder if given, else next to the mask.
def roi_output_path(tif_path, mask_path, roi_output_dir=None):
    folder = Path(roi_output_dir) if roi_output_dir is not None else Path(mask_path).parent
    return folder / f"{Path(tif_path).stem}_ROI.tif"


#Validate a binary/instance mask and never silently merge 0/255 objects.
def _prepare_labeled_mask(mask):
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


#Return True when a region reaches the array boundary or, for aligned stacks, the border of
#the real data (alignment padding marked 0 in the valid-voxel mask).
def _touches_image_edge(region, image_shape, valid_voxels=None):
    bbox_min = region.bbox[:region.image.ndim]
    bbox_max = region.bbox[region.image.ndim:]
    if any(start == 0 or end == size for start, end, size in zip(bbox_min, bbox_max, image_shape)):
        return True
    return valid_voxels is not None and _touches_roi_edge(region, valid_voxels)


def _touches_roi_edge(region, roi_mask):
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


#Split a drawn ROI stack into per-nucleus groups of (z, 2D mask) shapes.
def _group_roi_shapes(mask_3d):
    # napari Shapes.to_labels gives every shape its own ID (index + 1), not a nucleus ID,
    # so nuclei are recovered from geometry: one node per connected region of each ID on
    # each slice; nodes on different slices are linked when they overlap in XY, nodes on
    # the same slice when they touch (overlapping shapes are painted over each other).
    nodes = []
    for z in np.where(mask_3d.any(axis=(1, 2)))[0]:
        plane = mask_3d[z]
        for value in np.unique(plane[plane > 0]):
            components = label(plane == value)
            for component_id in range(1, int(components.max()) + 1):
                region = components == component_id
                ys, xs = np.nonzero(region)
                nodes.append((int(z), region, (ys.min(), ys.max(), xs.min(), xs.max())))

    parent = list(range(len(nodes)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for i, (z_i, mask_i, box_i) in enumerate(nodes):
        grown_i = None
        for j in range(i + 1, len(nodes)):
            z_j, mask_j, box_j = nodes[j]
            margin = 1 if z_i == z_j else 0
            if (box_i[0] > box_j[1] + margin or box_j[0] > box_i[1] + margin
                    or box_i[2] > box_j[3] + margin or box_j[2] > box_i[3] + margin):
                continue
            if z_i == z_j:
                if grown_i is None:
                    grown_i = binary_dilation(mask_i)
                linked = np.any(grown_i & mask_j)
            else:
                linked = np.any(mask_i & mask_j)
            if linked:
                parent[find(j)] = find(i)

    groups = {}
    for index, (z, mask, _) in enumerate(nodes):
        groups.setdefault(find(index), []).append((z, mask))

    ordered = []
    for group in groups.values():
        group.sort(key=lambda node: node[0])
        slices = [z for z, _ in group]
        duplicates = sorted({z for z in slices if slices.count(z) > 1})
        if duplicates:
            raise ValueError(
                f"ROI shapes on the same Z-slice ({', '.join(map(str, duplicates))}) overlap or touch, "
                "so they cannot be assigned to separate nuclei. Redraw them so they do not touch, "
                "or draw one shape per nucleus per slice."
            )
        first_z, first_mask = group[0]
        ordered.append(((first_z, int(np.flatnonzero(first_mask)[0])), group))
    ordered.sort(key=lambda item: item[0])
    return [group for _, group in ordered]


#Raise when two nucleus ROIs claim the same voxel.
def _assert_roi_groups_disjoint(group_masks):
    claimed = None
    for index, mask in enumerate(group_masks, start=1):
        if claimed is None:
            claimed = np.zeros(mask.shape, dtype=bool)
        if np.any(claimed & mask):
            raise ValueError(
                f"The interpolated ROI of nucleus {index} overlaps another nucleus ROI; "
                "the nuclei cannot be separated unambiguously."
            )
        claimed |= mask


def _interpolate_or_extrude_roi(mask_3d: np.ndarray) -> np.ndarray:
    """
    Process a user-drawn 3D ROI mask into a per-nucleus label volume (cell_id 1..N):
    - Shapes are grouped into nuclei by XY overlap across slices (see _group_roi_shapes).
    - A nucleus drawn on multiple Z-slices is interpolated between its own slices with
      signed distance fields; slices outside its [z_min, z_max] stay zero (cuts off noise).
    - A nucleus drawn on exactly 1 Z-slice is extruded across all Z slices
      (preserves the single-slice workflow).
    With one nucleus (ExM) the result equals the former binary implementation.
    """
    if mask_3d.ndim != 3 or mask_3d.max() == 0:
        return mask_3d

    group_masks = []
    for group in _group_roi_shapes(mask_3d):
        group_mask = np.zeros(mask_3d.shape, dtype=bool)
        if len(group) == 1:
            group_mask[:] = group[0][1]
        else:
            for (z0, m0), (z1, m1) in zip(group, group[1:]):
                # 2D EDT within one XY slice needs no `sampling`: XY pixels are isotropic.
                # Z enters only through the linear blend weight (slice index), not the EDT,
                # so this is not the anisotropic 3D EDT case.
                d0 = distance_transform_edt(m0) - distance_transform_edt(~m0)
                d1 = distance_transform_edt(m1) - distance_transform_edt(~m1)
                for z in range(z0, z1 + 1):
                    w1 = (z - z0) / (z1 - z0)
                    group_mask[z] |= ((1.0 - w1) * d0 + w1 * d1) > 0
        group_masks.append(group_mask)

    _assert_roi_groups_disjoint(group_masks)
    out = np.zeros(mask_3d.shape, dtype=np.int32)
    for cell_id, group_mask in enumerate(group_masks, start=1):
        out[group_mask] = cell_id
    return out


#Detect split silhouettes without deciding biological object identity.
def _max_components_per_z_slice(object_mask):
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

#Validate or infer the channel axis used to extract signal images.
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

def _squeeze_with_axes(array, axes):
    data = np.asarray(array)
    normalized_axes = str(axes).upper()
    if len(normalized_axes) != data.ndim:
        return np.squeeze(data), ""
    keep = [index for index, size in enumerate(data.shape) if size != 1]
    if len(keep) == data.ndim:
        return data, normalized_axes
    if not keep:
        return np.squeeze(data), ""
    return np.squeeze(data), "".join(normalized_axes[index] for index in keep)


def _read_tiff_with_axes(tif_path):
    with tifffile.TiffFile(tif_path) as tif:
        if len(tif.series) != 1:
            raise ValueError(f"TIFF must contain exactly one series: {Path(tif_path).name}")
        series = tif.series[0]
        return _squeeze_with_axes(series.asarray(), series.axes)


def _alignment_drift_csv(tif_path):
    path = Path(tif_path)
    stem = path.name
    for suffix in (".ome.tiff", ".ome.tif", ".tiff", ".tif"):
        if stem.lower().endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    return path.with_name(f"{stem}_drift.csv")


#Load the aligner's valid-voxel sidecar (1 = real data, 0 = zero padding from shifting).
def _load_valid_voxel_mask(tif_path, expected_shape):
    drift_csv = _alignment_drift_csv(tif_path)
    valid_path = drift_csv.with_name(drift_csv.name.replace("_drift.csv", "_alignment_valid.tif"))
    if not valid_path.is_file():
        if drift_csv.is_file():
            print(
                f"      [Alignment] Warning: {Path(tif_path).name} is aligned but has no valid-voxel mask "
                f"({valid_path.name}); re-align with the current ExQt so padding is excluded from statistics."
            )
        return None
    valid = np.asarray(tifffile.imread(valid_path)) > 0
    if valid.shape != tuple(expected_shape):
        raise ValueError(
            f"Alignment valid-voxel mask {valid_path.name} has shape {valid.shape}, "
            f"but the image has shape {tuple(expected_shape)}."
        )
    return valid


def validate_alignment_qc(tif_path):
    """Return the aligner's overall status ("PASS" / "REVIEW") for an aligned input, None when
    the input is not aligned; raise ValueError when the stack must be excluded."""
    drift_csv = _alignment_drift_csv(tif_path)
    if not drift_csv.is_file():
        return None
    try:
        table = pd.read_csv(drift_csv)
    except Exception as error:
        raise ValueError(f"Cannot read alignment QC for {Path(tif_path).name}: {error}") from error
    required = {"status", "step_shift_y_px", "step_shift_x_px", "step_magnitude_px"}
    if table.empty or not required.issubset(table.columns):
        raise ValueError(f"Alignment QC is incomplete for {Path(tif_path).name}: {drift_csv.name}")

    numeric = table[["step_shift_y_px", "step_shift_x_px", "step_magnitude_px"]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError(f"Alignment QC contains non-finite shifts for {Path(tif_path).name}")

    statuses = set(table["status"].astype(str))
    if "overall_status" in table.columns:
        overall_values = set(table["overall_status"].astype(str))
        if len(overall_values) != 1 or not overall_values <= {"PASS", "REVIEW", "FAIL"}:
            raise ValueError(f"Alignment QC has an invalid overall_status for {Path(tif_path).name}: {sorted(overall_values)}")
        overall = overall_values.pop()
    elif "FAIL" in statuses:
        # Drift CSV from an older ExQt: the overall verdict is unknown, so stay strict.
        raise ValueError(
            f"Alignment QC for {Path(tif_path).name} has FAIL steps and no overall_status "
            "(older aligner); the stack is excluded from statistics. Please re-align it."
        )
    else:
        overall = "REVIEW" if "REVIEW" in statuses else "PASS"

    if overall == "FAIL":
        raise ValueError(
            f"Alignment QC is FAIL for {Path(tif_path).name}; "
            "the stack is excluded from statistics."
        )
    if overall == "REVIEW":
        print(f"      [Alignment QC] {Path(tif_path).name} is REVIEW: analysed, but excluded from primary statistics.")
    return overall



#Read calibrated XY/Z spacing without guessing an unspecified TIFF unit.
def get_metadata_from_tif(tif_path):

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
            #Inch/cm ResolutionUnit is print DPI (e.g. default 72 dpi), not a microscope calibration; XResolution is trusted only with an explicit ImageJ length unit.
            imagej_length_unit = _physical_value_to_nm(1.0, imagej_unit) is not None
            if "pixel_size" not in meta and "XResolution" in tags and imagej_length_unit:
                numerator, denominator = tags["XResolution"].value
                if numerator > 0 and denominator > 0:
                    pixel_size = _physical_value_to_nm(denominator / numerator, imagej_unit)
                    if pixel_size is not None:
                        meta["pixel_size"] = pixel_size
                        sources["pixel_size"] = f"TIFF XResolution ({imagej_unit})"

            if not {"pixel_size", "z_step"}.intersection(meta):
                return None
            meta["sources"] = sources
            meta["complete"] = "pixel_size" in meta and "z_step" in meta
            return meta
        except (ET.ParseError, TypeError, ValueError, KeyError, IndexError):
            return None


def process_condensates(
    tif_path, mask_path, mode="3d", target_z_slice=None, expansion_factor=1.0,
    min_voxels=DEFAULT_SETTINGS["raw_min_voxels"], show_napari=True, pixel_size_nm=None, z_step_nm=None,
    signal_channel=1, dapi_channel=0, channel_axis=1, auto_roi=True,
    send_layer_func=None, request_roi_func=None, mode_a_enabled=False,
    mode_a_min_core_voxels=DEFAULT_SETTINGS["mode_a_min_core_voxels"], mode_a_exclude_split_slices=True,
    mode_a_z_split_min_component_voxels=MODE_A_Z_SPLIT_MIN_COMPONENT_VOXELS,
    mode_a_z_split_min_component_fraction=MODE_A_Z_SPLIT_MIN_COMPONENT_FRACTION,
    mode_a_z_split_pass_fraction=MODE_A_Z_SPLIT_PASS_FRACTION,
    mode_a_z_split_review_fraction=MODE_A_Z_SPLIT_REVIEW_FRACTION,
    detector_offset_adu=0.0,
    detector_offset_source="manual",
    roi_output_dir=None,
):
    #Process one raw/mask TIFF pair and return one row per valid object.

    tif_path = Path(tif_path)
    mask_path = Path(mask_path)
    alignment_status = validate_alignment_qc(tif_path) or "not_aligned"
    print(f"\n[1/5] Loading Pair: {tif_path.name} & {mask_path.name}")

    if send_layer_func:
        send_layer_func({"type": "clear_layers"})

    img_raw, raw_axes = _read_tiff_with_axes(tif_path)
    img_mask, mask_axes = _read_tiff_with_axes(mask_path)

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
    detector_offset_adu = _non_negative_finite(detector_offset_adu, "Detector offset")

    #Channel selection is isolated here - downstream processing works with a
    #single intensity volume regardless of the original TIFF layout.
    if img_raw.ndim == 4:
        ch_axis = raw_axes.index("C") if raw_axes.count("C") == 1 else _resolve_channel_axis(img_raw, expected_axis=channel_axis)

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
        raw_spatial = "".join(a for a in raw_axes if a != "C")
    elif img_raw.ndim == 3 and "C" in raw_axes:
        raise ValueError(
            f"{tif_path.name}: 2D multi-channel images (axes {raw_axes}) are not supported yet. "
            "Save the signal channel as a single-channel TIFF, or acquire a Z-stack."
        )
    else:
        img_dapi = img_raw
        img_intensity = img_raw
        raw_spatial = raw_axes

    if raw_spatial and set(raw_spatial) == {"Z", "Y", "X"} and raw_spatial != "ZYX":
        raw_perm = [raw_spatial.index(a) for a in "ZYX"]
        img_intensity = np.transpose(img_intensity, raw_perm)
        img_dapi = np.transpose(img_dapi, raw_perm)

    if mask_axes and set(mask_axes) == {"Z", "Y", "X"} and mask_axes != "ZYX":
        mask_perm = [mask_axes.index(a) for a in "ZYX"]
        img_mask = np.transpose(img_mask, mask_perm)
        mask_axes = "ZYX"
    elif mask_axes and set(mask_axes) == {"Y", "X"} and mask_axes != "YX":
        mask_perm = [mask_axes.index(a) for a in "YX"]
        img_mask = np.transpose(img_mask, mask_perm)
        mask_axes = "YX"

    if img_mask.shape != img_intensity.shape:
        if img_mask.shape == img_intensity.shape[::-1]:
            img_mask = np.transpose(img_mask)
            print("Mask was automatically transposed.")
        else:
            raise ValueError(f"CRITICAL ERROR: TIF shape {img_intensity.shape} does not match "
                            f"Mask shape {img_mask.shape}. Ensure your segmentation software outputs correct dimensions.")

    is_stack = img_intensity.ndim == 3
    valid_voxels = _load_valid_voxel_mask(tif_path, img_intensity.shape)
    print(f"      Mode: {mode}")

    if mode == "single_slice":
        if not is_stack:
            raise ValueError("single_slice mode needs a 3D (Z,Y,X) stack.")
        z_idx = (target_z_slice if target_z_slice is not None else _select_focus_slice(img_intensity))
        print(f"      Using Z-slice {z_idx} (auto-focus)")
        img_intensity = img_intensity[z_idx]
        img_dapi_process = img_dapi[z_idx]
        img_mask_process = img_mask[z_idx]
        if valid_voxels is not None:
            valid_voxels = valid_voxels[z_idx]
        is_3d = False

    elif mode == "2d":
        img_intensity = img_intensity if not is_stack else img_intensity.max(axis=0)
        img_dapi_process = img_dapi if not is_stack else img_dapi.max(axis=0)
        img_mask_process = img_mask if not is_stack else img_mask.max(axis=0)
        if valid_voxels is not None and is_stack:
            # A projected pixel is valid only if it is valid in every slice.
            valid_voxels = valid_voxels.all(axis=0)
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
    roi_source = "auto_fov" if (auto_roi or request_roi_func is None) else "manual"
    if auto_roi or request_roi_func is None:
        roi_mask = np.ones_like(img_mask_process, dtype=bool)
        extruded_mask = np.ones_like(img_mask_process, dtype=int)
    else:
        if send_layer_func:
            # In ExM, Z-step is typically larger than XY pixel size — pass scale so the ROI preview is physically accurate.
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
        if is_3d and extruded_mask.max() > 0:
            extruded_mask = _interpolate_or_extrude_roi(extruded_mask)
        roi_mask = extruded_mask > 0

    #Persist the ROI mask so auditors can verify which region was selected. With roi_output_dir (the run's output folder) every run keeps its own ROI and the input folder is not touched.
    roi_save_path = roi_output_path(tif_path, mask_path, roi_output_dir)
    roi_save_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        tifffile.imwrite(str(roi_save_path), extruded_mask.astype(np.uint8), compression="zlib")
        print(f"      [ROI] Saved ROI mask -> {roi_save_path.name}")
    except Exception as e:
        print(f"      [ROI] Warning: could not save ROI mask: {e}")

    print("\n[3/5] Applying provided Mask...")

    img_mask_process = img_mask_process * roi_mask

    labeled_mask = _prepare_labeled_mask(img_mask_process)

    print("[4/5] Extracting true signal metrics...")
    eff_pixel_size_nm = pixel_size_nm / expansion_factor
    eff_z_step_nm = z_step_nm / expansion_factor

    if mode_a_enabled and not is_3d:
        print("      [Radial FA Profiling] Skipped: radial FA is available only in 3D mode.")
        mode_a_enabled = False
    elif mode_a_enabled:
        print(
            "      [Radial FA Profiling] Enabled with effective sampling "
            f"(Z,Y,X)=({eff_z_step_nm:.3f}, {eff_pixel_size_nm:.3f}, {eff_pixel_size_nm:.3f}) nm"
        )

    #The detector offset is an explicit acquisition setting (default 0 = no subtraction).
    #It is never estimated from ROI pixels: a low percentile depends on ROI shape and noise.
    k_offset = detector_offset_adu

    #Precompute nucleoplasm (background within cell ROI across active Z-slices) mean intensity per cell_id.
    #Auto-ROI has no nucleus outline, so the "nucleoplasm" would include gel/background: no K there.
    cell_ids = np.unique(extruded_mask[extruded_mask > 0]) if roi_source != "auto_fov" else []
    nucleoplasm_means = {}
    for cid in cell_ids:
        bg_mask = (extruded_mask == cid) & (labeled_mask == 0)
        if valid_voxels is not None:
            bg_mask &= valid_voxels

        if is_3d:
            #Per-cell Z-bounds: only Z-slices where THIS cell has labeled objects
            cell_objects_mask = (extruded_mask == cid) & (labeled_mask > 0)
            if cell_objects_mask.any():
                z_occupied_cell = np.where(np.any(cell_objects_mask, axis=(1, 2)))[0]
                z_min_cell = int(z_occupied_cell.min())
                z_max_cell = int(z_occupied_cell.max())
                bg_mask[:z_min_cell, :, :] = False
                bg_mask[z_max_cell + 1:, :, :] = False

        if np.any(bg_mask):
            nucleoplasm_means[cid] = float(np.mean(img_intensity[bg_mask]))
        else:
            global_bg = (extruded_mask == cid) & (labeled_mask == 0)
            if valid_voxels is not None:
                global_bg &= valid_voxels
            nucleoplasm_means[cid] = float(np.mean(img_intensity[global_bg])) if np.any(global_bg) else np.nan



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

        bg_int = nucleoplasm_means.get(cell_id, np.nan)
        net_mean_int = max(mean_int - k_offset, 0.0)
        if roi_source == "auto_fov":
            net_bg_int = np.nan
            k_valid = False
            k_invalid_reason = "auto_roi_no_nucleoplasm"
        elif not np.isfinite(bg_int):
            net_bg_int = np.nan
            k_valid = False
            k_invalid_reason = "no_nucleoplasm"
        else:
            net_bg_int = bg_int - k_offset
            if net_bg_int <= 0:
                net_bg_int = np.nan
                k_valid = False
                k_invalid_reason = "negative_denominator"
            else:
                k_valid = True
                k_invalid_reason = ""
        part_coeff = round(net_mean_int / net_bg_int, 2) if np.isfinite(net_bg_int) and net_bg_int > 0 else np.nan


        row = {
            "filename": tif_path.name,
            "mode": mode,
            "is_3d": is_3d,
            "cell_id": cell_id,
            "roi_source": roi_source,
            "alignment_status": alignment_status,
            "object_id": region.label,
            "Z_px": z_px, "Y_px": y_px, "X_px": x_px,
            "mean_intensity": round(mean_int, 2),
            "max_intensity": round(region.intensity_max, 2),
            "integrated_density": round(region.area * mean_int, 2),
            "nucleoplasm_mean_intensity": round(bg_int, 2) if np.isfinite(bg_int) else np.nan,
            "partition_coefficient": part_coeff,
            "K_offset_adu": float(k_offset),
            "K_offset_method": "explicit_setting",
            "K_offset_source": str(detector_offset_source),
            "K_valid": k_valid,
            "K_invalid_reason": k_invalid_reason if not k_valid else "",
        }

        #Truncated objects (image/stack border, alignment padding, ROI boundary) are flagged in every mode; postprocessing keeps them out of primary.
        touches_edge = _touches_image_edge(region, labeled_mask.shape, valid_voxels)
        touches_roi_edge = _touches_roi_edge(region, roi_mask)
        row["touches_image_edge"] = touches_edge
        row["touches_roi_edge"] = touches_roi_edge

        row["applied_pixel_size_nm"] = float(pixel_size_nm)
        row["applied_z_step_nm"] = float(z_step_nm) if is_3d else float("nan")
        if is_3d:
            voxel_volume_bio_um3 = ((eff_pixel_size_nm**2) * eff_z_step_nm) / 1e9
            voxel_volume_gel_um3 = ((float(pixel_size_nm)**2) * float(z_step_nm)) / 1e9
            row["volume_px"] = region.area
            row["volume_bio_um3"] = float(region.area * voxel_volume_bio_um3)
            row["volume_gel_um3"] = float(region.area * voxel_volume_gel_um3)
            row["shape_metric_bio"] = row["volume_bio_um3"]
        else:
            pixel_area_bio_um2 = (eff_pixel_size_nm**2) / 1e6
            pixel_area_gel_um2 = (float(pixel_size_nm)**2) / 1e6
            row["area_px"] = region.area
            row["area_bio_um2"] = float(region.area * pixel_area_bio_um2)
            row["area_gel_um2"] = float(region.area * pixel_area_gel_um2)
            row["shape_metric_bio"] = row["area_bio_um2"]



        if mode_a_enabled:
            object_mask = region.image.astype(bool)
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
                intensity_image=region.image_intensity if hasattr(region, "image_intensity") else getattr(region, "intensity_image", None),
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
                # QC reasons found after the metrics were computed also void the class.
                metrics["condensate_class"] = UNCLASSIFIED_QC

            row.update({
                "A_object": metrics["A_object"],
                "A_shell": metrics["A_shell"],
                "A_middle": metrics["A_middle"],
                "A_core": metrics["A_core"],
                "Delta_A_middle_shell": metrics["delta_A_middle_shell"],
                "Delta_A_core_middle": metrics["delta_A_core_middle"],
                "Delta_A_core_shell": metrics["delta_A_core_shell"],
                "mean_intensity_object": metrics["mean_intensity_object"],
                "mean_intensity_shell": metrics["mean_intensity_shell"],
                "mean_intensity_middle": metrics["mean_intensity_middle"],
                "mean_intensity_core": metrics["mean_intensity_core"],
                "Delta_intensity_middle_shell": metrics["delta_intensity_middle_shell"],
                "Delta_intensity_core_middle": metrics["delta_intensity_core_middle"],
                "Delta_intensity_core_shell": metrics["delta_intensity_core_shell"],
                "condensate_class": metrics.get("condensate_class", "Unclassified"),
                "gradient_t_core_shell": metrics["gradient_t_core_shell"],
                "gradient_p_core_shell": metrics["gradient_p_core_shell"],
                "gradient_significant": metrics["gradient_significant"],
                "gradient_test": metrics["gradient_test"],
                "gradient_alpha": metrics["gradient_alpha"],
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
                "principal_std_1_nm": metrics.get("principal_std_1_nm", np.nan),
                "principal_std_2_nm": metrics.get("principal_std_2_nm", np.nan),
                "principal_std_3_nm": metrics.get("principal_std_3_nm", np.nan),
                "null_FA_object": metrics.get("null_FA_object", np.nan),
                "null_FA_shell": metrics.get("null_FA_shell", np.nan),
                "null_FA_middle": metrics.get("null_FA_middle", np.nan),
                "null_FA_core": metrics.get("null_FA_core", np.nan),
                "null_delta_FA_core_shell": metrics.get("null_delta_FA_core_shell", np.nan),
                "null_valid": metrics.get("null_valid", False),
                "null_invalid_reason": metrics.get("null_invalid_reason", ""),
                "delta_FA_excess": metrics.get("delta_FA_excess", np.nan),
            })

        objects_data.append(row)

    print(f"      Found {len(objects_data)} valid condensates inside ROI.")

    print("[5/5] Generating Preview")
    if send_layer_func:
        z_scale = z_step_nm / pixel_size_nm if is_3d else 1.0
        layer_scale = (z_scale, 1, 1) if is_3d else None

        img_kwargs: dict = {"colormap": "gray", "blending": "additive"}
        if layer_scale is not None:
            img_kwargs["scale"] = layer_scale

        lbl_kwargs_roi: dict = {"opacity": 0.2}
        lbl_kwargs_seg: dict = {"opacity": 0.6, "blending": "additive"}
        if layer_scale is not None:
            lbl_kwargs_roi["scale"] = layer_scale
            lbl_kwargs_seg["scale"] = layer_scale

        send_layer_func({"type": "image",  "name": f"Raw Signal ({tif_path.name})", "data": img_intensity,  "kwargs": img_kwargs})
        send_layer_func({"type": "labels", "name": "ROI Boundaries",                "data": roi_mask.astype(int), "kwargs": lbl_kwargs_roi})
        send_layer_func({"type": "labels", "name": "Segmentation Mask",             "data": labeled_mask,  "kwargs": lbl_kwargs_seg})

        #Only objects that made it into the table (min_voxels, cell_id checks) are marked.
        measured_ids = {row["object_id"] for row in objects_data}
        measured = [region for region in props if region.label in measured_ids]
        if measured:
            coords = [region.centroid for region in measured]
            if is_3d:
                sizes = [max((3 * region.area / (4 * math.pi))**(1/3) * 2.0, 3.0) for region in measured]
            else:
                sizes = [max(math.sqrt(region.area / math.pi) * 2.0, 3.0) for region in measured]
            points_kwargs = {"size": sizes, "symbol": "disc", "face_color": "yellow" if is_3d else "cyan", "out_of_slice_display": False}
            if layer_scale is not None:
                #Same Z scale as the image/labels, otherwise 3D points sit at the wrong height.
                points_kwargs["scale"] = layer_scale
            send_layer_func({"type": "points", "name": "Detected Condensates", "data": coords, "kwargs": points_kwargs})

    return pd.DataFrame(objects_data)