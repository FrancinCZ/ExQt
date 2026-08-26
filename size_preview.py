from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

from Batch import _prepare_labeled_mask, _resolve_channel_axis, _select_focus_slice

#Return a finite positive float used by the calibration formula
def _positive_number(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a positive number; received {value!r}.") from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and greater than zero; received {value!r}.")
    return result

#Convert one saved mask to the dimensionality used by Batch.py.
def _processing_mask(mask, mode, raw_path=None, signal_channel=1, channel_axis=1):
    mask = np.squeeze(np.asarray(mask))

    if mode == "3d":
        if mask.ndim != 3:
            raise ValueError(f"3d mode needs a 3D mask; received shape {mask.shape}.")
        return mask

    if mode == "2d":
        if mask.ndim == 2:
            return mask
        if mask.ndim == 3:
            return mask.max(axis=0)
        raise ValueError(f"2d mode needs a 2D image or 3D stack; received shape {mask.shape}.")

    if mode != "single_slice":
        raise ValueError(f"Unsupported process mode: {mode!r}.")
    if mask.ndim != 3:
        raise ValueError(f"single_slice mode needs a 3D mask; received shape {mask.shape}.")
    if raw_path is None:
        raise ValueError("single_slice preview needs the matching source TIF.")

    raw = np.squeeze(tifffile.imread(raw_path))
    if raw.ndim == 4:
        resolved_axis = _resolve_channel_axis(raw, expected_axis=channel_axis)
        channel_count = raw.shape[resolved_axis]
        if not 0 <= int(signal_channel) < channel_count:
            raise ValueError(
                f"Signal channel {signal_channel} is outside channel axis {resolved_axis} "
                f"with {channel_count} channels."
            )
        signal = np.take(raw, int(signal_channel), axis=resolved_axis)
    else:
        signal = raw

    signal = np.squeeze(signal)
    if signal.ndim != 3:
        raise ValueError(
            f"single_slice mode needs a 3D signal stack after channel selection; "
            f"received shape {signal.shape}."
        )
    if mask.shape != signal.shape:
        if mask.shape == signal.shape[::-1]:
            mask = np.transpose(mask)
        else:
            raise ValueError(
                f"Source TIF shape {signal.shape} does not match mask shape {mask.shape}."
            )
    return mask[_select_focus_slice(signal)]


def measure_mask_sizes(
    mask,
    mode,
    expansion_factor,
    pixel_size_nm,
    z_step_nm,
    min_voxels=5,
    raw_path=None,
    signal_channel=1,
    channel_axis=1,
):

    expansion_factor = _positive_number(expansion_factor, "Expansion factor")
    pixel_size_nm = _positive_number(pixel_size_nm, "Pixel size XY")
    z_step_nm = _positive_number(z_step_nm, "Z-step")
    min_voxels = int(min_voxels)
    if min_voxels < 1:
        raise ValueError("Raw noise filter must be at least one pixel/voxel.")

    processed_mask = _processing_mask(
        mask,
        mode,
        raw_path=raw_path,
        signal_channel=signal_channel,
        channel_axis=channel_axis,
    )
    labeled_mask = _prepare_labeled_mask(processed_mask)
    is_3d = mode == "3d"
    effective_xy_nm = pixel_size_nm / expansion_factor
    if is_3d:
        effective_z_nm = z_step_nm / expansion_factor
        scale = ((effective_xy_nm ** 2) * effective_z_nm) / 1e9
        unit = "µm³"
        raw_unit = "voxels"
    else:
        scale = (effective_xy_nm ** 2) / 1e6
        unit = "µm²"
        raw_unit = "pixels"

    records = []
    for region in regionprops(labeled_mask):
        raw_size = int(region.area)
        if raw_size < min_voxels:
            continue
        records.append(
            {
                "object_id": int(region.label),
                "raw_size": raw_size,
                "raw_unit": raw_unit,
                "biological_size": round(raw_size * scale, 5),
                "unit": unit,
            }
        )
    return records


def collect_size_preview(
    input_folder,
    mode,
    expansion_factor,
    pixel_size_nm,
    z_step_nm,
    min_voxels=5,
    signal_channel=1,
    channel_axis=1,
):
    folder = Path(input_folder)
    if not folder.exists() or not folder.is_dir():
        raise ValueError("Input folder does not exist or is not a directory.")

    raw_files = [
        path for path in sorted(folder.glob("*.tif"))
        if "Mask" not in path.name and "Final" not in path.name
    ]
    if not raw_files:
        raise ValueError("No source TIF files were found in the selected folder.")

    all_records = []
    for raw_path in raw_files:
        mask_path = folder / f"{raw_path.stem}_Mask.tif"
        if not mask_path.exists():
            raise ValueError(f"Missing matching mask for {raw_path.name}: {mask_path.name}")
        mask = tifffile.imread(mask_path)
        records = measure_mask_sizes(
            mask,
            mode=mode,
            expansion_factor=expansion_factor,
            pixel_size_nm=pixel_size_nm,
            z_step_nm=z_step_nm,
            min_voxels=min_voxels,
            raw_path=raw_path,
            signal_channel=signal_channel,
            channel_axis=channel_axis,
        )
        for record in records:
            record["filename"] = raw_path.name
            record["mask_filename"] = mask_path.name
        all_records.extend(records)

    columns = [
        "filename", "mask_filename", "object_id", "raw_size", "raw_unit",
        "biological_size", "unit",
    ]
    return pd.DataFrame(all_records, columns=columns)
