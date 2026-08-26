from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

from Batch import (
    _positive_finite,
    _prepare_labeled_mask,
    _resolve_channel_axis,
    _select_focus_slice,
)


PREVIEW_COLUMNS = [
    "filename",
    "mask_filename",
    "mode",
    "object_id",
    "raw_size",
    "biological_size",
    "unit",
]


def _empty_preview_table():
    """Return an empty result with the same schema as a populated preview."""
    return pd.DataFrame(columns=PREVIEW_COLUMNS)


def _load_processing_mask(source_path, mask_path, mode, signal_channel):
    """Load a source/mask pair and reproduce Batch.process_condensates geometry."""
    image = np.squeeze(tifffile.imread(source_path))
    mask = np.squeeze(tifffile.imread(mask_path))

    if image.ndim == 4:
        channel_axis = _resolve_channel_axis(image, expected_axis=1)
        channel_count = image.shape[channel_axis]
        if not 0 <= signal_channel < channel_count:
            raise ValueError(
                f"Signal channel {signal_channel} is outside channel axis "
                f"{channel_axis} with {channel_count} channels for {source_path.name}."
            )
        image = np.take(image, signal_channel, axis=channel_axis)

    if image.ndim not in (2, 3):
        raise ValueError(
            f"Unsupported source shape {image.shape} for {source_path.name}; "
            "the selected signal image must be 2D or 3D."
        )

    if mask.shape != image.shape:
        if mask.shape == image.shape[::-1]:
            mask = np.transpose(mask)
        else:
            raise ValueError(
                f"TIF shape {image.shape} does not match mask shape {mask.shape} "
                f"for {source_path.name}."
            )

    is_stack = image.ndim == 3
    if mode == "single_slice":
        if not is_stack:
            raise ValueError(
                f"single_slice mode needs a 3D (Z,Y,X) stack; "
                f"{source_path.name} has shape {image.shape}."
            )
        focus_slice = _select_focus_slice(image)
        return mask[focus_slice]

    if mode == "2d":
        return mask if not is_stack else mask.max(axis=0)

    if mode == "3d":
        if not is_stack:
            raise ValueError(
                f"3d mode needs a 3D (Z,Y,X) stack; "
                f"{source_path.name} has shape {image.shape}."
            )
        return mask

    raise ValueError(
        f"Unsupported processing mode {mode!r}; expected '3d', '2d', or 'single_slice'."
    )


def collect_size_preview(
    input_folder,
    mode="3d",
    expansion_factor=1.0,
    pixel_size_nm=None,
    z_step_nm=None,
    min_voxels=5,
    signal_channel=1,
):

    folder = Path(input_folder)
    if not folder.is_dir():
        raise ValueError(f"Input folder does not exist or is not a directory: {folder}")

    if mode not in {"3d", "2d", "single_slice"}:
        raise ValueError(
            f"Unsupported processing mode {mode!r}; expected '3d', '2d', or 'single_slice'."
        )

    pixel_size_nm = _positive_finite(pixel_size_nm, "Pixel size XY")
    expansion_factor = _positive_finite(expansion_factor, "Expansion factor")
    if mode == "3d":
        z_step_nm = _positive_finite(z_step_nm, "Z-step")

    try:
        min_voxels = int(min_voxels)
    except (TypeError, ValueError) as error:
        raise ValueError("Raw minimum pixels/voxels must be a positive integer.") from error
    if min_voxels < 1:
        raise ValueError("Raw minimum pixels/voxels must be at least 1.")

    try:
        signal_channel = int(signal_channel)
    except (TypeError, ValueError) as error:
        raise ValueError("Signal channel must be a non-negative integer.") from error
    if signal_channel < 0:
        raise ValueError("Signal channel must be a non-negative integer.")

    source_files = sorted(
        path
        for path in folder.glob("*.tif")
        if "Mask" not in path.name and "Final" not in path.name
    )
    if not source_files:
        return _empty_preview_table()

    rows = []
    effective_xy_nm = pixel_size_nm / expansion_factor
    if mode == "3d":
        effective_z_nm = z_step_nm / expansion_factor
        size_per_element = (effective_xy_nm**2 * effective_z_nm) / 1e9
        unit = "µm³"
    else:
        size_per_element = effective_xy_nm**2 / 1e6
        unit = "µm²"

    for source_path in source_files:
        mask_path = folder / f"{source_path.stem}_Mask.tif"
        if not mask_path.is_file():
            raise ValueError(f"Missing matching mask for {source_path.name}: {mask_path.name}")

        processing_mask = _load_processing_mask(
            source_path,
            mask_path,
            mode,
            signal_channel,
        )
        labeled_mask = _prepare_labeled_mask(processing_mask)

        for region in regionprops(labeled_mask):
            if region.area < min_voxels:
                continue
            rows.append(
                {
                    "filename": source_path.name,
                    "mask_filename": mask_path.name,
                    "mode": mode,
                    "object_id": int(region.label),
                    "raw_size": int(region.area),
                    "biological_size": round(region.area * size_per_element, 5),
                    "unit": unit,
                }
            )

    if not rows:
        return _empty_preview_table()
    return pd.DataFrame.from_records(rows, columns=PREVIEW_COLUMNS)
