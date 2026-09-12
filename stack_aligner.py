"""Small XY-drift aligner for ExQt Z-stacks."""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.registration import phase_cross_correlation
import tifffile


@dataclass(frozen=True)
class AlignmentConfig:
    upsample_factor: int = 10
    max_step_px: float = 8.0
    min_texture_std: float = 1e-6
    highpass_sigma: float = 8.0
    smooth_sigma: float = 1.0
    edge_taper: bool = True
    min_post_correlation: float = 0.20
    max_correlation_loss: float = 0.05
    max_bidirectional_disagreement_px: float = 0.75
    max_fail_fraction: float = 0.20
    expand_canvas: bool = True
    min_residual_correlation: float = 0.20




@dataclass(frozen=True)
class AlignmentResult:
    step_shifts_yx: np.ndarray
    cumulative_shifts_yx: np.ndarray
    metrics: pd.DataFrame
    status: str
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class AlignedFile:
    source: Path
    output: Path | None
    mask_output: Path | None
    drift_plot: Path
    drift_csv: Path
    alignment: AlignmentResult


@dataclass(frozen=True)
class AlignmentBatch:
    records: pd.DataFrame
    summary_csv: Path


def _prepare_plane(plane: np.ndarray, config: AlignmentConfig) -> np.ndarray | None:
    image = np.asarray(plane, dtype=np.float64)
    finite = np.isfinite(image)
    if not finite.any():
        return None
    image = np.where(finite, image, np.median(image[finite]))
    low, high = np.percentile(image, (1, 99))
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return None
    image = np.clip(image, low, high)
    if config.highpass_sigma > 0:
        image = image - ndi.gaussian_filter(image, config.highpass_sigma)
    if config.smooth_sigma > 0:
        image = ndi.gaussian_filter(image, config.smooth_sigma)
    std = float(np.std(image))
    if not np.isfinite(std) or std < config.min_texture_std:
        return None
    image = (image - float(np.mean(image))) / std
    if config.edge_taper and min(image.shape) >= 4:
        image = image * np.outer(np.hanning(image.shape[0]), np.hanning(image.shape[1]))
    return image


def _correlation(reference: np.ndarray, moving: np.ndarray, shift_yx: np.ndarray) -> float:
    shifted = ndi.shift(
        moving, tuple(float(value) for value in shift_yx), order=1,
        mode="constant", cval=np.nan, prefilter=False,
    )
    valid = np.isfinite(reference) & np.isfinite(shifted)
    if np.count_nonzero(valid) < 16:
        return np.nan
    a = reference[valid].astype(float)
    b = shifted[valid].astype(float)
    a -= a.mean()
    b -= b.mean()
    denominator = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / denominator) if denominator > 0 else np.nan


def _notify_progress(callback: Callable | None, pct: int, text: str) -> None:
    if callback is None:
        return
    try:
        callback(pct, text)
    except TypeError:
        callback(f"[{pct}%] {text}")


def estimate_xy_drift(
    reference_stack: np.ndarray,
    config: AlignmentConfig | None = None,
    progress_callback: Callable | None = None,
) -> AlignmentResult:
    """Estimate one Y/X translation per adjacent Z pair in a (Z,Y,X) stack."""
    config = config or AlignmentConfig()
    stack = np.asarray(reference_stack)
    if stack.ndim != 3 or stack.shape[0] < 2:
        raise ValueError("reference_stack must have shape (Z,Y,X) with at least two planes")
    if config.upsample_factor < 1 or config.max_step_px <= 0:
        raise ValueError("upsample_factor must be >= 1 and max_step_px must be positive")
    if not 0 <= config.min_post_correlation <= 1:
        raise ValueError("min_post_correlation must be between 0 and 1")
    if not 0 <= config.max_correlation_loss <= 1:
        raise ValueError("max_correlation_loss must be between 0 and 1")
    if config.max_bidirectional_disagreement_px <= 0:
        raise ValueError("max_bidirectional_disagreement_px must be positive")
    if not 0.0 < config.max_fail_fraction <= 1.0:
        raise ValueError("max_fail_fraction must be between 0 (exclusive) and 1 (inclusive)")

    shifts = np.zeros((stack.shape[0], 2), dtype=float)
    cumulative = np.zeros_like(shifts)
    records, reasons = [], []
    saw_review = False
    prepared = [_prepare_plane(plane, config) for plane in stack]
    total_steps = stack.shape[0] - 1

    for z in range(1, stack.shape[0]):
        if progress_callback and (z % 2 == 0 or z == total_steps or z == 1):
            plane_pct = int((z / total_steps) * 100)
            _notify_progress(progress_callback, plane_pct, f"Calculating drift: slice {z}/{total_steps}")
        previous, current = prepared[z - 1], prepared[z]
        status, reason, phase_error = "PASS", "", np.nan
        before = after = np.nan
        shift = np.zeros(2, dtype=float)
        reverse_shift = np.full(2, np.nan, dtype=float)
        bidirectional_disagreement = np.nan
        if previous is None or current is None:
            status, reason = "FAIL", "insufficient_texture"
        else:
            before = _correlation(previous, current, np.zeros(2))
            try:
                estimated, phase_error, _ = phase_cross_correlation(
                    previous, current, upsample_factor=config.upsample_factor,
                    normalization="phase",
                )
                shift = np.asarray(estimated[-2:], dtype=float)
                reverse_estimated, _, _ = phase_cross_correlation(
                    current, previous, upsample_factor=config.upsample_factor,
                    normalization="phase",
                )
                reverse_shift = np.asarray(reverse_estimated[-2:], dtype=float)
                bidirectional_disagreement = float(np.linalg.norm(shift + reverse_shift))
                after = _correlation(previous, current, shift)
            except (ValueError, FloatingPointError, RuntimeError) as error:
                status, reason = "FAIL", f"registration_error:{type(error).__name__}"

            if status == "PASS" and not np.all(np.isfinite(shift)):
                status, reason = "FAIL", "non_finite_shift"
            elif status == "PASS" and np.linalg.norm(shift) > config.max_step_px:
                status, reason = "FAIL", "step_above_limit"
            elif status == "PASS" and not (isinstance(after, (int, float, np.floating)) and np.isfinite(after)):
                status, reason = "FAIL", "non_finite_correlation"
            elif status == "PASS" and not np.all(np.isfinite(reverse_shift)):
                status, reason = "REVIEW", "non_finite_reverse_shift"
            elif status == "PASS" and bidirectional_disagreement > config.max_bidirectional_disagreement_px:
                status, reason = "REVIEW", "inconsistent_bidirectional_shift"
            elif status == "PASS" and after < config.min_post_correlation:
                status, reason = "REVIEW", "weak_post_correlation"
            elif status == "PASS" and isinstance(before, (int, float, np.floating)) and np.isfinite(before) and after < before - config.max_correlation_loss:
                status, reason = "REVIEW", "correlation_loss"
            elif status == "PASS" and isinstance(after, (int, float, np.floating)) and np.isfinite(after) and after < config.min_residual_correlation:

                status, reason = "REVIEW", "residual_misalignment_after_translation"

        if status == "REVIEW":
            saw_review = True
        if status == "FAIL":
            #Freeze shift at zero for this step; whether this makes the overallresult a hard FAIL is decided after all steps based on fail fraction.
            shift = np.zeros(2, dtype=float)
        shifts[z] = shift
        cumulative[z] = cumulative[z - 1] + shift
        if reason:
            reasons.append(f"Z{z - 1}->Z{z}: {reason}")
        records.append({
            "reference_z": z - 1, "moving_z": z,
            "step_shift_y_px": float(shift[0]), "step_shift_x_px": float(shift[1]),
            "step_magnitude_px": float(np.linalg.norm(shift)),
            "cumulative_shift_y_px": float(cumulative[z, 0]),
            "cumulative_shift_x_px": float(cumulative[z, 1]),
            "correlation_before": before, "correlation_after": after,
            "phase_error": float(phase_error), "reverse_shift_y_px": float(reverse_shift[0]),
            "reverse_shift_x_px": float(reverse_shift[1]), "bidirectional_disagreement_px": bidirectional_disagreement,
            "status": status, "reason": reason,
        })

    statuses = [r["status"] for r in records]
    reasons_list = [r["reason"] for r in records]

    #Find the contiguous runs of insufficient_texture at both ends.
    first_non_texture = next(
        (i for i, (s, rr) in enumerate(zip(statuses, reasons_list))
         if not (s == "FAIL" and rr == "insufficient_texture")),
        len(statuses),
    )
    last_non_texture = next(
        (len(statuses) - 1 - i
         for i, (s, rr) in enumerate(zip(reversed(statuses), reversed(reasons_list)))
         if not (s == "FAIL" and rr == "insufficient_texture")),
        -1,
    )

    if last_non_texture < first_non_texture:
        #Every step is insufficient_texture — no sample signal at all.
        overall = "FAIL"
    else:
        mid_records = records[first_non_texture: last_non_texture + 1]
        mid_fails = sum(1 for r in mid_records if r["status"] == "FAIL")
        mid_steps = len(mid_records)
        if mid_fails == 0:
            overall = "REVIEW" if saw_review else "PASS"
        elif mid_fails / mid_steps <= config.max_fail_fraction:
            overall = "REVIEW"
        else:
            overall = "FAIL"

    return AlignmentResult(shifts, cumulative, pd.DataFrame(records), overall, tuple(reasons))


def extract_reference_stack(data: np.ndarray, axes: str, channel_index: int = 0) -> np.ndarray:
    """Return the selected channel in unambiguous (Z,Y,X) order."""
    array, axes = np.asarray(data), axes.upper()
    if len(axes) != array.ndim or len(set(axes)) != len(axes):
        raise ValueError(f"Axes {axes!r} do not describe array shape {array.shape}")
    if set(axes) - {"T", "C", "Z", "Y", "X"} or not set("ZYX") <= set(axes):
        raise ValueError("axes must contain only T/C/Z/Y/X and include Z/Y/X")
    if "T" in axes:
        index = axes.index("T")
        if array.shape[index] != 1:
            raise ValueError("Time series are not supported by this light aligner")
        array, axes = np.take(array, 0, axis=index), axes.replace("T", "")
    if "C" in axes:
        index = axes.index("C")
        if not 0 <= channel_index < array.shape[index]:
            raise IndexError(f"channel_index {channel_index} is outside the available channels")
        array, axes = np.take(array, channel_index, axis=index), axes.replace("C", "")
    elif channel_index != 0:
        raise ValueError("channel_index must be zero when no C axis is present")
    return np.transpose(array, [axes.index(axis) for axis in "ZYX"])


def apply_xy_shifts(
    data: np.ndarray,
    axes: str,
    cumulative_shifts_yx: np.ndarray,
    interpolation_order: int = 1,
    expand_canvas: bool = False,
) -> np.ndarray:
    """Apply per-Z Y/X shifts while retaining the original axis order and dtype."""
    array, axes = np.asarray(data), axes.upper()
    if len(axes) != array.ndim or len(set(axes)) != len(axes):
        raise ValueError(f"axes {axes!r} do not describe array shape {array.shape}")
    if set(axes) - {"T", "C", "Z", "Y", "X"} or not set("ZYX") <= set(axes):
        raise ValueError(f"axes {axes!r} must contain only T/C/Z/Y/X and include Z/Y/X")
    shifts = np.asarray(cumulative_shifts_yx, dtype=float)
    z_axis = axes.index("Z")
    if shifts.shape != (array.shape[z_axis], 2) or not np.all(np.isfinite(shifts)):
        raise ValueError("cumulative_shifts_yx must be finite and have shape (Z,2)")
    if interpolation_order not in (0, 1, 3):
        raise ValueError("interpolation_order must be 0, 1, or 3")

    working = np.moveaxis(array, z_axis, 0)
    y_axis = axes.replace("Z", "").index("Y")
    x_axis = axes.replace("Z", "").index("X")

    if expand_canvas:
        min_y, max_y = min(0.0, float(shifts[:, 0].min())), max(0.0, float(shifts[:, 0].max()))
        min_x, max_x = min(0.0, float(shifts[:, 1].min())), max(0.0, float(shifts[:, 1].max()))
        pad_top = int(np.ceil(abs(min_y))) if min_y < 0 else 0
        pad_bottom = int(np.ceil(max_y)) if max_y > 0 else 0
        pad_left = int(np.ceil(abs(min_x))) if min_x < 0 else 0
        pad_right = int(np.ceil(max_x)) if max_x > 0 else 0

        plane_shape = list(working[0].shape)
        plane_shape[y_axis] += (pad_top + pad_bottom)
        plane_shape[x_axis] += (pad_left + pad_right)
        output = np.zeros((working.shape[0], *plane_shape), dtype=working.dtype)

        in_slices = [slice(None)] * working[0].ndim
        in_slices[y_axis] = slice(pad_top, pad_top + working[0].shape[y_axis])
        in_slices[x_axis] = slice(pad_left, pad_left + working[0].shape[x_axis])

        for z, shift in enumerate(shifts):
            padded_plane = np.zeros(tuple(plane_shape), dtype=working.dtype)
            padded_plane[tuple(in_slices)] = working[z]
            full_shift = np.zeros(working[z].ndim, dtype=float)
            full_shift[y_axis] = shift[0]
            full_shift[x_axis] = shift[1]
            output[z] = ndi.shift(
                padded_plane, full_shift, order=interpolation_order,
                mode="constant", cval=0, prefilter=interpolation_order > 1,
            )
    else:
        output = np.empty_like(working)
        for z, shift in enumerate(shifts):
            full_shift = np.zeros(working[z].ndim, dtype=float)
            full_shift[y_axis] = shift[0]
            full_shift[x_axis] = shift[1]
            output[z] = ndi.shift(
                working[z], full_shift, order=interpolation_order,
                mode="constant", cval=0, prefilter=interpolation_order > 1,
            )
    return np.moveaxis(output, 0, z_axis)



def _read_tiff(path: Path) -> tuple[np.ndarray, str, str | None]:
    with tifffile.TiffFile(path) as tif:
        if len(tif.series) != 1:
            raise ValueError(f"{path.name}: expected one TIFF series")
        series = tif.series[0]
        data, axes, ome = series.asarray(), series.axes.upper(), tif.ome_metadata
    if axes == "IYX" and data.ndim == 3:
        axes = "ZYX"
    if not axes:
        if data.ndim == 3:
            #A plain 3-D stack has an unambiguous Z/Y/X interpretation.
            axes = "ZYX"
        else:
            raise ValueError(f"{path.name}: missing TIFF axes; save it as OME-TIFF with explicit axes")
    elif set(axes) - {"T", "C", "Z", "Y", "X"}:
        raise ValueError(f"{path.name}: unsupported TIFF axes {axes!r}; expected T/C/Z/Y/X")
    if not np.issubdtype(data.dtype, np.number):
        raise ValueError(f"{path.name}: TIFF data must be numeric")
    if not np.isfinite(data).all():
        raise ValueError(f"{path.name}: TIFF contains NaN or infinite values")
    return data, axes, ome


def _ome_physical_metadata(ome: str | None) -> dict[str, object]:
    if not ome:
        return {}
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(ome)
        pixels = next((node for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "Pixels"), None)
    except ET.ParseError:
        return {}
    if pixels is None:
        return {}
    metadata = {key: pixels.attrib[key] for key in (
        "PhysicalSizeX", "PhysicalSizeXUnit", "PhysicalSizeY", "PhysicalSizeYUnit",
        "PhysicalSizeZ", "PhysicalSizeZUnit", "TimeIncrement", "TimeIncrementUnit",
    ) if key in pixels.attrib}
    return metadata


def _extract_physical_metadata(source_path: Path | None = None, ome: str | None = None) -> dict[str, object]:
    """Preserve physical calibration from OME XML or fall back to ImageJ / TIFF tags."""
    metadata = _ome_physical_metadata(ome) if ome else {}
    if "PhysicalSizeX" in metadata and "PhysicalSizeZ" in metadata:
        return metadata
    if source_path is not None:
        try:
            from Batch import get_metadata_from_tif
            batch_meta = get_metadata_from_tif(source_path)
            if batch_meta:
                if "pixel_size" in batch_meta and "PhysicalSizeX" not in metadata:
                    px_um = float(batch_meta["pixel_size"]) / 1000.0
                    metadata["PhysicalSizeX"] = px_um
                    metadata["PhysicalSizeXUnit"] = "µm"
                    metadata["PhysicalSizeY"] = px_um
                    metadata["PhysicalSizeYUnit"] = "µm"
                if "z_step" in batch_meta and "PhysicalSizeZ" not in metadata:
                    z_um = float(batch_meta["z_step"]) / 1000.0
                    metadata["PhysicalSizeZ"] = z_um
                    metadata["PhysicalSizeZUnit"] = "µm"
        except Exception:
            pass
    return metadata


def _canonicalize_ome_axes(data: np.ndarray, axes: str) -> tuple[np.ndarray, str]:
    """Ensure array axes end with YX as required by OME-TIFF writers."""
    if axes.endswith("YX"):
        return data, axes
    axis_set = set(axes)
    if axis_set == {"Z", "Y", "X"}:
        target = "ZYX"
    elif axis_set == {"C", "Z", "Y", "X"}:
        target = "CZYX"
    elif axis_set == {"T", "C", "Z", "Y", "X"}:
        target = "TCZYX"
    elif axis_set == {"Y", "X"}:
        target = "YX"
    else:
        return data, axes
    perm = [axes.index(a) for a in target]
    return np.transpose(data, perm), target


def _atomic_write_tiff(
    path: Path,
    data: np.ndarray,
    axes: str,
    ome: str | None = None,
    source_path: Path | None = None,
) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    data_to_write, write_axes = _canonicalize_ome_axes(data, axes)
    physical_meta = _extract_physical_metadata(source_path=source_path, ome=ome)
    metadata = {"axes": write_axes, **physical_meta}
    try:
        tifffile.imwrite(temporary, data_to_write, ome=True, metadata=metadata)
        written, written_axes, _ = _read_tiff(temporary)
        if written.shape != data_to_write.shape or written_axes != write_axes:
            raise IOError(
                f"TIFF verification failed for {path.name}: expected {write_axes}/{data_to_write.shape}, "
                f"got {written_axes}/{written.shape}"
            )
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _plot_drift(path: Path, alignment: AlignmentResult) -> None:
    z = np.arange(len(alignment.cumulative_shifts_yx))
    figure, (cumulative_axis, step_axis) = plt.subplots(2, 1, figsize=(7, 6), sharex=True, constrained_layout=True)
    cumulative_axis.plot(z, alignment.cumulative_shifts_yx[:, 0], "o-", label="Y cumulative")
    cumulative_axis.plot(z, alignment.cumulative_shifts_yx[:, 1], "o-", label="X cumulative")
    cumulative_axis.axhline(0, color="0.6", linewidth=0.8)
    cumulative_axis.set_ylabel("Cumulative correction (px)")
    cumulative_axis.set_title(f"XY drift — {alignment.status}")
    cumulative_axis.grid(alpha=0.25)
    cumulative_axis.legend()
    step_z = np.arange(1, len(alignment.step_shifts_yx))
    step_axis.plot(step_z, alignment.step_shifts_yx[1:, 0], ".-", label="Y step")
    step_axis.plot(step_z, alignment.step_shifts_yx[1:, 1], ".-", label="X step")
    step_axis.axhline(0, color="0.6", linewidth=0.8)
    step_axis.set(xlabel="Moving Z plane", ylabel="Step correction (px)")
    step_axis.grid(alpha=0.25)
    step_axis.legend()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def _drift_stem(name: str) -> str:
    for suffix in (".ome.tiff", ".ome.tif", ".tiff", ".tif"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def align_tiff_stack(
    source_path,
    output_folder,
    reference_channel=0,
    mask_path=None,
    align_mask=True,
    config=None,
    progress_callback: Callable | None = None,
) -> AlignedFile:
    """Align one TIFF; write QC always, write images for PASS and REVIEW, block only on FAIL."""
    config = config or AlignmentConfig()
    source = Path(source_path)
    output_folder = Path(output_folder)
    if source.resolve().parent == output_folder.resolve():
        raise ValueError("output_folder must be different from the source folder")
    data, axes, ome = _read_tiff(source)
    reference = extract_reference_stack(data, axes, reference_channel)
    alignment = estimate_xy_drift(reference, config, progress_callback=progress_callback)


    mask = mask_axes = mask_ome = None
    if mask_path and align_mask:
        mask, mask_axes, mask_ome = _read_tiff(Path(mask_path))
        if mask.dtype.kind not in "biu" or np.any(mask < 0):
            raise ValueError(f"{Path(mask_path).name}: mask must contain non-negative integer labels")
        mask_reference = extract_reference_stack(mask, mask_axes, 0)
        if mask_reference.shape != reference.shape:
            raise ValueError(f"{Path(mask_path).name}: mask Z/Y/X shape does not match the source")

    stem = _drift_stem(source.name)
    drift_csv = output_folder / f"{stem}_drift.csv"
    drift_plot = output_folder / f"{stem}_drift.png"
    output = output_folder / source.name if alignment.status in ("PASS", "REVIEW") else None
    mask_output = output_folder / Path(mask_path).name if output and mask is not None else None
    for existing in (output, mask_output, drift_csv, drift_plot):
        if existing is not None and existing.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {existing}")

    if output is not None:
        try:
            _atomic_write_tiff(
                output,
                apply_xy_shifts(data, axes, alignment.cumulative_shifts_yx, expand_canvas=config.expand_canvas),
                axes,
                ome,
                source_path=source,
            )
            if mask_output is not None:
                _atomic_write_tiff(
                    mask_output,
                    apply_xy_shifts(mask, mask_axes, alignment.cumulative_shifts_yx, interpolation_order=0, expand_canvas=config.expand_canvas),
                    mask_axes,
                    mask_ome,
                    source_path=Path(mask_path),
                )
        except Exception:
            for partial in (output, mask_output):
                if partial is not None and partial.exists():
                    partial.unlink()
            raise
    output_folder.mkdir(parents=True, exist_ok=True)
    alignment.metrics.to_csv(drift_csv, index=False)
    _plot_drift(drift_plot, alignment)
    if progress_callback:
        action = "written" if output is not None else "blocked"
        _notify_progress(progress_callback, 100, f"{source.name}: {alignment.status} ({action})")
    return AlignedFile(source, output, mask_output, drift_plot, drift_csv, alignment)


def align_tiff_folder(
    input_folder,
    output_folder,
    reference_channel=0,
    align_masks=True,
    config=None,
    progress_callback: Callable | None = None,
) -> AlignmentBatch:
    """Align all source TIFFs in a folder; matching *_Mask files follow automatically."""
    input_folder, output_folder = Path(input_folder), Path(output_folder)
    if input_folder.resolve() == output_folder.resolve():
        raise ValueError("output_folder must be different from input_folder")
    if output_folder.exists() and any(output_folder.iterdir()):
        raise FileExistsError(f"Output folder is not empty; choose a new folder: {output_folder}")
    sources = sorted(
        path for path in input_folder.iterdir()
        if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}
        and not path.stem.lower().endswith("_mask")
    )
    if not sources:
        raise FileNotFoundError(f"No TIFF files found in {input_folder}")
    rows = []
    n_sources = len(sources)
    for i, source in enumerate(sources):
        def stack_progress(plane_pct: int, plane_text: str) -> None:
            overall_pct = int(((i + (plane_pct / 100.0)) / n_sources) * 100)
            _notify_progress(progress_callback, overall_pct, f"[{i + 1}/{n_sources}] {plane_text}")

        mask = next((source.with_name(f"{source.stem}_Mask{suffix}") for suffix in (".tif", ".tiff") if source.with_name(f"{source.stem}_Mask{suffix}").is_file()), None)
        try:
            result = align_tiff_stack(source, output_folder, reference_channel, mask, align_masks, config, progress_callback=stack_progress)
            rows.append({"source_file": source.name, "status": result.alignment.status, "output_file": result.output.name if result.output else "", "mask_output": result.mask_output.name if result.mask_output else "", "drift_plot": result.drift_plot.name, "drift_csv": result.drift_csv.name, "reason": "; ".join(result.alignment.reasons)})
        except Exception as error:
            rows.append({"source_file": source.name, "status": "ERROR", "output_file": "", "mask_output": "", "drift_plot": "", "drift_csv": "", "reason": str(error)})
            if progress_callback:
                _notify_progress(progress_callback, int(((i + 1) / n_sources) * 100), f"Error in {source.name}: {error}")
    records = pd.DataFrame(rows)
    summary = output_folder / "alignment_summary.csv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    records.to_csv(summary, index=False)
    if progress_callback:
        _notify_progress(progress_callback, 100, "Alignment batch complete.")
    return AlignmentBatch(records, summary)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Light XY drift aligner for ExQt TIFF Z-stacks")
    parser.add_argument("input_folder")
    parser.add_argument("output_folder")
    parser.add_argument("--reference-channel", type=int, default=0)
    parser.add_argument("--no-masks", action="store_true")
    args = parser.parse_args()
    batch = align_tiff_folder(args.input_folder, args.output_folder, args.reference_channel, not args.no_masks, progress_callback=print)
    print(f"Saved summary: {batch.summary_csv}")
