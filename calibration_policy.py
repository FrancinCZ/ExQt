# Last significant change: 2026-09-23 (explicit calibration source, plausibility limits)
from __future__ import annotations

import math


# Same limits as the Advanced Settings spinboxes; anything outside is an error, not a value.
CALIBRATION_LIMITS_NM = {
    "pixel_size_nm": (1.0, 2000.0),
    "z_step_nm": (1.0, 5000.0),
}
CALIBRATION_SOURCE_GUI = "gui"
CALIBRATION_SOURCE_TIFF = "tiff_per_file"
CALIBRATION_SOURCES = (CALIBRATION_SOURCE_GUI, CALIBRATION_SOURCE_TIFF)


#Report whether detected TIFF calibrations are consistent within one batch.
def summarize_calibrations(metadata_by_file: dict) -> dict:
    pairs = {
        (
            round(float(values["pixel_size_nm"]), 6),
            round(float(values["z_step_nm"]), 6),
        )
        for values in metadata_by_file.values()
        if "pixel_size_nm" in values and "z_step_nm" in values
    }
    return {
        "calibration_count": len(pairs),
        "calibrations": sorted(pairs),
        "has_mismatch": len(pairs) > 1,
    }


#Raise ValueError when a calibration lies outside the supported acquisition range.
def validate_calibration(pixel_size_nm, z_step_nm) -> None:
    for name, value in (("pixel_size_nm", pixel_size_nm), ("z_step_nm", z_step_nm)):
        low, high = CALIBRATION_LIMITS_NM[name]
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Calibration {name}={value!r} is not a number.") from None
        if not math.isfinite(number) or not low <= number <= high:
            raise ValueError(
                f"Calibration {name}={number:g} nm is outside the plausible range {low:g}–{high:g} nm."
            )


#Return the calibration applied to one file and where it came from.
def resolve_file_calibration(
    filename, calibration_source, gui_pixel_size_nm, gui_z_step_nm, detected_metadata_by_file
) -> dict:
    if calibration_source not in CALIBRATION_SOURCES:
        raise ValueError(f"Unknown calibration_source {calibration_source!r}; expected one of {CALIBRATION_SOURCES}.")

    applied = {
        "pixel_size_nm": float(gui_pixel_size_nm),
        "z_step_nm": float(gui_z_step_nm),
        "source": "gui",
    }
    if calibration_source == CALIBRATION_SOURCE_TIFF:
        detected = (detected_metadata_by_file or {}).get(filename) or {}
        if "pixel_size_nm" in detected and "z_step_nm" in detected:
            applied = {
                "pixel_size_nm": float(detected["pixel_size_nm"]),
                "z_step_nm": float(detected["z_step_nm"]),
                "source": "tiff",
            }
        else:
            applied["source"] = "gui_fallback_no_tiff_metadata"

    validate_calibration(applied["pixel_size_nm"], applied["z_step_nm"])
    return applied


#True when TIFF calibration ≈ GUI calibration / ExF, i.e. the file was likely calibrated
#to pre-expansion size and applying it would divide by ExF twice.
def looks_pre_expansion(
    pixel_size_nm, z_step_nm, gui_pixel_size_nm, gui_z_step_nm, expansion_factor, tolerance=0.10
) -> bool:
    if expansion_factor is None or float(expansion_factor) <= 1.0 + tolerance:
        return False
    expected_xy = float(gui_pixel_size_nm) / float(expansion_factor)
    expected_z = float(gui_z_step_nm) / float(expansion_factor)
    return (
        abs(float(pixel_size_nm) - expected_xy) <= tolerance * expected_xy
        and abs(float(z_step_nm) - expected_z) <= tolerance * expected_z
    )
