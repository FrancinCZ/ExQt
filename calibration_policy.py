
from __future__ import annotations


def summarize_calibrations(metadata_by_file: dict) -> dict:
    #Report whether detected TIFF calibrations are consistent within one batch
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
