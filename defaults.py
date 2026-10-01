"""Default settings shared by the GUI, the batch analysis, the size preview and reporting.

Change a default here only (VALIDATION_PROTOCOL §4): a value written in two places once
silently overrode the calibration just by opening a dialog.
"""

DEFAULT_SETTINGS = {
    "adv_pixel_size": 58.0,
    "adv_z_step": 250.0,
    "adv_detector_offset": 0.0,
    "adv_signal_ch": 1,
    "adv_dapi_ch": 0,
    "raw_min_voxels": 5,
    "plot_min_size": 0.0001,
    "plot_max_size": 2.0,
    "mode_a_enabled": False,
    "mode_a_min_core_voxels": 20,
    "mode_a_exclude_split_slices": True,
    "mode_a_z_split_min_component_voxels": 20,
    "mode_a_z_split_min_component_fraction": 0.10,
    "mode_a_z_split_pass_fraction": 0.10,
    "mode_a_z_split_review_fraction": 0.30,
}
