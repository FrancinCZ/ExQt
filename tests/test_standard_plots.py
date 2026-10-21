# Standard plots (supervisor feedback 2026-10-17): every run is drawn with the same panels and
# fixed axes so two runs can be compared by eye, and intensity is shown as K_part (nucleoplasm = 1)
# instead of raw ADU, which is not comparable between images.
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from defaults import DEFAULT_SETTINGS
from postprocessing import build_standard_figure, standard_plot_frame

MIN_SIZE, MAX_SIZE = 0.0001, 2.0


def _frame(scale=1.0, n=12):
    rng = np.random.default_rng(0)
    rows = []
    for i in range(n):
        rows.append(dict(
            filename="a.tif" if i < n // 2 else "b.tif", cell_id=1, object_id=i + 1, mode="3d", is_3d=True,
            volume_px=600, volume_bio_um3=0.001 * scale * (i + 1),
            mean_intensity=1000.0 * scale, integrated_density=1e5 * scale,
            partition_coefficient=1.5 + 0.2 * i * scale, K_valid=True,
            A_object=float(rng.uniform(0.1, 0.6)),
            touches_image_edge=False, touches_roi_edge=False,
        ))
    return pd.DataFrame(rows)


def test_frame_keeps_only_primary_objects_and_masks_invalid_k():
    df = _frame()
    df.loc[0, "touches_image_edge"] = True
    df.loc[1, "volume_bio_um3"] = 5.0
    df.loc[2, "K_valid"] = False
    plot_df = standard_plot_frame(df, MIN_SIZE, MAX_SIZE)
    assert len(plot_df) == len(df) - 2
    assert np.isnan(plot_df.loc[plot_df["object_id"] == 3, "partition_coefficient"]).all()


def test_no_raw_intensity_axes():
    fig = build_standard_figure(_frame(), MIN_SIZE, MAX_SIZE, title="run")
    labels = " ".join(ax.get_xlabel() + " " + ax.get_ylabel() for ax in fig.axes)
    assert "ADU" not in labels and "a.u." not in labels and "Mean Intensity" not in labels
    assert "K_part" in labels or "K$_{part}$" in labels


def test_axes_are_fixed_and_identical_for_different_runs():
    fig_a = build_standard_figure(_frame(scale=1.0), MIN_SIZE, MAX_SIZE, title="a")
    fig_b = build_standard_figure(_frame(scale=3.0), MIN_SIZE, MAX_SIZE, title="b")
    checked = {"size": 0, "k": 0, "fa": 0}
    for ax_a, ax_b in zip(fig_a.axes, fig_b.axes):
        assert ax_a.get_title() == ax_b.get_title()
        if ax_a.get_xscale() == "log" and "volume" in ax_a.get_xlabel().lower():
            assert np.allclose(ax_a.get_xlim(), (MIN_SIZE, MAX_SIZE))
            assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())
            checked["size"] += 1
        if "K" in ax_a.get_ylabel():
            assert np.allclose(ax_a.get_ylim(), (DEFAULT_SETTINGS["plot_k_min"], DEFAULT_SETTINGS["plot_k_max"]))
            assert np.allclose(ax_a.get_ylim(), ax_b.get_ylim())
            checked["k"] += 1
        if "FA" in ax_a.get_ylabel():
            assert np.allclose(ax_a.get_ylim(), (0.0, 1.0))
            checked["fa"] += 1
    assert checked == {"size": 3, "k": 2, "fa": 2}


def test_k_baseline_line_at_one():
    fig = build_standard_figure(_frame(), MIN_SIZE, MAX_SIZE, title="run")
    k_axes = [ax for ax in fig.axes if "K" in ax.get_ylabel() and "volume" in ax.get_xlabel().lower()]
    assert k_axes
    assert any(np.allclose(line.get_ydata(), 1.0) for line in k_axes[0].get_lines())


def test_title_reports_objects_and_nuclei():
    fig = build_standard_figure(_frame(), MIN_SIZE, MAX_SIZE, title="run")
    assert "12 objects" in fig._suptitle.get_text() and "2 nuclei" in fig._suptitle.get_text()


def test_fa_panels_explain_missing_fa():
    fig = build_standard_figure(_frame().drop(columns=["A_object"]), MIN_SIZE, MAX_SIZE, title="run")
    texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
    assert "Radial FA Profiling" in texts


def test_generate_plots_writes_png(tmp_path):
    from postprocessing import generate_plots

    csv = tmp_path / "run_Output_Batch_3d.csv"
    _frame().to_csv(csv, index=False)
    generate_plots(str(csv), min_size=MIN_SIZE, max_size=MAX_SIZE)
    assert (tmp_path / "run_Output_Batch_3d_Analysis_Plots.png").is_file()
