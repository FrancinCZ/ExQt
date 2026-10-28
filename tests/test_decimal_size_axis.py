# Supervisor feedback 2026-10-27: size axes stay logarithmic but show plain decimals (0.01, 0.1, 1)
# instead of powers of ten (10^-2, 10^-1, 10^0).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from postprocessing import build_condition_comparison_figure, build_standard_figure
from partitioning_plots import _four_panel_figure

MIN_SIZE, MAX_SIZE = 0.005, 2.0


def _frame():
    rows = [dict(filename="a.tif", cell_id=1 + i % 2, object_id=i + 1, mode="3d", is_3d=True, volume_px=600,
                 volume_bio_um3=0.01 * (i + 1), mean_intensity=1000.0, integrated_density=1e5,
                 partition_coefficient=1.2 + 0.1 * i, K_valid=True, A_object=0.3,
                 touches_image_edge=False, touches_roi_edge=False, primary_qc_valid=True, run_id="Control")
            for i in range(10)]
    return pd.DataFrame(rows)


#Tick labels of every logarithmic volume axis in the figure (x and y).
def _labels(fig):
    fig.canvas.draw()
    out = []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            if axis.get_scale() == "log" and "volume" in axis.get_label_text().lower():
                out.append([t.get_text() for t in axis.get_ticklabels() if t.get_text()])
    plt.close(fig)
    return out


def _assert_decimal(axes_labels, expected_axes, expected=("0.01", "0.1", "1")):
    assert len(axes_labels) == expected_axes
    for labels in axes_labels:
        assert not any("10^" in t or "^{" in t for t in labels), labels
        assert set(expected) <= set(labels), labels


def test_standard_figure_size_axes_are_decimal():
    fig = build_standard_figure(_frame(), MIN_SIZE, MAX_SIZE, title="run")
    _assert_decimal(_labels(fig), 3)


def test_condition_figure_size_axes_are_decimal():
    fig = build_condition_comparison_figure(_frame(), MIN_SIZE, MAX_SIZE)
    _assert_decimal(_labels(fig), 3)


def test_pdf_summary_size_axes_are_decimal():
    #Autoscaled to the points (0.01–0.1 µm³): under two decades, so 0.02 and 0.05 are labelled too.
    fig = _four_panel_figure(_frame(), "run")
    _assert_decimal(_labels(fig), 3, expected=("0.01", "0.02", "0.05", "0.1"))
