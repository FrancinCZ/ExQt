"""Statistics with the nucleus (ROI) as the unit (feedback 2026-11-25).

Every metric is first computed per nucleus from its QC-passed objects; a sample is then described
by the nuclei: mean ± SD for means, median (IQR) for medians. Class shares are the mean over
nuclei of each nucleus's percentage, so every nucleus weighs the same however many objects it has.
A nucleus is (run, image, cell_id); a nucleus without QC-passed objects counts with 0 condensates.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from partitioning_plots import MAIN_CLASS_ORDER

SHAPE_CLASSES = ("Elongated (High FA)", "Globular (Low FA)")
PROFILE_CLASSES = ("Core-Enriched", "No Significant Gradient", "Shell-Enriched")
OTHER_CLASS = "No shape / profile class"
ALL_NUCLEI = "All nuclei"

#(label, column in the per-nucleus table, how the nuclei are summarised)
SUMMARY_METRICS = [
    ("Condensates per nucleus (ROI)", "n_condensates", "mean"),
    ("Mean partitioning (K_part)", "mean_K", "mean"),
    ("Median partitioning (K_part)", "median_K", "median"),
    ("Mean particle mean intensity [ADU]", "mean_intensity", "mean"),
    ("Mean nucleoplasm intensity [ADU]", "nucleoplasm_intensity", "mean"),
    ("Mean particle volume [µm³]", "mean_volume", "mean"),
    ("Median particle volume [µm³]", "median_volume", "median"),
    ("Mean FA (object)", "mean_FA", "mean"),
]


def _bool(series):
    if series.dtype == bool:
        return series
    return series.map(lambda v: v.strip().lower() in ("true", "1", "yes") if isinstance(v, str)
                      else (bool(v) if pd.notna(v) else False))


def _size_column(frame):
    return next((c for c in ("volume_bio_um3", "area_bio_um2", "shape_metric_bio") if c in frame.columns), None)


def class_share_columns():
    """Per-nucleus percentage columns, in table order."""
    shape = [f"Total {name}" for name in SHAPE_CLASSES]
    profile = [f"Total {name}" for name in PROFILE_CLASSES]
    return shape + profile + list(MAIN_CLASS_ORDER) + [OTHER_CLASS]


def _class_shares(classes):
    classes = classes.astype(str)
    n = len(classes)
    shares = {}
    for name in SHAPE_CLASSES:
        shares[f"Total {name}"] = classes.str.startswith(name).sum()
    for name in PROFILE_CLASSES:
        shares[f"Total {name}"] = classes.str.endswith(name).sum()
    for name in MAIN_CLASS_ORDER:
        shares[name] = (classes == name).sum()
    shares[OTHER_CLASS] = (~classes.isin(MAIN_CLASS_ORDER)).sum()
    return {key: (100.0 * value / n if n else np.nan) for key, value in shares.items()}


def per_nucleus_table(report_df, group_cols=()):
    """One row per nucleus from all objects of a run (rows with primary_qc_valid)."""
    group_cols = [c for c in group_cols if c in report_df.columns]
    key_cols = group_cols + [c for c in ("run_id", "filename", "cell_id") if c in report_df.columns and c not in group_cols]
    frame = report_df.copy()
    primary_mask = _bool(frame["primary_qc_valid"]) if "primary_qc_valid" in frame.columns else pd.Series(True, index=frame.index)
    size_col = _size_column(frame)
    frame["_size"] = pd.to_numeric(frame[size_col], errors="coerce") if size_col else np.nan
    frame["_k"] = pd.to_numeric(frame.get("partition_coefficient"), errors="coerce")
    if "K_valid" in frame.columns:
        frame.loc[~_bool(frame["K_valid"]), "_k"] = np.nan
    frame["_fa"] = pd.to_numeric(frame.get("A_object"), errors="coerce")
    frame["_intensity"] = pd.to_numeric(frame.get("mean_intensity"), errors="coerce")
    frame["_nucleoplasm"] = pd.to_numeric(frame.get("nucleoplasm_mean_intensity"), errors="coerce")
    class_col = next((c for c in ("particle_class", "condensate_class") if c in frame.columns), None)

    rows = []
    for key, rows_of in frame.groupby(key_cols, sort=True, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        primary = rows_of[primary_mask.loc[rows_of.index]]
        row = dict(zip(key_cols, key))
        row.update({
            "n_condensates": len(primary),
            "mean_K": primary["_k"].mean(),
            "median_K": primary["_k"].median(),
            "mean_intensity": primary["_intensity"].mean(),
            "nucleoplasm_intensity": rows_of["_nucleoplasm"].mean(),
            "mean_volume": primary["_size"].mean(),
            "median_volume": primary["_size"].median(),
            "mean_FA": primary["_fa"].mean(),
        })
        row.update(_class_shares(primary[class_col]) if class_col else {})
        rows.append(row)
    return pd.DataFrame(rows)


def _describe(values, how):
    values = pd.to_numeric(values, errors="coerce").dropna()
    if values.empty:
        return {"n": 0, "mean": np.nan, "sd": np.nan, "median": np.nan, "q1": np.nan, "q3": np.nan}
    return {
        "n": len(values),
        "mean": values.mean(),
        "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
        "median": values.median(),
        "q1": values.quantile(0.25),
        "q3": values.quantile(0.75),
    }


def _fmt(value, digits):
    return "–" if not np.isfinite(value) else f"{value:.{digits}f}"


def _digits(label, column):
    if column == "n_condensates":
        return 1
    if "volume" in column:
        return 4
    if "intensity" in column:
        return 1
    if column.endswith("FA"):
        return 2
    return 2


def _text(stats, how, digits, percent=False):
    unit = "%" if percent else ""
    if stats["n"] == 0:
        return "–"
    if how == "median":
        return f"{_fmt(stats['median'], digits)}{unit} ({_fmt(stats['q1'], digits)}–{_fmt(stats['q3'], digits)})"
    if stats["n"] == 1:
        return f"{_fmt(stats['mean'], digits)}{unit}"
    return f"{_fmt(stats['mean'], digits)} ± {_fmt(stats['sd'], digits)}{unit}"


def summary_tables(nuclei, column_col=None, columns=None):
    """Side-by-side summary: one column per value of column_col (or per nucleus when column_col is
    None and columns=='nuclei'), rows = metrics. Returns (text table, numeric long table)."""
    if column_col is None:
        groups = [(ALL_NUCLEI, nuclei)]
    else:
        order = columns or list(dict.fromkeys(nuclei[column_col]))
        groups = [(name, nuclei[nuclei[column_col] == name]) for name in order]

    text = {}
    numbers = []
    metric_rows = [("Nuclei (ROIs)", None, "count")] + SUMMARY_METRICS + [
        (name, name, "mean") for name in class_share_columns() if name in nuclei.columns
    ]
    for group_name, group in groups:
        column = {}
        for label, source, how in metric_rows:
            if how == "count":
                column[label] = str(len(group))
                numbers.append({"column": group_name, "metric": label, "n": len(group)})
                continue
            if source not in group.columns:
                continue
            stats = _describe(group[source], how)
            percent = source in class_share_columns()
            column[label] = _text(stats, how, 1 if percent else _digits(label, source), percent)
            numbers.append({"column": group_name, "metric": label, "summary": "median (IQR)" if how == "median" else "mean ± SD", **stats})
        text[group_name] = column
    table = pd.DataFrame(text)
    table.index.name = "Metric"
    return table, pd.DataFrame(numbers)


def nucleus_columns_table(nuclei, label_cols=("filename", "cell_id")):
    """Each nucleus as its own column (its values, no spread), followed by 'All nuclei'."""
    columns = {}
    for _, row in nuclei.iterrows():
        name = " / ".join(
            (f"Cell {row[c]}" if c == "cell_id" else str(row[c]).rsplit(".", 1)[0])
            for c in label_cols if c in nuclei.columns
        )
        column = {"Nuclei (ROIs)": "1"}
        for label, source, _ in SUMMARY_METRICS:
            if source in nuclei.columns:
                value = row[source]
                column[label] = "–" if not np.isfinite(value) else (
                    f"{int(value)}" if source == "n_condensates" else _fmt(value, _digits(label, source)))
        for name_ in class_share_columns():
            if name_ in nuclei.columns:
                column[name_] = "–" if not np.isfinite(row[name_]) else f"{row[name_]:.1f}%"
        columns[name] = column
    table = pd.DataFrame(columns)
    all_column = summary_tables(nuclei)[0][ALL_NUCLEI]
    table[ALL_NUCLEI] = all_column.reindex(table.index if len(table.index) else all_column.index)
    table.index.name = "Metric"
    return table


def table_figure(table, title, note=""):
    """A text table as a figure (PNG / PDF page)."""
    n_rows, n_cols = table.shape
    value_width = max(2.2, 0.085 * max([len(str(c)) for c in table.columns] + [len(str(v)) for v in table.values.ravel()] + [0]))
    width = min(5.0 + value_width * n_cols, 26)
    fig, ax = plt.subplots(figsize=(width, 0.9 + 0.3 * (n_rows + 1)))
    ax.axis("off")
    cell_text = [[str(index)] + [str(v) for v in row] for index, row in zip(table.index, table.values)]
    header = ["Metric"] + [str(c) for c in table.columns]
    mpl_table = ax.table(cellText=cell_text, colLabels=header, loc="upper center", cellLoc="center",
                         colWidths=[4.6 / width] + [value_width / width] * n_cols)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(8.5)
    mpl_table.scale(1.0, 1.25)
    for (row, col), cell in mpl_table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#2c3e50")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        if col == 0:
            cell.get_text().set_ha("left")
    ax.set_title(title, fontsize=12, fontweight="bold")
    if note:
        fig.text(0.5, 0.01, note, ha="center", va="bottom", fontsize=8, style="italic", color="#555555", wrap=True)
    return fig


TABLE_NOTE = ("Unit = nucleus (ROI). Means: mean ± SD over nuclei; medians: median (IQR) over nuclei. "
              "Class shares: mean over nuclei of each nucleus's % of QC-passed objects.")


def class_totals_figure(nuclei, column_col, order, colors):
    """Total Elongated vs Globular and Total Core / No gradient / Shell: bar = mean over nuclei,
    error bar = SD, dots = single nuclei."""
    fig, axes = plt.subplots(1, 2, figsize=(max(12, 3.5 * len(order) + 4), 6))
    rng = np.random.default_rng(0)
    for ax, names, title in (
        (axes[0], SHAPE_CLASSES, "Total Elongated vs Globular"),
        (axes[1], PROFILE_CLASSES, "Total Core vs No Gradient vs Shell"),
    ):
        columns = [f"Total {name}" for name in names]
        width = 0.8 / len(columns)
        for i, condition in enumerate(order):
            group = nuclei[nuclei[column_col] == condition]
            for j, column in enumerate(columns):
                x = i + (j - (len(columns) - 1) / 2) * width
                values = pd.to_numeric(group[column], errors="coerce").dropna()
                mean = values.mean() if len(values) else 0.0
                sd = values.std(ddof=1) if len(values) > 1 else 0.0
                #A percentage stays within 0–100 %, so the error bar is cut there.
                yerr = [[min(sd, mean)], [min(sd, 100.0 - mean)]]
                ax.bar(x, mean, width=width * 0.9, color=colors[condition], alpha=0.35 + 0.65 * (j == 0),
                       edgecolor=colors[condition], yerr=yerr, capsize=4)
                ax.scatter(x + rng.uniform(-width * 0.2, width * 0.2, len(values)), values, s=18, color="white",
                           edgecolors="#333333", linewidths=0.8, zorder=3)
        ticks, labels = [], []
        for i, condition in enumerate(order):
            for j, name in enumerate(names):
                ticks.append(i + (j - (len(columns) - 1) / 2) * width)
                labels.append(name.split(" (")[0])
        ax.set_xticks(ticks, labels=labels, rotation=45, ha="right", fontsize=8)
        for i, condition in enumerate(order):
            n = int((nuclei[column_col] == condition).sum())
            ax.text(i, 104, f"{condition}\nn = {n} ROI", ha="center", va="bottom", fontsize=9, fontweight="bold")
            if i:
                ax.axvline(i - 0.5, color="#999999", linestyle=":", linewidth=1)
        ax.set_ylim(0, 100)
        ax.set_ylabel("QC-passed objects in class (%)")
        ax.set_title(title, fontsize=12, fontweight="bold", pad=40)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Bars: mean over nuclei ± SD; dots: single nuclei", fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig
