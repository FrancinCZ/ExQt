import numpy as np
import pandas as pd
from pathlib import Path


def generate_excel_stats(csv_filename="Final_Output_HybridH5Batch_3d.csv"):
    csv_path = Path(csv_filename)
    
    if not csv_path.exists():
        print(f"Error '{csv_filename}' was not found in the current directory. Please check the file name and try again.")
        return

    folder_name = csv_path.resolve().parent.name
    

    output_excel = f"{folder_name}_Detailed_Stats.xlsx"
    output_csv = f"{folder_name}_All_Condensates_With_Diameters.csv"

    df = pd.read_csv(csv_path)

    if "volume_bio_um3" in df.columns:
        is_3d = True
        df["equivalent_diameter_um"] = ((6 * df["volume_bio_um3"]) / np.pi) ** (1 / 3)
        size_col = "volume_bio_um3"
        
    elif "area_bio_um2" in df.columns:
        is_3d = False
        df["equivalent_diameter_um"] = 2 * np.sqrt(df["area_bio_um2"] / np.pi)
        size_col = "area_bio_um2"
        
    elif "shape_metric_bio" in df.columns:
        is_3d = df["is_3d"].iloc[0] if "is_3d" in df.columns else True
        size_col = "shape_metric_bio"
        if is_3d:
            df["equivalent_diameter_um"] = ((6 * df["shape_metric_bio"]) / np.pi) ** (1 / 3)
        else:
            df["equivalent_diameter_um"] = 2 * np.sqrt(df["shape_metric_bio"] / np.pi)
    else:
        raise KeyError("In CSV the is no ('volume_bio_um3') or ('area_bio_um2')!")


    df["equivalent_diameter_nm"] = df["equivalent_diameter_um"] * 1000.0

    summary_stats = (
        df.groupby("filename")
        .agg(
            condensate_count=("object_id", "count"),
            mean_size=(size_col, "mean"),
            median_size=(size_col, "median"),
            mean_diameter_nm=("equivalent_diameter_nm", "mean"),
            median_diameter_nm=("equivalent_diameter_nm", "median"),
            mean_intensity=("mean_intensity", "mean"),
            integrated_density_sum=("integrated_density", "sum"),
        )
        .reset_index()
    )


    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:

        df.to_csv(output_csv, index=False)
        
        summary_stats.to_excel(
            writer, sheet_name="Per_Cell_Summary", index=False
        )
        df.to_excel(writer, sheet_name="All_Condensates", index=False)

    print(f"\n[Statistika úspěšně vygenerována! Režim: {'3D' if is_3d else '2D'}]")
    print(f" ➜ Excel: {output_excel}")
    print(f" ➜ CSV:   {output_csv}")


if __name__ == "__main__":
    generate_excel_stats()