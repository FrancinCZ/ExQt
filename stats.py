import numpy as np
import pandas as pd


def generate_excel_stats():
    df = pd.read_csv("Final_Output_Batch.csv")


    if "volume_bio_nm3" in df.columns:
        is_3d = True
  
        df["equivalent_diameter_nm"] = ((6 * df["volume_bio_nm3"]) / np.pi) ** (
            1 / 3
        )
        size_col = "volume_bio_nm3"
        size_unit = "nm³"
    elif "area_bio_nm2" in df.columns:
        is_3d = False

        df["equivalent_diameter_nm"] = 2 * np.sqrt(
            df["area_bio_nm2"] / np.pi
        )
        size_col = "area_bio_nm2"
        size_unit = "nm²"
    else:
        raise KeyError("V CSV chybí sloupec 'volume_bio_nm3' i 'area_bio_nm2'!")

  
    df["equivalent_diameter_um"] = df["equivalent_diameter_nm"] / 1000.0


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


    with pd.ExcelWriter(
        "ExM_Analysis_Detailed_Stats.xlsx", engine="openpyxl"
    ) as writer:
        df.to_csv("Final_Output_Batch_With_Diameters.csv", index=False)
        summary_stats.to_excel(
            writer, sheet_name="Per_Cell_Summary", index=False
        )
        df.to_excel(writer, sheet_name="All_Condensates", index=False)

    print(f"Statistika úspěšně vygenerována! Režim: {'3D' if is_3d else '2D'}")


if __name__ == "__main__":
    generate_excel_stats()