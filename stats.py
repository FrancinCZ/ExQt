import numpy as np
import pandas as pd
from scipy.stats import skew

def generate_excel_stats():


    df = pd.read_csv("Final_Output_Batch.csv")

    df["equivalent_diameter_nm"] = ((6 * df["volume_bio_nm3"]) / np.pi) ** (1/3)


    metrics = {
        "Biological Volume (nm³)": df["volume_bio_nm3"],
        "Equivalent Spherical Diameter (nm)": df["equivalent_diameter_nm"],
        "Volume (voxels)": df["volume_voxels"],
        "Mean Intensity (a.u.)": df["mean_intensity"],
        "Integrated Density (a.u.)": df["integrated_density"],
    }

    stats_data = []

    for name, series in metrics.items():

        series = series.dropna()
        
        count_val = len(series)
        mean_val = series.mean()
        std_val = series.std()
        sem_val = std_val / np.sqrt(count_val) if count_val > 0 else 0
        median_val = series.median()
        

        if count_val > 0:
            q75, q25 = np.percentile(series, [75, 25])
            iqr_val = q75 - q25
            min_val = series.min()
            max_val = series.max()
            skew_val = skew(series)
        else:
            q75 = q25 = iqr_val = min_val = max_val = skew_val = np.nan

        stats_data.append(
            {
                "Metric": name,
                "Count (N)": count_val,
                "Mean": mean_val,
                "Std. Deviation (SD)": std_val,
                "Std. Error (SEM)": sem_val,
                "Median": median_val,
                "Interquartile Range (IQR)": iqr_val,
                "Min": min_val,
                "Max": max_val,
                "Skewness": skew_val,
            }
        )

    summary_df = pd.DataFrame(stats_data)


    output_xlsx = "ExM_3D_Summary_Statistics.xlsx"
    summary_df.to_excel(output_xlsx, index=False)

    print(f"File saved as '{output_xlsx}'")

if __name__ == "__main__":
    generate_excel_stats()