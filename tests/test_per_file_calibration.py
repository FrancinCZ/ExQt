import numpy as np
import pytest
import tifffile
from pathlib import Path
from size_preview import collect_size_preview
from Batch import process_condensates



def test_collect_size_preview_per_file_calibration(tmp_path):
    ome_1 = '<?xml version="1.0" encoding="UTF-8"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"><Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" Type="uint16" SizeX="20" SizeY="20" SizeZ="5" SizeC="1" SizeT="1" PhysicalSizeX="0.2" PhysicalSizeXUnit="um" PhysicalSizeY="0.2" PhysicalSizeYUnit="um" PhysicalSizeZ="0.3" PhysicalSizeZUnit="um"></Pixels></Image></OME>'
    ome_2 = '<?xml version="1.0" encoding="UTF-8"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"><Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" Type="uint16" SizeX="20" SizeY="20" SizeZ="5" SizeC="1" SizeT="1" PhysicalSizeX="0.1" PhysicalSizeXUnit="um" PhysicalSizeY="0.1" PhysicalSizeYUnit="um" PhysicalSizeZ="0.1" PhysicalSizeZUnit="um"></Pixels></Image></OME>'


    raw_data = np.ones((5, 20, 20), dtype=np.uint16) * 100
    mask_data = np.zeros((5, 20, 20), dtype=np.uint8)
    mask_data[1:3, 5:7, 5:7] = 1


    tifffile.imwrite(tmp_path / "sample1.tif", raw_data, description=ome_1)
    tifffile.imwrite(tmp_path / "sample1_Mask.tif", mask_data)

    tifffile.imwrite(tmp_path / "sample2.tif", raw_data, description=ome_2)
    tifffile.imwrite(tmp_path / "sample2_Mask.tif", mask_data)

    df = collect_size_preview(
        input_folder=tmp_path,
        mode="3d",
        expansion_factor=1.0,
        pixel_size_nm=500.0,
        z_step_nm=500.0,
        min_voxels=1,
        signal_channel=0,
    )

    assert len(df) == 2
    row1 = df[df["filename"] == "sample1.tif"].iloc[0]
    row2 = df[df["filename"] == "sample2.tif"].iloc[0]


    assert row1["raw_size"] == 8
    assert row2["raw_size"] == 8

    np.testing.assert_allclose(row1["biological_size"], 0.096, rtol=1e-3)
    np.testing.assert_allclose(row2["biological_size"], 0.008, rtol=1e-3)



def test_process_condensates_records_applied_calibration(tmp_path):
    raw_data = np.ones((5, 20, 20), dtype=np.uint16) * 100
    mask_data = np.zeros((5, 20, 20), dtype=np.uint8)
    mask_data[1:3, 5:7, 5:7] = 1

    raw_path = tmp_path / "test_raw.tif"
    mask_path = tmp_path / "test_raw_Mask.tif"
    tifffile.imwrite(raw_path, raw_data)
    tifffile.imwrite(mask_path, mask_data)

    df = process_condensates(
        tif_path=raw_path,
        mask_path=mask_path,
        mode="3d",
        expansion_factor=2.0,
        min_voxels=1,
        auto_roi=True,
        pixel_size_nm=123.4,
        z_step_nm=567.8,
        signal_channel=0,
        dapi_channel=0,
    )

    assert not df.empty
    assert "applied_pixel_size_nm" in df.columns
    assert "applied_z_step_nm" in df.columns
    assert df["applied_pixel_size_nm"].iloc[0] == pytest.approx(123.4)
    assert df["applied_z_step_nm"].iloc[0] == pytest.approx(567.8)
