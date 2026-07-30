import tifffile
import os
import glob


input_folder = r"C:\Users\franc\Desktop\SON + PI45P2"
output_folder = r"C:\Users\franc\Desktop\SON"


os.makedirs(output_folder, exist_ok=True)


for filepath in glob.glob(os.path.join(input_folder, "*.tif")):
    filename = os.path.basename(filepath)
    print(f"Zpracovávám: {filename}")
    

    image = tifffile.imread(filepath)
    

    channel_son = image[:, 0, :, :]
    channel_pip2 = image[:, 1, :, :]
    

    name_without_ext = os.path.splitext(filename)[0]
    out_son = os.path.join(output_folder, f"{name_without_ext}_SON.tif")
    out_pip2 = os.path.join(output_folder, f"{name_without_ext}_PIP2.tif")
    

    tifffile.imwrite(out_son, channel_son)
    tifffile.imwrite(out_pip2, channel_pip2)

print("Hotovo! Všechny kanály jsou rozdělené.")