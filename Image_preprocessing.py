import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import remove_small_objects

# 📁 Paths & Parameters
input_folder = '/Users/Vishwa/Downloads/Master_Thesis/iMGLs_images/Lena_GBA1mut Mics'
output_folder = '/Users/Vishwa/Downloads/Master_Thesis/Image_quantification/Output_07_08_14'
min_nucleus_size = 100
min_speck_size = 10

os.makedirs(output_folder, exist_ok=True)

# 🔁 Loop over .nd2 files
for filename in os.listdir(input_folder):
    if not filename.endswith('.nd2'):
        continue

    filepath = os.path.join(input_folder, filename)
    image_name = filename.replace('.nd2', '')
    print(f"🧪 Processing: {image_name}")

    try:
        with ND2Reader(filepath) as nd2:
            nd2.iter_axes = 't'
            nd2.bundle_axes = 'zcyx'
            stack = np.array(nd2[0])  # Shape: [z, c, y, x]

            if stack.ndim != 4 or stack.shape[1] < 4:
                print(f"⚠️ Skipping {filename}: unexpected shape {stack.shape}")
                continue

            # 📊 Max projection per channel
            dapi_proj = np.max(stack[:, 0], axis=0)
            lamp1_proj = np.max(stack[:, 1], axis=0)
            iba1_proj = np.max(stack[:, 2], axis=0)
            mki67_proj = np.max(stack[:, 3], axis=0)

            # 🧠 DAPI nuclei segmentation
            dapi_thresh = threshold_otsu(dapi_proj)
            nuclei_mask = dapi_proj > dapi_thresh
            nuclei_mask = remove_small_objects(nuclei_mask, min_size=min_nucleus_size)
            nuclei_label = label(nuclei_mask)
            nuclei_props = regionprops(nuclei_label, intensity_image=mki67_proj)

            # 🌟 LAMP1 specks
            lamp_thresh = threshold_otsu(lamp1_proj)
            lamp_mask = lamp1_proj > lamp_thresh
            lamp_mask = remove_small_objects(lamp_mask, min_size=min_speck_size)
            lamp_label = label(lamp_mask)
            lamp_props = regionprops(lamp_label, intensity_image=lamp1_proj)

            # 📈 Adaptive MKI67 threshold
            mki67_intensities = [nuc.mean_intensity for nuc in nuclei_props]
            adaptive_thresh = np.percentile(mki67_intensities, 90)

            # 📊 Quantify each nucleus
            per_nucleus_data = []
            for nucleus in nuclei_props:
                nuc_id = nucleus.label
                minr, minc, maxr, maxc = nucleus.bbox
                cy, cx = nucleus.centroid

                speck_count = 0
                total_speck_area = 0
                total_speck_intensity = 0

                for speck in lamp_props:
                    sy, sx = speck.centroid
                    if minr <= sy <= maxr and minc <= sx <= maxc:
                        speck_count += 1
                        total_speck_area += speck.area
                        total_speck_intensity += speck.mean_intensity

                iba1_roi = iba1_proj[minr:maxr, minc:maxc]
                iba1_intensity = np.mean(iba1_roi)

                mki67_status = nucleus.mean_intensity > adaptive_thresh

                per_nucleus_data.append({
                    'Image': image_name,
                    'Nucleus_ID': nuc_id,
                    'Centroid_Y': cy,
                    'Centroid_X': cx,
                    'Speck_Count': speck_count,
                    'Total_Speck_Area': total_speck_area,
                    'Mean_Speck_Intensity': total_speck_intensity / speck_count if speck_count else 0,
                    'IBA1_Intensity': iba1_intensity,
                    'MKI67_Intensity': nucleus.mean_intensity,
                    'MKI67_Positive': int(mki67_status)
                })

            # 💾 Save CSV
            df = pd.DataFrame(per_nucleus_data)
            csv_path = os.path.join(output_folder, image_name + '_nucleus_data.csv')
            df.to_csv(csv_path, index=False)
            print(f"✅ Saved CSV: {csv_path}")

            # 🖼️ Create clean overlay (DAPI + LAMP1)
            overlay_rgb = np.stack([
                dapi_proj / np.max(dapi_proj),
                lamp1_proj / np.max(lamp1_proj),
                np.zeros_like(dapi_proj)
            ], axis=-1)

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(overlay_rgb)
            ax.set_title(image_name)
            ax.axis('off')
            fig.tight_layout()
            overlay_path = os.path.join(output_folder, image_name + '_overlay.png')
            fig.savefig(overlay_path, dpi=300)
            plt.close(fig)
            print(f"🖼️ Saved overlay: {overlay_path}")

    except Exception as e:
        print(f"❌ Error with {filename}: {e}")
        continue

print(f"🎉 All images processed. Results saved in '{output_folder}'")
