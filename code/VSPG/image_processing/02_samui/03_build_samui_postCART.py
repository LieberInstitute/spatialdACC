# module load loopy/1.0.0-next.8

from pathlib import Path
import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
from pyhere import here
import json

import scanpy as sc
import numpy as np
import pandas as pd
import glob

from rasterio import Affine
from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

spot_diameter_m = 55e-6 # 5-micrometer diameter for Visium spot
img_channels = ['AF', 'DAPI', 'GFAP', 'NeuN', 'OLIG2', 'TMEM119']
default_channels = {'blue': 'DAPI', 'green': 'NeuN', 'yellow': 'TMEM119', 'red': 'GFAP', 'magenta': 'OLIG2', 'cyan': 'AF'}
inten_features = ['gfap', 'neun', 'olig2', 'tmem119', 'area']
celltypes = ['cell_type']
#spe_path = here('processed-data', 'spot_deconvo', 'shared_utilities', 'spe.h5ad')

img_dir = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/VistoSeg/'
coord_path =  here('processed-data', 'VSPG', 'image_processing', '03_CART', 'DLPFC_CART', '{}_cell_metrics.csv')
#JSON_path = here('processed-data','01_spaceranger','spaceranger_if_2023-06-29_KMay061223', '{}', 'outs', 'spatial','scalefactors_json.json')
JSON_path = here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/01_spaceranger_IF/', '{}', 'outs', 'spatial','scalefactors_json.json') )
OUT_dir = here('processed-data', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', '{}')

################################################################################
#   Read in sample info and clean
################################################################################
tif_files = glob.glob(os.path.join(img_dir, "*1.tif"))
# Read the SLURM_ARRAY_TASK_ID from environment variables
task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))  # Default to 1 if not set
sample_id = os.path.basename(tif_files[task_id - 1]).replace('.tif', '')

sample_ids = pd.Series(['V10B01-087_A1', 'V10B01-087_B1', 'V10B01-087_C1', 'V10B01-087_D1'], dtype=object)
brnums = pd.Series(['Br2720_Ant_IF', 'Br6432_Ant_IF', 'Br6522_Ant_IF', 'Br8667_Post_IF'], dtype=object)
idx = sample_ids[sample_ids == sample_id].index[0]
brnum = brnums.iloc[idx]

# Print both
print(f"sample_id: {sample_id}")
print(f"brnum: {brnum}")

#   Update paths for this sample ID
out_dir = Path(str(OUT_dir).format(sample_id))
json_path = Path(str(JSON_path).format(brnum))
img_path = Path(tif_files[task_id - 1])
coord_path = Path(str(coord_path).format(sample_id))
out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, json_path, img_path]])

################################################################################
#   Read in scale-factors info
################################################################################

#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   extract coords for cells
################################################################################

coords=pd.read_csv(coord_path)
coords.index = coords.index.astype(str)
################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id, path = out_dir)
this_sample.add_coords(coords[['x', 'y']],name = "coords", mPerPx = m_per_px, size = 10e-6)

#   Add the IF image for this sample
this_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px, defaultChannels = default_channels)
this_sample.add_csv_feature(coords[inten_features], name = "meanIntensities", coordName = "coords", dataType = "quantitative")
this_sample.add_csv_feature(coords[celltypes], name = "celltypes", coordName = "coords", dataType = "categorical")
this_sample.set_default_feature(group = "celltypes", feature = "cell_type")
this_sample.write()
