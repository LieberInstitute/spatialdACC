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

#   Names of continuous features expected to be columns in the observation data (colData) of the AnnData
# spe_cont_features = ['NDAPI', 'CNDAPI', 'PDAPI']


#spe_path = here('processed-data', 'spot_deconvo', 'shared_utilities', 'spe.h5ad')
IMG_path = here('processed-data', 'VSPG', 'image_processing', 'samui_input', '{}.tif')
img_dir = here('processed-data', 'VSPG', 'image_processing', 'samui_input')
coord_path =  here('processed-data', 'VSPG', 'image_processing', 'counts', 'df_unfiltered', '{}_df.csv')
JSON_path = here('processed-data','01_spaceranger','spaceranger_if_2023-06-29_KMay061223', '{}', 'outs', 'spatial','scalefactors_json.json')
OUT_dir = here('processed-data', 'VSPG', 'image_processing', 'samui_manualAnnotation', '{}')

################################################################################
#   Read in sample info and clean
################################################################################
tif_files = glob.glob(os.path.join(img_dir, "*_mask1_cp_masks.png.tif"))
# Read the SLURM_ARRAY_TASK_ID from environment variables
slurm_array_task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))  # Default to 1 if not set
sample_path = tif_files[slurm_array_task_id - 1]
sample_path = tif_files[4 - 1]
sample_id = os.path.basename(sample_path).replace('_mask1_cp_masks.png.tif', '')

#   Update paths for this sample ID
out_dir = Path(str(OUT_dir).format(sample_id))
json_path = Path(str(JSON_path).format(sample_id))
img_path = Path(str(IMG_path).format(sample_id))
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
this_sample.add_coords(coords[['x', 'y']],name = "coords", mPerPx = m_per_px, size = 5e-6)

#   Add the IF image for this sample
this_sample.add_image( tiff = img_path, channels = img_channels, scale = m_per_px, defaultChannels = default_channels)
this_sample.add_csv_feature(coords[inten_features], name = "meanIntensities", coordName = "coords", dataType = "quantitative")
this_sample.set_default_feature(group = "meanIntensities", feature = "gfap")
this_sample.write()
