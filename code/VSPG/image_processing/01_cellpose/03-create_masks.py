print("in python")
import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
import numpy as np
from cellpose import models, io
from cellpose.io import imread
import pyhere
import pandas as pd
from pathlib import Path
import sys
import glob
print("imported all libraries")

model_path = pyhere.here('code','VSPG','image_processing','01_cellpose','CP_20230819_163525')
#model_path = pyhere.here('code','VSPG','image_processing','CP_20220415_152031')
cell_diameter = None
channel = 1 # DAPI
print("loaded Pretrained model based on 'cyto' and iteratively refined in the GUI")

#   Visium-IF images
img_dir = pyhere.here('processed-data', 'Images', 'VistoSeg', 'Capture_areas', 'if-images')
mask_dir = pyhere.here('processed-data', 'VSPG', 'image_processing', 'cellpose_masks') #if using maddy models
Path(mask_dir).mkdir(parents=True, exist_ok=True)
print("made output directories")

#   Determine paths to IF images; read in just the one for this sample
tif_files = glob.glob(os.path.join(img_dir, "*.tif"))
# Read the SLURM_ARRAY_TASK_ID from environment variables
slurm_array_task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))  # Default to 1 if not set
sample_path = tif_files[slurm_array_task_id - 1]
# Extract the sample ID without the .tif extension
sample_id = os.path.basename(sample_path).replace('.tif', '')
print(f'reading {sample_id}')

img = imread(pyhere.here(img_dir, sample_id + '.tif'))[channel, :, :]
model = models.CellposeModel(gpu = True, pretrained_model = model_path.as_posix())
masks, flows, styles = model.eval(img, diameter=cell_diameter)
print("processed images wih pretrained model")

#   Save masks
mask_npy = str(pyhere.here(mask_dir, f'{sample_id}_mask1.npy'))
np.save(mask_npy, masks)
print("saved npy array")
#   Save PNG version of the masks to visually inspect results
mask_png = str(pyhere.here(mask_dir, f'{sample_id}_mask1.png'))
io.save_to_png(img, masks, flows, mask_png)
print("saved png")
