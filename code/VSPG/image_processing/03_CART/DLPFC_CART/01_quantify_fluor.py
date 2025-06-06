#   Given masks from cellpose (from segmenting the DAPI channel of Visium IF
#   images), quantify mean fluorescence in each non-DAPI channel within each
#   nucleus (dilated to include a region around each nucleus), and save a pandas
#   DataFrame with these values for each nucleus.

import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy.spatial import KDTree
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import json

################################################################################
#   Variable definitions
################################################################################

img_path = pyhere.here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/VistoSeg/', '{}.tif')
mask_path = pyhere.here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/02-cellpose/masks', '{}_mask.npy')
spot_path = pyhere.here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/01_spaceranger_IF', '{}', 'outs', 'spatial', 'tissue_positions_list.csv')
scale_path = pyhere.here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/01_spaceranger_IF', '{}','outs', 'spatial','scalefactors_json.json')

out_df_path = pyhere.here('processed-data', 'VSPG', 'image_processing' ,'03_CART', 'DLPFC_CART')
Path(out_df_path).mkdir(parents=True, exist_ok=True)
plot_dir = pyhere.here('plots','VSPG', 'image_processing','03_CART', 'DLPFC_CART')
Path(plot_dir).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

names = {0: "junk", 1: "dapi", 2: "gfap", 3: "neun", 4: "olig2", 5: "tmem119"}
cell_types = {"neun": "neuron","olig2": "oligo","tmem119": "micro"}

#   Nucleus area, in number of pixels, below which a cell is ignored. Tested by
#   eye
area_threshold = 60
dilation_radius = 15
dilation_chunk_size = 6

plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots

################################################################################
#   Functions
################################################################################

#   Perform a dilation-like transformation on 'img' by total size
#   'dilation_radius'. The whole transformation actually involves a series
#    of dilations by size [chunk_size / 2] followed by inversion. The idea
#   is to expand each mask equally without bias towards masks with larger-valued
#   labels (in case of overlap, which frequently happens)
def balanced_dilation(img, dilation_radius, chunk_size, verbose = False):
    assert chunk_size % 2 == 0, 'chunk_size must be even'
    assert 2 * dilation_radius % chunk_size == 0, 'cannot break this radius into chunks'
    num_chunks =  int(2 * dilation_radius / chunk_size)
    dilation_chunked = int(dilation_radius / num_chunks)
    assert num_chunks % 2 == 1, 'must use an odd number of chunks'
    #   We'll use -1 * MAX_VALUE as a placeholder, assumed to be smaller than
    #   all elements in 'img'. Check this assumption
    MAX_VALUE = 2 ** 31 - 1
    assert(np.all(img < MAX_VALUE))
    expanded_masks = img.copy().astype(np.int32)
    for i in range(num_chunks):
        if verbose:
            print(f'Dilating by {dilation_chunked} pixels...')
        #   Make sure zero-valued elements are always treated as the smallest
        #   value possible for dilation
        zero_indices = expanded_masks == 0
        expanded_masks[zero_indices] = -1 * MAX_VALUE
        expanded_masks = ndimage.grey_dilation(
            expanded_masks,
            size = (dilation_chunked, dilation_chunked)
        )
        #   Return "zero-valued" elements to a true value of 0
        zero_indices = expanded_masks == -1 * MAX_VALUE
        expanded_masks[zero_indices] = 0
        if i < num_chunks - 1:
            if verbose:
                print('Inverting...')
            expanded_masks *= -1
    return expanded_masks.astype(img.dtype)


################################################################################
#   Analysis
################################################################################

# os.environ.get['SLURM_ARRAY_TASK_ID'] = '1'

rng = default_rng()


#-------------------------------------------------------------------------------
#   Read in sample info and adjust paths for this particular sample ID
#-------------------------------------------------------------------------------

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly
sample_ids = ['V10B01-087_A1', 'V10B01-087_B1', 'V10B01-087_C1', 'V10B01-087_D1']
slurm_array_task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1)) 
sample_id_img = sample_ids[slurm_array_task_id - 1]
img_path = str(img_path).format(sample_id_img)
mask_path = str(mask_path).format(sample_id_img)

#   Path to JSON from spaceranger including spot size for this sample
brnums = ['Br2720_Ant_IF', 'Br6432_Ant_IF', 'Br6522_Ant_IF', 'Br8667_Post_IF']
brnum = brnums[slurm_array_task_id -1]
spot_path = str(spot_path).format(brnum)
json_path = str(scale_path).format(brnum)
with open(json_path) as f: 
    json_data = json.load(f)

#-------------------------------------------------------------------------------
#   Read in spot data
#-------------------------------------------------------------------------------

#   Read in all spot data
raw = pd.read_csv(spot_path,header=None,names=["barcode", "included", "row", "col", "x", "y"],)

#   Take only spots that overlap tissue
raw = raw.iloc[raw.included[raw.included == 1].index].reset_index().drop(columns=["included", "index"])

#-------------------------------------------------------------------------------
#   Quantify mean fluorescence for each channel at each nucleus
#-------------------------------------------------------------------------------

#   Load multi-channel image and masks from segmenting DAPI channel
imgs = tifffile.imread(img_path)
masks = np.load(mask_path,allow_pickle=True)
#masks = dat['masks']

its = {
    names[i]: regionprops_table(
        masks, intensity_image=imgs[i], properties=["intensity_mean"]
       # masks, intensity_image=imgs[i], properties=["intensity_mean"]
    )["intensity_mean"]
    for i in range(2, 6)
}

#   Create a table containing the centroids and areas of each mask
#   (nucleus), and add this info to the intensities table
#   general = regionprops_table(masks, properties=["centroid", "area"])
general = regionprops_table(masks, properties=["centroid", "area"])
its["area"] = general["area"]
its["x"] = general["centroid-0"]
its["y"] = general["centroid-1"]

df = pd.DataFrame(its)

df.rename(
    {
        'x': 'y',
        'y': 'x',
        'Unnamed: 0': 'id'
    },
    axis = 1, inplace = True
)
df.to_csv(pyhere.here(out_df_path,str(sample_id_img + '_df.csv')))

#   Dilate the original masks
print(f'Dilating original masks by radius {dilation_radius} and chunk size {dilation_chunk_size}.')
expanded_masks = balanced_dilation(masks, dilation_radius, dilation_chunk_size)


#   Quantify the mean image fluorescence intensity at each nucleus
#   identified by segmenting the DAPI channel. This is done for each
#   (non-lipofuscin, non-DAPI) channel
its = {
    names[i]: regionprops_table(
        expanded_masks, intensity_image=imgs[i], properties=["intensity_mean"]
       # masks, intensity_image=imgs[i], properties=["intensity_mean"]
    )["intensity_mean"]
    for i in range(2, 6)
}

#   Create a table containing the centroids and areas of each mask
#   (nucleus), and add this info to the intensities table
#   general = regionprops_table(masks, properties=["centroid", "area"])
general = regionprops_table(expanded_masks, properties=["centroid", "area"])
its["area"] = general["area"]
its["x"] = general["centroid-0"]
its["y"] = general["centroid-1"]

df = pd.DataFrame(its)
df.rename(
    {
        'x': 'y',
        'y': 'x',
        'Unnamed: 0': 'id'
    },
    axis = 1, inplace = True
)
df.to_csv(pyhere.here(out_df_path,str(sample_id_img + '_expanded_df.csv')))

#-------------------------------------------------------------------------------
#   Exploratory plot: show the distribution of masks over spots
#-------------------------------------------------------------------------------

#   Plot mask spatial distribution vs. spot distribution; there should be
#   quite a bit of overlap
plt.clf()
plt.scatter(raw["x"], raw["y"], 2)
plt.scatter(df["x"], df["y"], 2)
plt.savefig(
    os.path.join(
        plot_dir, f'mask_spot_overlap_{sample_id_img}.{plot_file_type}'
    )
)

