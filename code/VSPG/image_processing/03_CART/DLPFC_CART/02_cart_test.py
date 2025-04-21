#   Given mean fluorescence intensities for each nucleus, classify cell types
#   present in each spot.

import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import pickle
from scipy.spatial import KDTree
import json

################################################################################
#   Paths
################################################################################


df_path = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'DLPFC_CART', 'df_cells', '{}' + '_df.csv')
spot_path = pyhere.here('processed-data', '01_spaceranger', 'spaceranger_if_2023-06-29_KMay061223', '{}', 'outs', 'spatial', 'tissue_positions.csv')
scale_path = pyhere.here('processed-data', '01_spaceranger', 'spaceranger_if_2023-06-29_KMay061223', '{}', 'outs', 'spatial','scalefactors_json.json')
clusters_path = pyhere.here('processed-data', 'VSPG', 'image_processing','03_CART', 'DLPFC_CART', 'df_cells', '{}' + '_clusters.csv')
cells_path = pyhere.here('processed-data', 'VSPG', 'image_processing','03_CART', 'DLPFC_CART', 'df_cells', '{}' + '_cell_metrics.csv')

#   Trained DecisionTreeClassifier path
#model_path = pyhere.here('processed-data', 'VSPG', 'image_processing', 'decision_tree.pkl')
model_path = Path('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/02-cellpose/decision_tree.pkl')

# os.environ['SLURM_ARRAY_TASK_ID'] = '1'
sample_ids = ['V12N28-333_A1', 'V12N28-333_B1', 'V12N28-333_C1', 'V12N28-333_D1']
slurm_array_task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1)) 
sample_id = sample_ids[slurm_array_task_id - 1]
spot_path = str(spot_path).format(sample_id)
json_path = str(scale_path).format(sample_id)
df_path = str(df_path).format(sample_id)
clusters_path = str(clusters_path).format(sample_id)
Path(clusters_path).parents[0].mkdir(parents=True, exist_ok=True)
cells_path = str(cells_path).format(sample_id)
Path(cells_path).parents[0].mkdir(parents=True, exist_ok=True)

################################################################################
#   Analysis
################################################################################


#-------------------------------------------------------------------------------
#   Call cell types (classify nuclei)
#-------------------------------------------------------------------------------
with open(model_path, 'rb') as f:
    model = pickle.load(f)

df = pd.read_csv(df_path)
#df.rename({'Unnamed: 0': 'id'}, axis = 1, inplace = True)
#df.index = df['id']
#x = df.drop(['id', 'x', 'y'], axis = 1)
df.rename({'Unnamed: 0': 'id','neun': 'NeuN','olig2': 'OLIG2',  'tmem119': 'TMEM119', 'gfap': 'GFAP'}, axis = 1, inplace = True)
feature_names = model.feature_names_in_
x = df[feature_names]

df['cell_type'] = model.predict(x)

#   Ensure cell-type names match with what we expect
cell_types = {"NeuN": "neuron","OLIG2": "oligo","TMEM119": "microglia","GFAP": "astrocyte", "other":"other"}
assert set(df['cell_type']) == set(list(cell_types.values()))

#   Save
df.to_csv(cells_path, float_format="%.3f")

#-------------------------------------------------------------------------------
#   Count cells per spot and print some related statistics
#-------------------------------------------------------------------------------
raw = pd.read_csv(spot_path,header=0,names=["barcode", "included", "row", "col", "y", "x"],)

#   Take only spots that overlap tissue
raw = raw.iloc[raw.included[raw.included == 1].index].reset_index().drop(columns=["included", "index"])

# Build KD tree for nearest neighbor search.
kd = KDTree(raw[["x", "y"]])
#   For each mask, assign a distance to the nearest spot ('dist') and index of that spot in 'df' ('idx'). Add this info to 'df'
dist, idx = kd.query(df[["x", "y"]])
dist = pd.DataFrame({"dist": dist, "idx": idx})
df = pd.concat([df, dist], axis=1)

with open(json_path) as f: 
    json_data = json.load(f)
    
frac_kept = round(100 * np.sum(df.dist < json_data['spot_diameter_fullres'] / 2) / len(df.dist), 1)
print(f'Keeping {frac_kept}% of masks, which were within spots covered by tissue.')
df = df[df.dist < json_data['spot_diameter_fullres'] / 2]

#   Count cell types in each spot
for cell_type in df['cell_type'].unique():
    raw[cell_type] = df[df['cell_type'] == cell_type].reset_index(drop=True).groupby('idx')['idx'].count().astype(int)

#   Some spots will have no cells (and NaN values for this reason). Also count the total number of cells per spot
raw.fillna(0, inplace=True)
raw['n_cells'] = raw[[x for x in df['cell_type'].unique()]].sum(axis=1)

#   Correct the data type
for column_name in list(df['cell_type'].unique()) + ['n_cells']:
    raw[column_name] = raw[column_name].astype(int)

#   Print number of cells of each type
for cell_type in df['cell_type'].unique():
    print(f'Total number of {cell_type} cells: {raw[cell_type].sum()}')


#   Print some quick stats about the cell counts per spot
prop_nonzero = round(100 * np.count_nonzero(raw['n_cells']) / raw.shape[0], 1)
print(f"Percentage of spots with at least one cell: {prop_nonzero}%")
print(f"Mean number of cells per spot: {round(np.mean(raw['n_cells']), 3)}")
print(f"Standard deviation (num cells per spot): {round(np.std(raw['n_cells']), 2)}")
print(f"Max number of cells per spot: {np.max(raw['n_cells'])}")

#-------------------------------------------------------------------------------
#   Export spot-level table as a 'clusters.csv' file
#-------------------------------------------------------------------------------

#   Make compatible with spatialLIBD 'clusters.csv' format
raw.index = raw['barcode'] + '_' + sample_id
raw.index.name = 'key'
raw.drop(['row', 'col', 'x', 'y', 'barcode'], axis = 1, inplace = True)

#   Save
raw.to_csv(clusters_path, float_format="%.3f")