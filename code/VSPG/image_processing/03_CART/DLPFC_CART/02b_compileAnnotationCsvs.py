#   Standardize and clean manual cell-type labels and re-write CSVs for
#   immediate downstream use
#
#   1. Drop NA labels and use consistent names for cell types

import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

import pandas as pd
import numpy as np
import pyhere
from pathlib import Path

################################################################################
#   Variable definitions
################################################################################

df_path = pyhere.here("processed-data", "VSPG", "image_processing", "03_CART", "DLPFC_CART", "{}" + "_expanded_df.csv")
manual_label_path = pyhere.here('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/02-cellpose/',  "{}", "manual_labels_clean.csv")
manual_label_path_out = pyhere.here('processed-data', 'VSPG', 'image_processing', '02_samui', 'samui_manualAnnotation', "annotations_DLPFC.csv")

#   Define expected cell-type labels (after cleaning), expected columns
#   in the fluorescence-intensity data frame, and expected number of counts
#   for each cell-type label
expected_df_cols = ['id','gfap', 'neun', 'olig2', 'tmem119','area', 'y', 'x']
#expected_label_counts = 30


################################################################################
#   Clean data
################################################################################

#   Read in the list of sample IDs for the ID data
sample_ids = pd.Series(['V10B01-087_A1', 'V10B01-087_B1', 'V10B01-087_C1', 'V10B01-087_D1'], dtype=object)
brnums = ['Br2720_Ant_IF', 'Br6432_Ant_IF', 'Br6522_Ant_IF', 'Br8667_Post_IF']
df = pd.DataFrame()

for sample_id, brnum in zip(sample_ids, brnums):
    #   Determine paths for this sample
    this_df_path = str(df_path).format(sample_id)
    this_manual_label_path = str(manual_label_path).format(brnum)
    #   Read in all required CSVs
    this_df = pd.read_csv(this_df_path, names = expected_df_cols, header=0)
    this_manual_labels = pd.read_csv(this_manual_label_path) 
    this_manual_labels = this_manual_labels[~this_manual_labels.duplicated(subset='id', keep=False)]
    #   Check the columns are all as expected
    #assert all([x == y for x, y in zip(this_df.columns.tolist(), expected_df_cols)])
    #   Fix indices (index rows by cell ID)
    this_manual_labels.index = this_manual_labels['id']
    this_df.index = this_df['id']
    this_df['label'] = this_manual_labels['label']
    this_df['label_sample'] = this_df['label'] + '_' + sample_id
    df = pd.concat([df, this_df.dropna()])
#   Verify we have the correct amount of cells
#    assert(df.shape[0] == len(sample_ids) * expected_num_labels * num_cell_types)
label_counts = df['label'].value_counts()
#label_counts        
#other     120
#neuron    120
#astro     120
#micro     120
#oligo     120
#Name: count, dtype: int64
    
    #   Write a clean copy of both sets of manual labels
df.to_csv(manual_label_path_out, index = False)

