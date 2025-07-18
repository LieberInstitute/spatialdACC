import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

import pandas as pd
import numpy as np
import copy
import pickle

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import tree

import pyhere
from pathlib import Path

import matplotlib.pyplot as plt
import graphviz

# sample_id = "V12D07-332_D1"
# df_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', sample_id + '_df.csv')
# predictions_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', '{}' + 'CART_predictions.csv')
dataset_path = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', 'annotation_dataset_final_expandedP.pkl')
tree_path = pyhere.here('plots', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', 'decision_tree_final_expandedP.pdf')
model_out_path = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', 'decision_tree_final_expandedP.pkl')

random_seed = 0

#   Number of cells per type to predict and output to Loopy-compatible CSV for inspection
num_examples_check = 10

################################################################################
#   Functions
################################################################################

#   Remove splits that don't differentiate classes. Credit to:
#   https://github.com/scikit-learn/scikit-learn/issues/10810#issuecomment-409028102
def prune(tree):
    tree = copy.deepcopy(tree)
    dat = tree.tree_
    nodes = range(0, dat.node_count)
    ls = dat.children_left
    rs = dat.children_right
    classes = [[list(e).index(max(e)) for e in v] for v in dat.value]
    leaves = [(ls[i] == rs[i]) for i in nodes]
    LEAF = -1
    for i in reversed(nodes):
        if leaves[i]:
            continue
        if leaves[ls[i]] and leaves[rs[i]] and classes[ls[i]] == classes[rs[i]]:
            ls[i] = rs[i] = LEAF
            leaves[i] = True
    return tree

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Read in sample IDs and manual-annotation dataset
#-------------------------------------------------------------------------------

with open(dataset_path, 'rb') as f:
    x_train, x_test, y_train, y_test = pickle.load(f)

#-------------------------------------------------------------------------------
#   Train the DecisionTreeClassifier
#-------------------------------------------------------------------------------

#   Instantiate, fit, and save the CART. Good hyperparameters were determined in
#   10-explore_models.py (via cross-validation, with final test accuracy of
#   99.0%) and are used here. Note that we still reserve part of the data as
#   test data, since we want to be able to make statements about the test
#   performance of the exact tree used for inference, and trees can
#   qualitatively change with the introduction of new data (if we used all
#   manual labels for inference).
model = tree.DecisionTreeClassifier(
    criterion = 'gini', 
    max_depth = 5,
    min_samples_leaf = 10, 
    random_state = random_seed,
    min_samples_split = 2,
    ccp_alpha = 0.005263157894736842,
    splitter = 'best'
)

model.fit(x_train, y_train)

with open(model_out_path, 'wb') as f:
    pickle.dump(model, f)

#   Compute training and test accuracy on the model
acc_train = round(100 * model.score(x_train, y_train), 1)
acc_test = round(100 * model.score(x_test, y_test), 1)
print(f'CART training accuracy: {acc_train}%.') #CART training accuracy: 88.6%.
print(f'CART test accuracy: {acc_test}%.') #CCART test accuracy: 86.0%.

#CART training accuracy: 89.2%.
#CART test accuracy: 90.5%.
#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = model.predict(x_train)
labels_test = model.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
#Training report:
#               precision    recall  f1-score   support
#
#        DAPI       0.83      0.53      0.65       234
#        GFAP       0.92      0.93      0.92       298
#        NeuN       0.96      0.97      0.96       300
#       OLIG2       0.77      0.97      0.86       299
#     TMEM119       0.96      0.96      0.96       291
#
#    accuracy                           0.89      1422
#   macro avg       0.89      0.87      0.87      1422
#weighted avg       0.89      0.89      0.88      1422
#

#Training report:
#               precision    recall  f1-score   support
#
#        DAPI       0.83      0.53      0.65       196
#        GFAP       0.95      0.90      0.92       298
#        NeuN       0.92      0.98      0.95       300
#       OLIG2       0.79      0.96      0.87       299
#     TMEM119       0.96      0.96      0.96       291
#
#    accuracy                           0.89      1384
#   macro avg       0.89      0.87      0.87      1384
#weighted avg       0.89      0.89      0.89      1384

print('Test report:\n', classification_report(y_test, labels_test))

#Test report:
#               precision    recall  f1-score   support
#
#        DAPI       0.75      0.47      0.57        58
#        GFAP       0.90      0.92      0.91        75
#        NeuN       0.90      0.96      0.93        75
#       OLIG2       0.76      0.95      0.85        75
#     TMEM119       0.96      0.92      0.94        73
#
#    accuracy                           0.86       356
#   macro avg       0.85      0.84      0.84       356
#weighted avg       0.86      0.86      0.85       356

#Test report:
#               precision    recall  f1-score   support
#
#        DAPI       0.85      0.60      0.71        48
#        GFAP       0.95      0.95      0.95        75
#        NeuN       0.94      0.96      0.95        75
#       OLIG2       0.83      0.93      0.88        75
#     TMEM119       0.93      0.97      0.95        73
#
#    accuracy                           0.90       346
#   macro avg       0.90      0.88      0.89       346
#weighted avg       0.90      0.90      0.90       346
#-------------------------------------------------------------------------------
#   Write a small CSV of model predictions for import to loopy to visually
#   verify model performance
#-------------------------------------------------------------------------------

#   Read in the original unfiltered cells for a single sample
# expected_df_cols = ['id','NeuN', 'TMEM119', 'GFAP', 'OLIG2','area', 'y', 'x']
# this_df = pd.read_csv(df_path, names = expected_df_cols, header=0)
# this_df.index = this_df['id']
# this_df = this_df.loc[:, ['NeuN','TMEM119','GFAP','OLIG2','area']]
#
# #   Row IDs were not preserved during concatenation earlier, so we'll manually
# #   drop rows in 'this_df' that are present in 'x_train' or 'x_test'
# x = pd.concat([x_train, x_test])
# this_df = this_df.merge(
#     x, how='left', indicator=True,
#     on = ['NeuN','TMEM119','GFAP','OLIG2','area']
# )
# this_df = this_df[this_df['_merge'] == 'left_only'].drop('_merge', axis=1)
# this_df['label'] = model.predict(this_df)
#
# #   For each cell type, pick [num_examples_check] cells randomly. Then form a
# #   dataframe with these cells
# small_df = pd.DataFrame()
# for cell_type in this_df['label'].unique():
#     indices = np.random.choice(
#         this_df[this_df['label'] == cell_type].index,
#         size = num_examples_check, replace = False
#     )
#     small_df = pd.concat([small_df, this_df.loc[indices]])
#
# #   Make compatible with import to Loopy and write to CSV
# small_df['id'] = small_df.index
# small_df = small_df.loc[:, ['id', 'label']]
# small_df.to_csv(str(predictions_path).format(sample_id), index = False)

#-------------------------------------------------------------------------------
#   Plot the (simplified) decision tree visually (save to PDF)
#-------------------------------------------------------------------------------
x = pd.concat([x_train, x_test])
model = prune(model)
_ = tree.plot_tree(
    model, class_names = model.classes_, feature_names = x.columns,
    filled = True, rounded = True
)
plt.savefig(tree_path)