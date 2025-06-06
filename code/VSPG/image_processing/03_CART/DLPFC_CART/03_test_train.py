import os
os.chdir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import tree, svm

import pyhere
from pathlib import Path
import pickle

df_path = pyhere.here('processed-data', 'VSPG', 'image_processing', '02_samui', 'samui_manualAnnotation', 'annotations_DLPFC.csv')
#df_path_out = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', 'annotation_dataset.pkl')
df_path_out = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'DLPFC_CART', 'annotation_dataset_final_expanded.pkl')

num_cell_types = 5
test_proportion = 0.2 # for training/test split
random_seed = 0

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Preprocess and gather fluorescence data + original manual cell-type labels
#-------------------------------------------------------------------------------

df = pd.read_csv(df_path,header=0)

#   Define the inputs (features we want the model to access) and outputs to the
#   model
x = df.loc[:, ['gfap', 'neun', 'olig2', 'tmem119', 'area', 'label_sample']]
y = df['label']

#   Split data into training and test sets (80%: 20%), evenly stratified across
#   classes and sample ID
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = test_proportion, random_state = random_seed,
    stratify = x['label_sample']
)


#   Remove the column we only used to properly stratify the data
x_train.drop('label_sample', axis = 1, inplace = True)
x_test.drop('label_sample', axis = 1, inplace = True)


#-------------------------------------------------------------------------------
#   Combine the two different sets of manual annotation into one dataset and
#   write to disk
#-------------------------------------------------------------------------------

perc_train = round(100 * x_train.shape[0] / (x_train.shape[0] + x_test.shape[0]), 2)
print(f'Using {x_train.shape[0]} training and {x_test.shape[0]} test examples ({perc_train}% training).')
Using 480 training and 120 test examples (80.0% training).
#   Write the dataset to disk for later use
with open(df_path_out, 'wb') as f:
    pickle.dump((x_train, x_test, y_train, y_test), f)

#-------------------------------------------------------------------------------
#   Try a decision tree classifier
#-------------------------------------------------------------------------------
tuned_parameters = [
    {
        'criterion': ['gini', 'entropy'],
        'max_depth': [2, 3, 4, 5, 6, 8, 10],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 5, 10, 20, 30],
        'ccp_alpha': np.linspace(0, 0.05, 20),
        'splitter':['best', 'random']
    }
]

grid = GridSearchCV(
    tree.DecisionTreeClassifier(
        random_state = random_seed, splitter = 'best', class_weight = None
    ),
    tuned_parameters, cv=5, scoring = 'accuracy'
)
grid.fit(x_train, y_train)

#   Compute training and test accuracy on the best model
acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'CART training accuracy: {acc_train}%.') #CART training accuracy: 97.7%.
print(f'CART test accuracy: {acc_test}%.') #CART test accuracy: 95.0%.
print(f'Best params: {grid.best_params_}') #{'ccp_alpha': 0.002631578947368421, 'criterion': 'gini', 'max_depth': 6, 'min_samples_leaf': 2, 'min_samples_split': 2, 'splitter': 'best'}

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
Training report:
               precision    recall  f1-score   support

       astro       0.97      0.94      0.95        96
       micro       1.00      0.96      0.98        96
      neuron       0.99      1.00      0.99        96
       oligo       1.00      0.99      0.99        96
       other       0.93      1.00      0.96        96

    accuracy                           0.98       480
   macro avg       0.98      0.98      0.98       480
weighted avg       0.98      0.98      0.98       480

print('Test report:\n', classification_report(y_test, labels_test))
Test report:
               precision    recall  f1-score   support

       astro       0.88      0.96      0.92        24
       micro       0.96      1.00      0.98        24
      neuron       1.00      0.96      0.98        24
       oligo       1.00      0.92      0.96        24
       other       0.92      0.92      0.92        24

    accuracy                           0.95       120
   macro avg       0.95      0.95      0.95       120
weighted avg       0.95      0.95      0.95       120