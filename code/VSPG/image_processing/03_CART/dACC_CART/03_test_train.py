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

df_path = pyhere.here('processed-data', 'VSPG', 'image_processing', '02_samui', 'samui_manualAnnotation', 'annotations_Atharv_processed.csv')
df_path_out = pyhere.here('processed-data', 'VSPG', 'image_processing', '03_CART', 'dACC_CART', 'annotation_dataset.pkl')

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

 #df['label'] = df['label'].replace('AF', 'DAPI')
 #df['label_sample'] = df['label_sample'].str.replace('AF_', 'DAPI_', regex=False)
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
#Using 390 training and 98 test examples (79.92% training).
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
print(f'CART training accuracy: {acc_train}%.') #CART training accuracy: 74.1%.
print(f'CART test accuracy: {acc_test}%.') #CART test accuracy: 64.3%.
print(f'Best params: {grid.best_params_}') #Best params: {'ccp_alpha': 0.0, 'criterion': 'gini', 'max_depth': 6, 'min_samples_leaf': 20, 'min_samples_split': 2, 'splitter': 'best'}

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))

#Training report:
#               precision    recall  f1-score   support
#
#          AF       0.57      0.60      0.59        65
#        DAPI       0.57      0.45      0.50        65
#        GFAP       0.84      0.78      0.81        65
#        NeuN       0.84      0.91      0.87        65
#       OLIG2       0.72      0.79      0.75        66
#     TMEM119       0.87      0.92      0.89        64
#
#    accuracy                           0.74       390
#   macro avg       0.74      0.74      0.74       390
#weighted avg       0.73      0.74      0.74       390

print('Test report:\n', classification_report(y_test, labels_test))

#Test report:
#               precision    recall  f1-score   support
#
#          AF       0.47      0.53      0.50        17
#        DAPI       0.40      0.25      0.31        16
#        GFAP       0.71      0.75      0.73        16
#        NeuN       0.82      0.88      0.85        16
#       OLIG2       0.72      0.76      0.74        17
#     TMEM119       0.76      0.81      0.79        16
#
#    accuracy                           0.66        98
#   macro avg       0.65      0.66      0.65        98
#weighted avg       0.65      0.66      0.65        98

#-------------------------------------------------------------------------------
#   Try logistic regression
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'logisticregression__C': np.logspace(-2, 4, 20)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    LogisticRegression(random_state = random_seed, max_iter=10000)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'Logistic regression training accuracy: {acc_train}%.')
print(f'Logistic regression test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

#-------------------------------------------------------------------------------
#   Try SVM (linear and non-linear kernels)
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'svc__kernel': ['linear', 'rbf'],
        'svc__gamma': np.logspace(-5, -1, 10),
        'svc__C': np.logspace(0, 4, 5)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    svm.SVC(random_state = random_seed)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'SVM training accuracy: {acc_train}%.')
print(f'SVM test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))
