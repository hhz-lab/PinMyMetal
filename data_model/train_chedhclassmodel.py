import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier
from sklearn.ensemble import VotingClassifier
#from scipy.stats import uniform
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold

from sklearn.pipeline import Pipeline
import itertools
from sklearn.metrics import confusion_matrix, recall_score, classification_report, roc_curve, auc
from sklearn.preprocessing import LabelBinarizer
from math import sqrt
from joblib import dump, load
from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.metrics import make_scorer, f1_score


train_file = 'classmodel_train_set'
test_file = 'classmodel_test_set'

# Read training set and test set data
train_data = pd.read_csv(train_file)
test_data = pd.read_csv(test_file)

train_data = train_data.dropna()
test_data = test_data.dropna()

# Identify the feature and label columns
features_to_drop = ['metalid', 'pdbid', 'residueid_ion', 'label_metal', 'source', 'ched_count']
features_columns = train_data.columns.drop(features_to_drop)

X_train = train_data[features_columns]
y_train = train_data['label_metal']
X_test = test_data[features_columns]
y_test = test_data['label_metal']


# Logistic Regression with best parameters
LR = LogisticRegression(solver='liblinear', C=1, penalty='l2', random_state=1)

# Balanced Random Forest Classifier with best parameters
RF = BalancedRandomForestClassifier(n_estimators=50, max_depth=10, random_state=1)

# Neural Network (MLPClassifier) with best parameters
MLP = MLPClassifier(solver='adam', hidden_layer_sizes=(100,), alpha=0.0001, activation='relu', random_state=1)

# Support Vector Machine with best parameters
SVM = SVC(kernel='rbf', gamma='scale', C=100, probability=True, random_state=1)

# Easy Ensemble Classifier with best parameters
Easy = EasyEnsembleClassifier(n_estimators=100, random_state=1)

ensemble_model = VotingClassifier(estimators=[
    ('LR', LR),
    ('RF', RF),
    ('MLP', MLP),
    ('SVM', SVM),
    ('Easy', Easy)
], voting='soft')

ensemble_model.fit(X_train, y_train)
dump(ensemble_model, 'chedh_classmodel.joblib')

