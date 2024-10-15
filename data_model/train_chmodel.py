import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier

from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.metrics import make_scorer, f1_score
from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.pipeline import Pipeline
import itertools
from sklearn.metrics import confusion_matrix, recall_score, classification_report, roc_curve, auc
from sklearn.preprocessing import LabelBinarizer
from math import sqrt
from joblib import dump, load
from keras.models import load_model, save_model

import keras
from keras.models import Sequential
from keras.layers import Dense, BatchNormalization, Dropout 
from tensorflow.keras.optimizers import SGD, Adam, Nadam, RMSprop
from keras import initializers
from keras import regularizers
from keras.layers import Dense, Activation, Dropout, LSTM
from keras.wrappers.scikit_learn import KerasClassifier


train_file = 'chmodel_train_set'
test_file = 'chmodel_test_set'

# Read training set and test set data
train_data = pd.read_csv(train_file)
test_data = pd.read_csv(test_file)

# Identify the feature and label columns
features_to_drop = ['metalid', 'pdbid', 'resi_type', 'exp_tag']
features_columns = train_data.columns.drop(features_to_drop)

X_train = train_data[features_columns]
y_train = train_data['exp_tag']
X_test = test_data[features_columns]
y_test = test_data['exp_tag']


# Logistic Regression
LR = LogisticRegression(solver='liblinear', C=0.1, penalty='l2', random_state=1)

# Decision Tree
DT = DecisionTreeClassifier(min_samples_split=5, min_samples_leaf=4, max_depth=30, random_state=1)

# Neural Network (MLPClassifier)
MLP = MLPClassifier(solver='adam', hidden_layer_sizes=(100,), alpha=0.001, activation='relu', random_state=1)

# Support Vector Machine
SVM = SVC(kernel='linear', gamma='auto', C=0.1, probability=True, random_state=1)

LR.fit(X_train, y_train)
DT.fit(X_train, y_train)
MLP.fit(X_train, y_train)
SVM.fit(X_train, y_train)

models = {'LR':LR,
        'DT':DT,
        'MLP':MLP,
        'SVM':SVM}

dump(models, 'ch_models.joblib')

# Keras Model
def create_keras_model(optimizer='rmsprop', learning_rate=0.001):
    model = Sequential()
    model.add(Dense(128, input_shape=(X_train.shape[1],), activation="relu"))
    model.add(BatchNormalization())
    model.add(Dropout(0.4))
    model.add(Dense(64, activation="relu"))
    model.add(BatchNormalization())
    model.add(Dropout(0.4))
    model.add(Dense(1, activation="sigmoid", kernel_regularizer=regularizers.l2(0.05)))
    opt = RMSprop(learning_rate=learning_rate)
    model.compile(loss="binary_crossentropy", optimizer=opt, metrics=["accuracy"])
    return model

FFNN = KerasClassifier(build_fn=create_keras_model, optimizer='rmsprop', learning_rate=0.001, epochs=200, batch_size=32, verbose=0)
FFNN.fit(X_train, y_train)
FFNN.model.save('ch_FFNN.h5')
