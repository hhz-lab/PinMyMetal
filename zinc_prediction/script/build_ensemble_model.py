import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler

from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
import itertools
from sklearn.metrics import confusion_matrix, recall_score, classification_report
from math import sqrt
from joblib import dump, load
from keras.models import load_model, save_model
from sklearn.model_selection import StratifiedShuffleSplit

import keras
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import BatchNormalization
from tensorflow.keras.optimizers import SGD, Adam, Nadam
from keras import initializers
from keras import regularizers
from sklearn.metrics import classification_report, roc_curve, auc
from sklearn.preprocessing import LabelBinarizer
from keras.layers import Dense, Activation, Dropout, LSTM
from keras.wrappers.scikit_learn import KerasClassifier

filename="mvc_train_data"

raw = pd.read_csv(filename)
raw = pd.DataFrame(raw)

#The data set is divided into training set and test set

from sklearn.model_selection import StratifiedShuffleSplit

features_data = raw['resi_type']
labels_data = raw["exp_tag"]

features = raw.loc[:, 'resname_a':'s21']
labels=raw['exp_tag']


split = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=42)
 
for train_index, test_index in split.split(features_data, labels_data):
    features_train, features_test = features.iloc[train_index], features.iloc[test_index]
    labels_train, labels_test = labels.iloc[train_index], labels.iloc[test_index]

X_train, X_test, y_train, y_test= features_train, features_test, labels_train, labels_test

def printing_Kfold_scores(x_train_data, y_train_data):
    fold = KFold(10,shuffle=False) 
    c_param_range = [0.01, 0.1, 1, 10, 100] 
    results_table = pd.DataFrame(index = range(len(c_param_range),2), columns = ['C_parameter','Mean recall score'],dtype=object)
    results_table['C_parameter'] = c_param_range
    
    j=0
    for c_param in c_param_range:
        
        recall_accs = []
        for iteration, indices in enumerate(fold.split(y_train_data), start=1):
            lr = LogisticRegression(C= c_param, penalty = 'l2') 
            lr.fit(x_train_data.iloc[indices[0]],y_train_data.iloc[indices[0]].values.ravel())
            
            y_pred_undersample = lr.predict(x_train_data.iloc[indices[1]].values)
            
            recall_acc = recall_score(y_train_data.iloc[indices[1]].values,y_pred_undersample)
            recall_accs.append(recall_acc)
        
        results_table.loc[j,'Mean recall score']= np.mean(recall_accs)
        j += 1
    
    results_table['Mean recall score'] = results_table['Mean recall score'].astype('float64')  
    best_c = results_table.loc[results_table['Mean recall score'].idxmax()]['C_parameter']
    return best_c

best_c = printing_Kfold_scores(features_train,labels_train)

clf1 = LogisticRegression(penalty='l2',C=best_c,random_state=1)
clf2 = DecisionTreeClassifier(criterion='entropy', random_state=1, splitter='random')
clf3 = MLPClassifier(hidden_layer_sizes=(100, 100))
clf4 = SVC(probability=True)

clf1.fit(X_train, y_train)
clf2.fit(X_train, y_train)
clf3.fit(X_train, y_train)
clf4.fit(X_train, y_train)

def create_keras_model():
    model = Sequential()
    model.add(Dense(128, input_shape=(61,), activation="relu"))
    model.add(BatchNormalization())
    model.add(Dropout(0.4))
    model.add(Dense(64, activation="relu"))
    model.add(BatchNormalization())
    model.add(Dropout(0.4))
    model.add(Dense(1, activation="sigmoid", kernel_regularizer=regularizers.l2(0.05)))
    model.summary()
    INIT_LR = 0.005
    EPOCHS = 300
    opt = Adam(lr=INIT_LR)
    model.compile(loss="binary_crossentropy", optimizer=opt, metrics=["accuracy"])
    model._estimator_type="classifier"
    return model

keras_classifier = KerasClassifier(build_fn=create_keras_model)
keras_classifier.fit(X_train, y_train)


models = {'model1':clf1,
        'model2':clf2,
        'model3':clf3,
        'model4':clf4}

dump(models, 'models.joblib')

keras_classifier.model.save('fcnn.h5')

