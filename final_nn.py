# -*- coding: utf-8 -*-
"""Final NN.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1m6rBMgURLNzYvtRhf1oX6gKa5e0XJ-jl
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

# for modeling
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.callbacks import EarlyStopping
from keras import backend as K
import keras

ratio = 1
df = pd.read_csv('binary_table.csv')
df.head(1)

def get_preds(deciles):
  df = pd.read_csv('binary_table.csv')
  imp = pd.read_csv('rf_imp.csv')
  imp.columns = ['factor','score']
  decile_size = int(len(imp['factor'])/10)
  imp = imp.head(deciles*decile_size)
  df = df[imp['factor'].append(pd.Series("GWAS"))]

  df_yes = df[df['GWAS'] == 1]
  df_no = df[df['GWAS'] == 0]
  to_take = int(np.round(ratio*df_yes.shape[0],0))

  taken = df_no.sample(n=to_take)
  df=df_yes.append(taken, ignore_index= True)

  X = df.drop(['GWAS','snp_maf','dist_nearest_gene','friends_ld07'], axis = 1)
  Y = df['GWAS']
  X = np.array(X)
  X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.8, random_state=20)

  model = Sequential()
  model.add(Dense(512, input_shape=(X.shape[1],), activation='relu'))
  model.add(Dense(256, activation='relu'))
  model.add(Dense(128, activation='relu'))
  model.add(Dense(64, activation='relu'))
  model.add(Dense(32, activation='relu'))
  model.add(Dense(16, activation='relu'))
  model.add(Dense(1, activation='sigmoid'))

  model.compile(optimizer='Adam', 
                loss='binary_crossentropy',
                metrics=['accuracy'])
                #metrics=[tf.keras.metrics.AUC()])

  K.set_value(model.optimizer.learning_rate, 0.0001)
  es = EarlyStopping(monitor='val_accuracy', 
                                    mode='max', # don't minimize the accuracy!
                                    patience=10,
                                    restore_best_weights=True)
  history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=0)
  y_pred = model.predict(X_test)
  print(roc_auc_score(y_test, y_pred))

  model = Sequential()
  model.add(Dense(256, input_shape=(X.shape[1],), activation='relu'))
  model.add(Dense(128, activation='relu'))
  model.add(Dense(64, activation='relu'))
  model.add(Dense(32, activation='relu'))
  model.add(Dense(16, activation='relu'))
  model.add(Dense(1, activation='sigmoid'))

  model.compile(optimizer='Adam', 
                loss='binary_crossentropy',
                metrics=['accuracy'])
                #metrics=[tf.keras.metrics.AUC()])

  K.set_value(model.optimizer.learning_rate, 0.0001)
  es = EarlyStopping(monitor='val_accuracy', 
                                    mode='max', # don't minimize the accuracy!
                                    patience=10,
                                    restore_best_weights=True)
  history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=0)
  y_pred = model.predict(X_test)
  print(roc_auc_score(y_test, y_pred))

  model = Sequential()
  model.add(Dense(128, input_shape=(X.shape[1],), activation='relu'))
  model.add(Dense(64, activation='relu'))
  model.add(Dense(32, activation='relu'))
  model.add(Dense(16, activation='relu'))
  model.add(Dense(1, activation='sigmoid'))

  model.compile(optimizer='Adam', 
                loss='binary_crossentropy',
                metrics=['accuracy'])
                #metrics=[tf.keras.metrics.AUC()])

  K.set_value(model.optimizer.learning_rate, 0.0001)
  es = EarlyStopping(monitor='val_accuracy', 
                                    mode='max', # don't minimize the accuracy!
                                    patience=10,
                                    restore_best_weights=True)
  history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=0)
  y_pred = model.predict(X_test)
  print(roc_auc_score(y_test, y_pred))

  model = Sequential()
  model.add(Dense(64, input_shape=(X.shape[1],), activation='relu'))
  model.add(Dense(32, activation='relu'))
  model.add(Dense(16, activation='relu'))
  model.add(Dense(8, activation='relu'))
  model.add(Dense(1, activation='sigmoid'))

  model.compile(optimizer='Adam', 
                loss='binary_crossentropy',
                metrics=['accuracy'])
                #metrics=[tf.keras.metrics.AUC()])

  K.set_value(model.optimizer.learning_rate, 0.0001)
  es = EarlyStopping(monitor='val_accuracy', 
                                    mode='max', # don't minimize the accuracy!
                                    patience=10,
                                    restore_best_weights=True)
  history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=0)
  y_pred = model.predict(X_test)
  print(roc_auc_score(y_test, y_pred))

  model = Sequential()
  model.add(Dense(32, input_shape=(X.shape[1],), activation='relu'))
  model.add(Dense(16, activation='relu'))
  model.add(Dense(8, activation='relu'))
  model.add(Dense(4, activation='relu'))
  model.add(Dense(1, activation='sigmoid'))

  model.compile(optimizer='Adam', 
                loss='binary_crossentropy',
                metrics=['accuracy'])
                #metrics=[tf.keras.metrics.AUC()])

  K.set_value(model.optimizer.learning_rate, 0.0001)
  es = EarlyStopping(monitor='val_accuracy', 
                                    mode='max', # don't minimize the accuracy!
                                    patience=10,
                                    restore_best_weights=True)
  history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=0)
  y_pred = model.predict(X_test)
  print(roc_auc_score(y_test, y_pred))

get_preds(10)

X = df.drop(['GWAS', 'dist_nearest_gene'], axis = 1)
Y = df['GWAS']
X = np.array(X)
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.8, random_state=20)

model = Sequential()
model.add(Dense(512, input_shape=(X.shape[1],), activation='relu'))
model.add(Dense(256, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(32, activation='relu'))
model.add(Dense(16, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.summary()

model.compile(optimizer='Adam', 
              loss='binary_crossentropy',
              metrics=['accuracy'])
              #metrics=[tf.keras.metrics.AUC()])

K.set_value(model.optimizer.learning_rate, 0.0001)
es = EarlyStopping(monitor='val_accuracy', 
                                   mode='max', # don't minimize the accuracy!
                                   patience=10,
                                   restore_best_weights=True)

history = model.fit(X_train,
                    y_train,
                    callbacks=[es],
                    epochs=100, # you can set this to a big number!
                    batch_size=10,
                    validation_split=0.2,
                    shuffle=True,
                    verbose=1)

y_pred = model.predict(X_test)

print(roc_auc_score(y_test, y_pred))