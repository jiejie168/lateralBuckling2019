__author__ = 'Jie'
# coding=utf8
'''
This linear machine learning method is originally from the coursera course "Applied ML in Python", week2.
assignment. It can be further used in my own research. 20/01/2020
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression, Lasso
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split
from sklearn.metrics.regression import r2_score


np.random.seed(0)
n = 15
x = np.linspace(0,10,n) + np.random.randn(n)/5
y = np.sin(x)+x/6 + np.random.randn(n)/10
X_train, X_test, y_train, y_test = train_test_split(x, y, random_state=0)

# y_predict=pd.DataFrame()
def answer_one():
    degrees=np.array([1,3,6,9])
    X_predict=np.linspace(0,10,100)
    X_predict=X_predict.reshape(-1,1)
    predicts=pd.DataFrame()

    for i in degrees:
        poly=PolynomialFeatures(degree=i)
        X_train_poly=poly.fit_transform(X_train.reshape(-1,1))
        X_predict_poly=poly.fit_transform(X_predict)

        linreg=LinearRegression().fit(X_train_poly,y_train)
        w=linreg.coef_
        w=w.transpose()
        b=linreg.intercept_

        y_predict=np.dot(X_predict_poly,w)+b
        predicts['degree'+str(i)]=y_predict
    predicts_final=predicts.values.transpose()
    return predicts_final

# predicts=answer_one(X_train, X_test, y_train, y_test)
# print (predicts[1])

def plot_one(degree_predictions):
    plt.figure(figsize=(10,5))
    plt.plot(X_train, y_train, 'o', label='training data', markersize=10)
    plt.plot(X_test, y_test, 'o', label='test data', markersize=10)
    for i,degree in enumerate([1,3,6,9]):
        plt.plot(np.linspace(0,10,100), degree_predictions[i], alpha=0.8, lw=2, label='degree={}'.format(degree))
    plt.ylim(-1,2.5)
    plt.legend(loc=4)

# plot_one(answer_one(X_train, X_test, y_train, y_test))
# plt.show()

def answer_two():
    R2_train=[]
    R2_test=[]

    for i in np.arange(0,10):
        poly=PolynomialFeatures(degree=i)
        X_train_poly=poly.fit_transform(X_train.reshape(-1,1))
        X_test_poly=poly.fit_transform(X_test.reshape(-1,1))

        linreg=LinearRegression().fit(X_train_poly,y_train)
        w=linreg.coef_
        w=w.transpose()
        b=linreg.intercept_

        y_predict_train=np.dot(X_train_poly,w)+b
        y_predict_test=np.dot(X_test_poly,w)+b

        r2_train=r2_score(y_train,y_predict_train)
        r2_test=r2_score(y_test,y_predict_test)
        R2_train.append(r2_train)
        R2_test.append(r2_test)

    R2_train=np.array(R2_train)
    R2_test=np.array(R2_test)
    return (R2_train,R2_test)

# R2_train,R2_test=answer_two(X_train, X_test, y_train, y_test,9)
# print (R2_train)
# print (R2_test)
# print (R2_test.shape)

def plot_r2():
    plt.figure(figsize=(10,5))
    R2_train,R2_test=answer_two(X_train, X_test, y_train, y_test,9)
    plt.plot(np.arange(0,10),R2_train,'-o',lw=2,label="train data")
    plt.plot(np.arange(0,10),R2_test,'-o',lw=2,label="test data")

    ax=plt.gca()
    ax.set_xlabel("Levels of polynomial features")
    ax.set_ylabel("R2")
    plt.legend(loc=2)

# plot_r2()
# plt.show()

def answer_four():
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import Lasso, LinearRegression
    from sklearn.metrics.regression import r2_score

    poly=PolynomialFeatures(degree=12)
    X_train_poly=poly.fit_transform(X_train.reshape(-1,1) )
    X_test_poly=poly.fit_transform(X_test.reshape(-1,1))
    linreg=LinearRegression().fit(X_train_poly,y_train)
    lasreg=Lasso(alpha=0.01,max_iter=10000).fit(X_train_poly,y_train)

    w1=linreg.coef_
    w1=w1.transpose()
    b1=linreg.intercept_
    y1_predict_test=np.dot(X_test_poly,w1)+b1
    r2_test1=r2_score(y_test,y1_predict_test)

    w2=lasreg.coef_
    w2=w2.transpose()
    b2=lasreg.intercept_
    y2_predict_test=np.dot(X_test_poly,w2)+b2
    r2_test2=r2_score(y_test,y2_predict_test)

    return (r2_test1,r2_test2)

r2_test1,r2_test2=answer_four(X_train, X_test, y_train, y_test)
print (r2_test1)
print (r2_test2)
