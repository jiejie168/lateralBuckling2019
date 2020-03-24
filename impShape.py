__author__ = 'Jie'
'''
This script is used for the production of imperfection shapes.
'''
import numpy as np
from matplotlib import pyplot as plt
from pandas import Series, DataFrame
import pandas as pd
from scipy.optimize import fminbound
from sympy import Symbol, Derivative
from sympy import *
import csv
import matplotlib.transforms
import matplotlib as mpl
import matplotlib.animation as animation
from matplotlib.widgets import  Slider
from matplotlib.widgets import TextBox

L=100
w0m=10
s=0.01
num=8.9868   # nL=8.986
x=np.linspace(-0.5,0.5,200)
w0=s/15.7*(1+8.9868**2/8-(num*x)**2/2-np.cos(num*x)/np.cos(8.9868/2))

w0_cos=s*np.cos(np.pi*x)

plt.plot(x,w0,'k-',linewidth=2)
plt.plot(x,w0_cos,linewidth=2)
ax=plt.gca()
ax.set_ylim(0,)
ax.xaxis.set_ticks(np.arange(-0.5,0.55,0.1))
ax.legend(['Hobbs\'shape','Cosinusoidal shape'],frameon=False,loc='upper left')

ax.set_xlabel("Normalized pipe length")

xmajorFormatter = mpl.ticker.FormatStrFormatter('%3.1f')
ymajorFormatter=mpl.ticker.FormatStrFormatter('%4.3f')
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.yaxis.set_major_formatter(ymajorFormatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position(('data',0))
ax.tick_params(left=False,right=False,top=False,bottom=True,direction='out',pad=4,labelsize=10)

ax.set_ylabel("Imperfection amplitude")
ax.set_yticklabels([])


plt.show()