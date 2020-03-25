__author__ = 'Jie'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm,lognorm
import random
from sympy import *

class FatigueCrackGrowth():
    '''
    Three different reliability methods are used for an optimization problem with respect to fitigue crack growth optimization.
    the failure probability is solved.
    '''
    def __init__(self,mu_xiegama=100,mu_a=0.01,mu_Kic=60,mu_c=1.2e-10,std_xiegama=10,std_a=0.005,std_Kic=6,std_c=1.2e-11,m=3.32,count=100000,Nc=3000):
        self.mu_xiegama=mu_xiegama
        self.mu_a=mu_a
        self.mu_Kic=mu_Kic
        self.mu_c=mu_c
        self.std_xiegama=std_xiegama
        self.std_a=std_a
        self.std_Kic=std_Kic
        self.std_c=std_c
        self.m=m
        self.count=count
        self.Nc=Nc

    def lognormToGauss(self,mu_x,std_x,x1):
        # x1 is the point that needed be transformed.
        std_y=np.sqrt(np.log(1+(std_x/mu_x)**2))  # x is lognormal distribution; y=lnx, is normal distribution
        mu_y=np.log(mu_x)-0.5*std_y**2
        dist=lognorm(s=std_y,loc=0,scale=np.exp(mu_y))
        # y=1/(np.sqrt(2*np.pi)*x*xiegama_y)*np.exp(-0.5*((np.log(x)-mu_y)/xiegama_y)**2)
        y1=dist.pdf(x1)                      # the value of y1 in the original distribution
        fai_u1= (np.log(x1)-mu_y)/std_y     # equivalent U1 in the standard normal distribution.
        pdf_u1=norm().pdf(fai_u1)            # the probability density of U1 in standard normal distribution
        std_eqv=pdf_u1/y1
        mu_eqv= x1-fai_u1*std_eqv
        return (mu_eqv,std_eqv)

    def extremIToGauss(self,alpha,delta,x1):
        # this is the transformation of points in a TYpe I Extreme value distribution
            pdf_x1=alpha*np.exp(-(x1-delta)*alpha-np.exp(-(x1-delta)*alpha))
            cdf_x1=np.exp(-np.exp(-alpha*(x1-delta)))
            cdf_inv=norm().ppf(cdf_x1)  #  compute the inverse value of CDF in the standard normal distribution, the cdf_P is used as the dependent value
            pdf_inv=norm().pdf(cdf_inv) # compute the pdf of value cdf_inv in the standard normal distribution.
            std_eqv=pdf_inv/pdf_x1
            mu_eqv=x1-cdf_inv*std_eqv
            return (mu_eqv,std_eqv)
    def fun_opt(self):
        pass

    def fun_opt_grad(self):
        pass