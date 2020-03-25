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
    def __init__(self,mu_xiegama=100,mu_a=0.01,mu_Kic=60,mu_c=1.2e-10,std_xiegama=10,std_a=0.005,std_Kic=6,std_c=1.2e-11,m=3.32,count=100000,Nc=4000):
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
        dist=lognorm(s=std_y,loc=0,scale=np.exp(mu_y))  # there is a transform between the lognorm function and the original expression
        # y=1/(np.sqrt(2*np.pi)*x*xiegama_y)*np.exp(-0.5*((np.log(x)-mu_y)/xiegama_y)**2)
        y1=dist.pdf(x1)                      # the value of y1 in the original distribution
        fai_inv= (np.log(x1)-mu_y)/std_y      # equivalent U1 in the standard normal distribution.
        pdf_faiinv=norm().pdf(fai_inv)            # the probability density of U1 in standard normal distribution
        std_eqv=pdf_faiinv/y1
        mu_eqv= x1-fai_inv*std_eqv
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
    def fun_opt(self,Kic,xiegama,ai,c):
        af=(Kic/(1.1215*xiegama))**2/np.pi
        frac=af**(1-self.m/2)-ai**(1-self.m/2)
        denomi=c*(1.1215*xiegama)**self.m*np.pi**(self.m/2)*(1-self.m/2)
        limit_state=self.Nc-frac/denomi
        return limit_state

    def limitState_1stDeriv(self):
        ## return the expression of first-order derivative, and the original expression
        xiegama,Kic,ai,c=Symbol('xiegama'),Symbol('Kic'),Symbol('ai'),Symbol('c')
        af=(Kic/(1.1215*xiegama))**2/np.pi
        frac=af**(1-self.m/2)-ai**(1-self.m/2)
        denomi=c*(1.1215*xiegama)**self.m*np.pi**(self.m/2)*(1-self.m/2)
        Nf=frac/denomi
        N_lim=self.Nc-Nf
        xiegama_prime=N_lim.diff(xiegama)
        Kic_prime=N_lim.diff(Kic)
        ai_prime=N_lim.diff(ai)
        c_prime=N_lim.diff(c)
        return (N_lim,xiegama_prime,Kic_prime,ai_prime,c_prime)

    def limitState_1stDeriv_value(self,Kic,xiegama,ai,c):
        # the values of the first order of limit state expression and the original expression
        N_lim,xiegama_prime,Kic_prime,ai_prime,c_prime=self.limitState_1stDeriv()
        N_lim_value=N_lim.evalf(subs={"Kic":Kic,'xiegama':xiegama,'ai':ai,'c':c})
        xiegama_prime_value=xiegama_prime.evalf(subs={"Kic":Kic,'xiegama':xiegama,'ai':ai,'c':c})
        Kic_prime_value=Kic_prime.evalf(subs={"Kic":Kic,'xiegama':xiegama,'ai':ai,'c':c})
        ai_prime_value=ai_prime.evalf(subs={"Kic":Kic,'xiegama':xiegama,'ai':ai,'c':c})
        c_prime_value=c_prime.evalf(subs={"Kic":Kic,'xiegama':xiegama,'ai':ai,'c':c})
        return (N_lim_value,xiegama_prime_value ,Kic_prime_value,ai_prime_value,c_prime_value)

    def fun_opt_grad(self):
        pass

    def FORM_HLRF(self):
        ### the first iteration, by using the mean value point as the inital design point.
        ## obtain the new design point X2. The first beta calculation is due to original definition

        #### transform the equivelant mus and xiegamas in the standard normal distribution.
        mu_xiegama_eqv,std_xiegama_eqv=self.lognormToGauss(self.mu_xiegama,self.std_xiegama,self.mu_xiegama)
        mu_Kic_eqv,std_Kic_eqv=self.lognormToGauss(self.mu_Kic,self.std_Kic,self.mu_Kic)
        mu_ai_eqv,std_ai_eqv=self.lognormToGauss(self.mu_a,self.std_a,self.mu_a)
        mu_c_eqv,std_c_eqv=self.lognormToGauss(self.mu_c,self.std_c,self.mu_c)

        # compute the first iteration
        mu_g=self.fun_opt(self.mu_Kic,self.mu_xiegama,self.mu_a,self.mu_c)  # the mu_P should not be updated in the first iteration
        N_lim_value,xiegama_prime_value ,Kic_prime_value,ai_prime_value,c_prime_value=self.limitState_1stDeriv_value(self.mu_Kic,self.mu_xiegama,self.mu_a,self.mu_c)
        N_lim_value,xiegama_prime_value ,Kic_prime_value,ai_prime_value,c_prime_value=float(N_lim_value),float(xiegama_prime_value),\
                                                                                      float(Kic_prime_value),float(ai_prime_value),float(c_prime_value)
        std_g=np.sqrt(xiegama_prime_value**2+Kic_prime_value**2+ai_prime_value**2+c_prime_value**2)

        beta_init=mu_g/std_g
        cos_alpha1=-xiegama_prime_value*self.std_xiegama/std_g
        cos_alpha2=-Kic_prime_value*self.std_Kic/std_g
        cos_alpha3=-ai_prime_value*self.std_a/std_g
        cos_alpha4=-c_prime_value*self.std_c/std_g

        xiegama_init=self.mu_xiegama+beta_init*self.std_xiegama*cos_alpha1
        Kic_init=self.mu_Kic+beta_init*self.std_Kic*cos_alpha2
        ai_init=self.mu_a+beta_init*self.std_a*cos_alpha3
        c_init=self.mu_c+beta_init*self.std_c*cos_alpha4
        uxiegama_init=(xiegama_init-self.mu_xiegama)/self.std_xiegama
        uKic_init=(Kic_init-self.mu_Kic)/self.std_Kic
        uai_init=(ai_init-self.mu_a)/self.std_a
        uc_init=(c_init-self.mu_c)/self.std_c

        ## loop from the new design point X2.
        ibxilong=0.001
        i=1
        flag=True

        while flag:
            ### correct the non-normal distribution as previous procedures

            mu_xiegama_eqv,std_xiegama_eqv=self.lognormToGauss(self.mu_xiegama,self.std_xiegama,xiegama_init)
            mu_Kic_eqv,std_Kic_eqv=self.lognormToGauss(self.mu_Kic,self.std_Kic,Kic_init)
            mu_ai_eqv,std_ai_eqv=self.lognormToGauss(self.mu_a,self.std_a,ai_init)
            mu_c_eqv,std_c_eqv=self.lognormToGauss(self.mu_c,self.std_c,c_init)

            # interation starts from the second step

            g_new=self.fun_opt(Kic_init,xiegama_init,ai_init,c_init)

            N_lim_value,xiegama_prime_value ,Kic_prime_value,ai_prime_value,c_prime_value=self.limitState_1stDeriv_value(Kic_init,xiegama_init,ai_init,c_init)
            N_lim_value,xiegama_prime_value ,Kic_prime_value,ai_prime_value,c_prime_value=float(N_lim_value),float(xiegama_prime_value),\
                                                                                      float(Kic_prime_value),float(ai_prime_value),float(c_prime_value)
            std_g=np.sqrt(xiegama_prime_value**2+Kic_prime_value**2+ai_prime_value**2+c_prime_value**2)

            beta_new=(g_new-xiegama_prime_value*std_xiegama_eqv*uxiegama_init-Kic_prime_value*std_Kic_eqv*uKic_init-
                      ai_prime_value*std_ai_eqv*uai_init-c_prime_value*std_c_eqv*uc_init)/std_g  # noted: the std value is updated.
            cos_alpha1=-xiegama_prime_value*std_xiegama_eqv/std_g
            cos_alpha2=-Kic_prime_value*std_Kic_eqv/std_g
            cos_alpha3=-ai_prime_value*std_ai_eqv/std_g
            cos_alpha4=-c_prime_value*std_c_eqv/std_g

            xiegama_new=mu_xiegama_eqv+beta_new*std_xiegama_eqv*cos_alpha1
            Kic_new=mu_Kic_eqv+beta_init*std_Kic_eqv*cos_alpha2
            ai_new=mu_ai_eqv+beta_init*std_ai_eqv*cos_alpha3
            c_new=mu_c_eqv+beta_init*std_c_eqv*cos_alpha4
            uxiegama_new=(xiegama_init-mu_xiegama_eqv)/std_xiegama_eqv    # update the X value in U-space
            uKic_new=(Kic_init-mu_Kic_eqv)/std_Kic_eqv
            uai_new=(ai_init-mu_ai_eqv)/std_ai_eqv
            uc_new=(c_init-mu_c_eqv)/std_c_eqv

            residual=np.abs((beta_new-beta_init)/beta_init)
            if residual<=ibxilong:
                flag=False
            beta_init=beta_new
            xiegama_init=xiegama_new
            Kic_init=Kic_new
            ai_init=ai_new
            c_init=c_new
            uxiegama_init=uxiegama_new
            uKic_init=uKic_new
            uai_init=uai_new
            uc_init=uc_new
            g_test=g_new
            i+=1
            print ("loop:{}".format(i))
            print ("g_new: {}".format(g_new))
            print ("the safety factor beta: {} with a residual error {}".format(beta_new,residual))
        print ("solution is converged")
        gauss=norm(loc=0,scale=1)
        Pf=gauss.cdf(-beta_new)
        return (Pf,beta_new,residual)

fatigueCrackGrowth=FatigueCrackGrowth()
Pf,beta_new,residual=fatigueCrackGrowth.FORM_HLRF()
print (Pf)