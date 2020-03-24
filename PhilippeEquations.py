__author__ = 'Jie'
# coding=utf8
'''
Calculate pipe properties based on the input file “pipeParameters_philippe.txt”.
Three types of buckling shapes of eigen values can be output,
including pipe with fixed rotation BC (connect to bottom line),
pipe with free rotation BC (connect to bottom line),
and pipe with fixed rotation BC ( connect to neutral axis).
'''

import numpy as np
from  matplotlib import pyplot as plt
import os, sys
from  Property_CrossSection import CSproperty
#from pipeCalculation import PipeProp
import math
from scipy.optimize import fsolve,root
class PipeData():
    '''
    this script is used for obtaining pipe data for calculations based on the file: "pipeParameters_philippe.txt"
    '''
    def philippeData(self):
        '''
        the data are only used for the analytical solution via Philippe's paper. The buckling length is fixed.
        '''
        with open("pipeParameters_philippe.txt","r") as philippeFile:
            data=[]
            for eachline in philippeFile:
                ele=eachline.split(':')
                data.append(ele[1].strip())
            od,t,L,L_b,k,E,alphaT=list(map(lambda x:float(x), data))
        return  od,t,L,L_b,k,E,alphaT

class PipeProp():
    '''
    This script calculates all the propertis of pipe such as I
    xiegamaE,xiegamaU,elonM=[float(el)  for el in input(" please input the three basic material parameters including
    xiegamaE,xiegamaU,ibxilongU:").split()] # MPa
    alpha is the thermal expansion coefficient for pipe steel
    E is the Young's modulus
    '''
    def __init__(self,L,L0,w0m,D,t,E,alpha,deltaT,deltaP):
        self.L=L
        self.L0=L0
        self.w0m=w0m
        self.D=D
        self.t=t
        self.E=E
        self.alpha=alpha
        self.deltaT=deltaT
        self.deltaP=deltaP

    def pipePro(self):
        d=self.D-2*self.t
        alpha1=d/self.D
        I=np.pi*self.D**4*(1-alpha1**4)/64
        J=2*I
        A=(1-alpha1**2)*np.pi*self.D**2/4
        ibxilong_temp=self.alpha*self.deltaT
        ibxilong_pressure=self.deltaP*self.D*(1-2*0.3)/(4*self.t*self.E)
        P0=self.E*A*(ibxilong_temp+ibxilong_pressure)
        L_expand=self.L*(ibxilong_temp+ibxilong_pressure)

        return I,J,A,P0,L_expand

class PipeTorsionPhilippe():
    '''
    the dimensionless parameter is defined via Philippe's paper.
    default data from Philippe: od=324mm,t=17mm,L_b=L=100m,stiffness density (k, N/mm2) : 5248e-6....
    '''
    def __init__(self,D,t,L_b,k,E,alpha_t,A,J,I):
        self.L_b=L_b
        self.D=D
        self.E=E
        self.k=k
        self.alpha_t=alpha_t
        self.J=J
        self.I=I
        self.A=A

        self.G=E/2.0/(1+0.3)
        self.alpha=I/A/(D**2/4)
        self.beta=self.G/E
        self.gama=A*L_b**2/I
        self.xi=k*(D**2/4)/E/A

    def dimensionalData(self):
        return self.G,self.alpha,self.beta,self.gama,self.xi

    def resultsNoTorsion_fun(self,n):
        '''
        case with simplify-supported bc with fixed torsion at the ends
        this is the derivative of Pc_norm for calculation of the limit points.
        '''
        #Pc_norml=n**2*np.pi*2/self.gama+self.alpha*self.gama*self.xi/(n**2*np.pi**2+self.gama*self.xi/2/self.beta)
        driv_Pc=2*n*np.pi**2/self.gama-self.alpha*self.gama*self.xi*2*n*np.pi**2/(n**2*np.pi**2+self.gama*self.xi/2/self.beta)**2
        return driv_Pc

    def resultsNoTorsion_cal(self,n):
        '''
        case with simplify-supported bc with fixed torsion at the ends
        this is the derivative of Pc_norm for calculation of the limit points.
        '''
        Pc_norml=n**2*np.pi**2/self.gama+self.alpha*self.gama*self.xi/(n**2*np.pi**2+self.gama*self.xi/2/self.beta)
        return Pc_norml

    def find_Pc_noTorsion(self):
        # n=fsolve(self.resultsNoTorsion_fun,guess_value)
        # r=root(self.resultsNoTorsion_fun,guess_value)
        # print ("The solving of P0 is: {}".format(r.success))
        # n=np.round(n)
        # Pc_norml=n**2*np.pi*2/self.gama+self.alpha*self.gama*self.xi/(n**2*np.pi**2+self.gama*self.xi/2/self.beta)
        n=np.arange(1,10)
        Pc_all=[]
        for i in n:
            Pc_temp=self.resultsNoTorsion_cal(i)
            Pc_all.append(Pc_temp)
        Pc_min=min(Pc_all)
        n_min=Pc_all.index(Pc_min)+1

        x=np.linspace(0,1,num=1000)
        w=np.sin(n_min*np.pi*x)

        fig,ax=plt.subplots()
        plt.plot(x,w,'-k',linewidth=3,label="Buckling shape with fixed torsion bottomline connection")
        ax=plt.gca()
        ax.set_xlabel(xlabel="The normalized pipe length L_b",fontsize=15)
        ax.set_ylabel(ylabel="Normalized lateral displacement (w_n)",fontsize=15)
        # plt.xlim(0,0.00065)
        # plt.ylim(0,)
        plt.legend(loc='best',prop = {'size':10})
        return Pc_min,n_min

    def resultsNoTorsionNeutral_fun(self,n):
        '''
        case with simplify-supported bc with fixed torsion at the ends, connection with the neutral axis
        this is the derivative of Pc_norm for calculation of the limit points.
        '''
        driv_Pc=2*n*np.pi**2/self.gama-2*self.alpha*self.gama*self.xi/n**3/np.pi**2
        return driv_Pc

    def resultsNoTorsionNeutral_cal(self,n):
        '''
        case with simplify-supported bc with fixed torsion at the ends, connection with the neutral axis
        '''
        temp1=n**2*np.pi**2/self.gama
        temp2=self.alpha*self.gama*self.xi/n**2/np.pi**2
        Pc_norml=temp1+temp2
        return Pc_norml

    def find_Pc_noTorsion_neutral(self):
        # n=fsolve(self.resultsNoTorsionNeutral_fun,guess_value)
        # r=root(self.resultsNoTorsionNeutral_fun,guess_value)
        # print ("The solving of P0 is: {}".format(r.success))
        # n=np.round(n)
        # Pc_norml=n**2*np.pi*2/self.gama+self.alpha*self.gama*self.xi/n**2/np.pi**2
        n=np.arange(1,10)
        Pc_all=[]
        for i in n:
            Pc_temp=self.resultsNoTorsionNeutral_cal(i)
            Pc_all.append(Pc_temp)
        #print (Pc_all)
        Pc_min=min(Pc_all)
        n_min=Pc_all.index(Pc_min)+1

        x=np.linspace(0,1,num=1000)
        w=np.sin(n_min*np.pi*x)

        fig,ax=plt.subplots()
        plt.plot(x,w,'-k',linewidth=3,label="Buckling shape with fixed torsion Neutral connection")
        ax=plt.gca()
        ax.set_xlabel(xlabel="The normalized pipe length L_b",fontsize=15)
        ax.set_ylabel(ylabel="Normalized lateral displacement (w_n)",fontsize=15)
        # plt.xlim(0,0.00065)
        # plt.ylim(0,)
        plt.legend(loc='best',prop = {'size':10})
        return Pc_min,n_min

    def resultsTorsion_fun(self,Pcn):

        a1=np.sqrt(4*self.gama*np.sqrt(2*self.xi*(2*self.alpha*self.beta-Pcn)/self.beta)+self.gama*(self.xi-2*self.beta*Pcn)/self.beta)/4
        b1=np.sqrt(4*self.gama*np.sqrt(2*self.xi*(2*self.alpha*self.beta-Pcn)/self.beta)-self.gama*(self.xi-2*self.beta*Pcn)/self.beta)/4

        a=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))+self.xi-2*self.beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))-self.xi+2*self.beta*Pcn)
        func=b*(-self.alpha*self.xi*self.gama**2+Pcn*self.gama*(a**2+b**2)+3*a**4+2*(a**2)*(b**2)-b**4)*np.sinh(a)+\
             a*(self.alpha*self.xi*self.gama**2+Pcn*self.gama*(a**2+b**2)+a**4-2*(a**2)*(b**2)-3*b**4)*np.sin(b)
        return np.array(func)

    def db_ab(self,Pcn):
        '''
        This function is only used for the coefficients a, b comparison.
        '''
        a=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))+self.xi-2*self.beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))-self.xi+2*self.beta*Pcn)
        a1=np.sqrt(4*self.gama*np.sqrt(2*self.xi*(2*self.alpha*self.beta-Pcn)/self.beta)+self.gama*(self.xi-2*self.beta*Pcn)/self.beta)/4
        b1=np.sqrt(4*self.gama*np.sqrt(2*self.xi*(2*self.alpha*self.beta-Pcn)/self.beta)-self.gama*(self.xi-2*self.beta*Pcn)/self.beta)/4
        return a,a1,b,b1

    def find_Pc_freeTorsion(self,guess_value):

        Pcn=fsolve(self.resultsTorsion_fun,guess_value)
        r=root(self.resultsTorsion_fun,guess_value)
        print ("The solving of P0 is: {}".format(r.success))
        a=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))+self.xi-2*self.beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))-self.xi+2*self.beta*Pcn)

        temp1=(2*a*b*np.cosh(a/2)*np.cos(b/2)+(a**2-b**2)*np.sinh(a/2)*np.sin(b/2))/(2*a*b*(np.cosh(a/2)**2*np.cos(b/2)**2+np.sinh(a/2)**2*np.sin(b/2)**2))
        temp2=((a**2-b**2)*np.cosh(a/2)*np.cos(b/2)-2*a*b*np.sinh(a/2)*np.sin(b/2))/(2*a*b*(np.cosh(a/2)**2*np.cos(b/2)**2+np.sinh(a/2)**2*np.sin(b/2)**2))

        x=np.linspace(-0.5,0.5,num=1000)
        w=1-temp1*np.cosh(a*x)*np.cos(b*x)+temp2*np.sinh(a*x)*np.sin(b*x)

        # fig1,ax1=plt.subplots()
        plt.figure(figsize=(8,6))
        plt.plot(x,w,'-k',linewidth=3,label="K={}$N/mm^2$".format(self.k))

        plt.scatter(0.321,1.305,marker='8',linewidths=10)  # the largest limit point on curves. only for Mode 6.

        ax=plt.gca()
        ax.annotate('(0.321,1.305)', xy=(0.321,1.305), xytext=(0.25,0.5),arrowprops=dict(arrowstyle='simple',facecolor='blue'))
        ax.spines['left'].set_position('center')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_xlabel("Normalized pipe buckle length $L_b$",fontsize=10)
        ax.set_ylabel("Normalized lateral displacement $Y$",fontsize=10)
        ax.set_title("Eigen mode of the case with free torsion")
        ax.axis([-0.5,0.5,0,1.8])
        # plt.xlim(0,0.00065)
        # plt.ylim(0,)
        plt.legend(loc='best',frameon=False,fontsize=10) #prop = {'size':10},
        return Pcn, temp1,temp2


def main ():
    pipeData=PipeData()
    od,t,L,L_b,k,E,alphaT=pipeData.philippeData()
    L0,w0m=0,0
    deltaT,deltaP=0,0
    print ("######################################################################################",flush=True)
    print ("Inputs pipe data for calculation (od,t,L,L_b,k,E,alphaT): \n",od,t,L,L_b,k,E,alphaT)
    print ("######################################################################################\n",flush=True)

    pipeProp=PipeProp(L,L0,w0m,od,t,E,alphaT,deltaT,deltaP)
    I,J,A,P0,L_expand=pipeProp.pipePro()
    pipeTorsionPhilippe=PipeTorsionPhilippe(od,t,L_b,k,E,alphaT,A,J,I)
    G,alpha,beta,gama,xi=pipeTorsionPhilippe.dimensionalData()

    print ("######################################################################################",flush=True)
    print ("the flexural second moment of inertia I is: {}mm4\n"
           "the torsional second moment of inertia J is:{}mm4\n"
           "the cross-section area of steel A is:{}mm2\n".format(I,J,A))

    print ("the calculated dimensionless parameters is:",end="",flush=True)
    print ("\nthe G:{}\nthe alpha:{}"
           "\nthe beta:{}\nthe gama:{}\nthe xi:{}".format(G,alpha,beta,gama,xi))
    print ("##################################################################################\n",flush=True)

    print ("######################################################################################",flush=True)
    print ("the calculated results for pipe with free torsion on ends on the bottom line is:",end="",flush=True)
    Pc_norm, temp1,temp2=pipeTorsionPhilippe.find_Pc_freeTorsion(0.00001)
    print ("Pc_norm is: {}".format(Pc_norm))
    Pc=Pc_norm*E*A
    print ("Pc is: {}".format(Pc))
    a,a1,b,b1=pipeTorsionPhilippe.db_ab(Pc_norm)
    print ("The calculated value of a is:{}".format(a))
    # print (a1)
    print ("The calculated value of b is: {}".format(b))
    # print (b1)
    print (temp1,temp2)
    print ("###################################################################################\n",flush=True)

    print ("######################################################################################",flush=True)
    print ("the calculated results for pipe with fixed torsion on ends on the bottom line is:",end="",flush=True)
    Pc_min,n_min=pipeTorsionPhilippe.find_Pc_noTorsion()
    print ("Pc_norm is : {}".format(Pc_min))
    Pc_min_1=Pc_min*E*A
    print ("Pc is: {}".format(Pc_min_1))
    print ("The corresponding half-wavelength is: {}".format(n_min))
    print ("###################################################################################\n",flush=True)

    print ("######################################################################################",flush=True)
    print ("the calculated results for pipe with fixed torsion on ends Neutral connection is:",end="",flush=True)
    Pc_min,n_min=pipeTorsionPhilippe.find_Pc_noTorsion_neutral()
    print ("Pc_norm is : {}".format(Pc_min))
    Pc_min_2=Pc_min*E*A
    print ("Pc is: {}".format(Pc_min_2))
    print ("The corresponding half-wavelength is: {}".format(n_min))
    print ("###################################################################################\n",flush=True)
    plt.show()

if __name__ == '__main__':
    main()