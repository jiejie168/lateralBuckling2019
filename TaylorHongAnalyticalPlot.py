__author__ = 'Jie'
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve,root
from pandas import DataFrame,Series

'''
This script is only used for the comparison of analytical solutions for pipes with model imperfection,
but buckled in modes in terms of Mode1, and Mode 3.
The analytical solutions are from Taylor, 1986, and Hong, 2015.
The solutions from Taylor is validated based on his own data. The solutions from HONG do not match with his
literature. It could because of the error in his inputs parameters. Based on the comparison results with Taylor, all
the scripts should be correct.

Besides, the effect of varying axial friction is studied through comparison. This case is only for the
pipes with Mode 1 imperfection, and Mode 1 buckling. It is possible to extend to the case for pipes with
Mode 1 imperfection, but Mode 3 buckling.  13/09/2019
'''

f_L,q,I,A,E,D,alpha=1.0,3.8,1.509e9,29920,2.06e5,650,1.1e-5  # original data from Taylor,1986

#f_L,q,I,A,E,D,t,alpha=0.4,3.8,1.506e8,12416,2.06e5,323.9,12.7,1.1e-5  # original data from Liu, 2014
u_fai=5
f_A=0.8
#L=125000

def buckleProTaylorImpf1(impfRatio,L,L0):
    '''
    ### this analytical solution is from Taylor,1986. imperfection mode 1, buckling in Mode 1 as well.
    #first consider the axial frictional force is constant.
    :param impfRatio:
    :param L:
    :param L0:
    :return:
    '''
    beta=L0/L
    R1=4.60314*(np.sin(4.4934*beta)+2.30157*(np.sin(4.4934*(1+beta))/(1/beta+1)+np.sin(4.4934*(1-beta))/(1/beta-1)))
    P=80.76*(E*I/L**2)*(1-(R1/75.60)*(beta)**2)    # buckling force within buckle
    P_e=np.pi**2*E*I/L**2*(1+q*f_L*L**3/np.pi**4/E/I)
    wm=2.407e-3*q*f_L*L**4/E/I
    w0m=impfRatio*L0

    temp1=(f_A*q*L/2)**2
    temp2=6.39064e-5*f_A*f_L**2*q**3*A*(L**7-L0**7)/E/I/I
    P0=P-f_A*q*L/2+0.5*(temp1+temp2)**0.5
    deltaT=P0/(A*E*alpha)
    return P,P0,deltaT,P_e

def buckleProTaylorImpf1_1(L,wm,L0,w0m):
    '''
    ### this analytical solution is from Taylor,1986. imperfection mode 1, buckling in Mode 1 as well.
    #first consider the axial frictional force is constant.
    The same with the former one, but with different inputs parameters.
    :param impfRatio:
    :param L:
    :param L0:
    :return:
    '''
    beta=L0/L
    R1=4.60314*(np.sin(4.4934*beta)+2.30157*(np.sin(4.4934*(1+beta))/(1/beta+1)+np.sin(4.4934*(1-beta))/(1/beta-1)))
    P=80.76*(E*I/L**2)*(1-(R1/75.60)*(beta)**2)    # buckling force within buckle
    P_e=np.pi**2*E*I/L**2*(1+q*f_L*L**3/np.pi**4/E/I)

    temp1=(f_A*q*L/2)**2
    temp2=6.39064e-5*f_A*f_L**2*q**3*A*(L**7-L0**7)/E/I/I
    P0=P-f_A*q*L/2+0.5*(temp1+temp2)**0.5
    deltaT=P0/(A*E*alpha)
    return P,P0,deltaT,P_e

# def func(x):    # nonlinear function for numerical solutions
#     u_s,P0=x[0],x[1]
#     return np.array([u_s-(P0-P)*L/(2*A*E)+7.9883e-6*
#                      (f_L*q/E/I)**2*(L**7-L0**7),
#     P0-P-(2*f_A*q*A*E*((-np.exp(25*u_s/u_fai)+1)/5-u_s))**0.5])

def buckleProTaylorImpf1_withVariationAxialFriction(impfRatio,L,L0,us_initial,P0_initial):
        '''
        this analytical solution is from Taylor,1986. imperfection mode 1, buckling in Mode 1 as well.
        #The axial frictional coeficient f_A here is considered as variation along the slide length.
        The function here provides only the specific results for each each. For variation of buckling length, please
        refer to the script in Jupyter notebook: taylorAnalytical-test.

        specified parameters: L (the pipe length),, impfRatio (the imperfection amplitude to length ratio)
        output parameters: P(buckling force within buckle),P0 (buckling force outside buckle),deltaT,
                      wm(the amplitude of pipe with L), xiegamaM (within buckle)
        tips: the initial guess for the first 100 iteration (when buckling length approching imperfection length), the
        initial u_s should be very small, such as -0.01. Then it should be increase to -50 or something
        '''
        beta=L0/L
        R1=4.60314*(np.sin(4.4934*beta)+2.30157*(np.sin(4.4934*(1+beta))/(1/beta+1)+np.sin(4.4934*(1-beta))/(1/beta-1)))
        P=80.76*(E*I/L**2)*(1-(R1/75.60)*(beta)**2)    # buckling force within buckle
        wm=2.407e-3*q*f_L*L**4/E/I
        w0m=impfRatio*L0
        u_fai=5
        def func(x):    # nonlinear function for numerical solutions
            u_s,P0=x[0],x[1]
            return np.array([u_s-(P0-P)*L/(2*A*E)+7.9883e-6*
                             (f_L*q/E/I)**2*(L**7-L0**7),
            P0-P-(2*f_A*q*A*E*((-np.exp(25*u_s/u_fai)+1)/5-u_s))**0.5])
        u_s,P0=fsolve(func,[us_initial,P0_initial])  #initial guess of solution is -50,1e6,
        r=root(func,[us_initial,P0_initial])         # another way for solutions, here only used for solution indications.
        print ("The solving of P0 is: {}".format(r.success))

        deltaT=P0/(A*E*alpha)
        M_m=-0.06938*f_L*q*(L**2-L0**2)
        xiegamaM=P/A+np.abs(M_m*D/2/I)
        u_s=(P0-P)*L/(2*A*E)-7.9883e-6*(f_L*q/E/I)**2*(L**7-L0**7)
        return P,P0,deltaT,wm,xiegamaM,u_s

def buckleProHongImpf1(impfRatio,L,L0):
    '''
    from Hong's, 2015 analytical solution for pipe with mode 3 imperfection.
    mode 1 imperfection, with mode 3 final buckling.
    :param impfRatio:
    :param L:  Be careful, for bulking mode 3 in Hong's,2015, this length is only the center part of the pipelines.
    the entire length of pipe is 2.59*L
    :param L0:
    :return:
    '''
    # impfRatio=w0m/L0
    #L0=np.power((impfRatio)*E*I/1.032e-2/f_L/q,1/3)
    beta=L0/L   # the L for mode 3 buckling is the length of arc in the center.
    R3=7.57*beta+7.81*np.sin(2.915*beta)+15.25*(np.sin(2.915*(1+beta))/(1+1/beta)+np.sin(2.915*(1-beta))/(-1+1/beta))
    P=33.99*(E*I/L**2)*(1-(R3/154.84)*(beta)**2)    # buckling force within buckle
    wm=1.032e-2*q*f_L*L**4/E/I
    w0m=impfRatio*L0

    temp1=(3.6213e-4*L**5-1.4556e-4*(impfRatio/wm/L)**(2/3)*(E*I*impfRatio/1.032e-2/f_L/q)**(5/3))
    temp2=(1+f_L*q*f_L*A*temp1/f_A/E/I/I)**0.5
    P0=33.99*(E*I/L**2)*(1-(R3/154.84)*(impfRatio/wm/L)**(2/3))+1.3*q*f_A*L*((temp2)-1)
    deltaT=P0/(A*E*alpha)
    u_s=1.3*(P0-P)*L/E/A-3.06e-4*(q*f_L/E/I)**2*L**7+1.23e-4*(q*f_L/E/I)**2*L0**7
    return P,P0,deltaT,wm,u_s

'''
#######################################################################################################################
## This part is simply to calculate the P,P0,deltaT using Taylor's,1986 solutions for individual impfRatio.
## Inputs: impfRatio,L,L0
## Outputs:  P,P0,deltaT
#######################################################################################################################
#L0=23000
# w0m=150
impfRatio=0.00652
L_max=250000

L0_taylor=np.power(impfRatio*E*I/2.407e-3/f_L/q,1.0/3)  #L0, w0m is fixed due to impfRatio
w0m=impfRatio*L0_taylor
L_taylor=np.linspace(L0_taylor,L_max,num=200)
#L_taylor=np.arange(L0_taylor,L_max,50)

P,P0,deltaT,P_e=buckleProTaylorImpf1(impfRatio,L_taylor,L0_taylor)
deltaT_are_NANs=np.isnan(deltaT)
P_are_NANs=np.isnan(P)
P0_are_NANs=np.isnan(P0)
deltaT[deltaT_are_NANs]=0
P[P_are_NANs]=0
P0[P0_are_NANs]=0
wm=2.407e-3*f_L*q*L_taylor**4/E/I
M_m=-0.06938*f_L*q*(L_taylor**2-L0_taylor**2)
xiegamaM=P/A+np.abs(M_m*D/2/I)
u_s=(P0-P)*L_taylor/(2*A*E)-7.9883e-6*(f_L*q/E/I)**2*(L_taylor**7-L0_taylor**7)
print ("###################################################################################################")
print ("Inputs pipe data for calculation (od,t,L_max,q,f_L,f_A,E,alpha): \n",
       D,t,L_max,q,f_L,f_A,E,alpha)
print ("Inputs imperfection data for calculation (impfRatio, only works for Taylor,1986 or Hong,2015): \n"
       "{}".format(impfRatio))
print ("################################################################################################\n")

print ("###################################################################################################")
print ("the calculated length of initial imperfection L0 is : {}mm".format(L0_taylor))
print ("the calculated amplitude of initial imperfection w0m is: {}mm".format(w0m))
print ("#################################################################################################\n")

print ("###################################################################################################")
print ("The first 5 calculated buckling length L is:{} mm".format(L_taylor[:10]))
print ("The corresponding buckling amplitude wm is:{} mm".format(wm[:10]))
print ("The corresponding elastic critical buckling force Pe is:{} N ".format(P_e[:10]))
print ("The corresponding Buckling force P is:{} N ".format(P[:10]))
print ("The corresponding preBuckling force P0 is: {} N".format(P0[:10]))
print("The corresponding temp increase deltaT is: {} oC ".format(deltaT[:10]))
print ("The corresponding axial stress xiegamaM is :{} MPa".format(xiegamaM[:10]))
print ("#################################################################################################\n")

print ("##################################################################################################")
print ("The last 5 calculated pipe buckling L is:{} mm".format(L_taylor[-5:]))
print ("The corresponding buckling amplitude wm is:{}mm".format(wm[-5:]))
print ("The corresponding elastic critical buckling force Pe is:{}N ".format(P_e[-5:]))
print ("The corresponding Buckling force P is:{}N ".format(P[-5:]))
print ("The corresponding preBuckling force P0 is: {}N".format(P0[-5:]))
print("The corresponding temp increase deltaT is: {} oC\n ".format(deltaT[-5:]))
print ("The corresponding axial stress xiegamaM is :{} MPa".format(xiegamaM[-5:]))
print ("#################################################################################################\n")

figure1,axes1=plt.subplots()
figure2,axes2=plt.subplots()
figure3,axes3=plt.subplots()
axes1.plot(wm,P,'-Dk',lw=2,label="imperf amplitude"+np.str(w0m),markevery=2,ms=6)
axes2.plot(wm,P0,'-Dk',lw=2,label="imperf amplitude"+np.str(w0m),markevery=2,ms=6)
axes1.set_xlabel(xlabel="Buckling amplitude wm (mm)",fontsize=15)
axes1.set_ylabel(ylabel="Buckling force within buckle P (N)",fontsize=15)

axes2.set_xlabel(xlabel="Buckling amplitude wm (mm)",fontsize=15)
axes2.set_ylabel(ylabel="Pre-buckling force N (N)",fontsize=15)

axes3.plot(wm,deltaT,'-Dk',lw=2,label="imperf amplitude"+np.str(w0m),markevery=2,ms=6)
axes3.set_xlabel(xlabel="Buckling amplitude wm (mm)",fontsize=15)
axes3.set_ylabel(ylabel="Temperature difference (oC)",fontsize=15)

axes1.set_xlim(0,L_max/100)
axes1.set_ylim(0,1.2e6)
axes2.set_xlim(0,L_max/100)
axes2.set_ylim(0,2e6)
axes3.set_xlim(0,L_max/100)
axes3.set_ylim(0,100)

axes1.legend()
axes2.legend()
axes3.legend()
plt.show()

'''
#####################################################################################################################
## Start with a series of initial imperfection with a specific impfRatio, and given the buckling length of pipe,
## the corresponding imperfection length, imperfection amplitude,P, P0, and deltaT are calculated.
## this script is only used for the validation of Taylor, 1986, and comparison with Hong,2015
#####################################################################################################################
ImpfRatio=[0.003,0.007,0.01]  #  the imperfection length and amplitude are then fixed. The pipelength is varied, larger
L_max=125000  # the maximum buckling length of L. for bulking mode 3 in Hong's,2015, the calculated length is only L_max/2.59
# ImpfRatio=[0.003,0.00652,0.01]  # original data from Liu, 2014.
# L_max=250000
num=len(ImpfRatio)
# than the imperfection length L0. Once L is known, the wm is known.
######################################################################################################################
#Calculate the P,P0,wm, and others under specific impfRatios based on Taylor analtyical with axial fricitial variation
######################################################################################################################
deltaT_all_taylor=[]
P_all_taylor=[]
P0_all_taylor=[]
wm_all_taylor=[]
xiegamaM_all_taylor=[]
u_sAll_taylor=[]
L_all_taylor=[]
for impfRatio in ImpfRatio:
    L0_taylor=np.power(impfRatio*E*I/2.407e-3/f_L/q,1.0/3)
    w0m=L0_taylor*impfRatio
    L_taylor=np.arange(L0_taylor,L_max,100)

    P,P0,deltaT,P_e=buckleProTaylorImpf1(impfRatio,L_taylor,L0_taylor)
    deltaT_are_NANs=np.isnan(deltaT)
    P_are_NANs=np.isnan(P)
    P0_are_NANs=np.isnan(P0)
    deltaT[deltaT_are_NANs]=0
    P[P_are_NANs]=0
    P0[P0_are_NANs]=0
    wm=2.407e-3*f_L*q*L_taylor**4/E/I
    M_m=-0.06938*f_L*q*(L_taylor**2-L0_taylor**2)
    xiegamaM=P/A+np.abs(M_m*D/2/I)
    u_s=(P0-P)*L_taylor/(2*A*E)-7.9883e-6*(f_L*q/E/I)**2*(L_taylor**7-L0_taylor**7)

    deltaT_all_taylor.append(deltaT)
    P_all_taylor.append(P)
    P0_all_taylor.append(P0)
    wm_all_taylor.append(wm)
    xiegamaM_all_taylor.append(xiegamaM)
    u_sAll_taylor.append(u_s)
    L_all_taylor.append(L_taylor)

######################################################################################################################
#Calculate the P,P0,wm, and others under specific impfRatios based on Taylor analtyical with axial fricitial variation
######################################################################################################################
# deltaT_all_temp={}
# P_all_temp={}
# P0_all_temp={}
# wm_all_temp={}
# xiegamaM_all_temp={}
# u_sAll_temp={}

deltaT_all_taylor1_1={}
P_all_taylor1_1={}
P0_all_taylor1_1={}
wm_all_taylor1_1={}
xiegamaM_all_taylor1_1={}
u_sAll_taylor1_1={}
L_all_taylor1_1={}
for i, impfRatio in enumerate(ImpfRatio):
    L0_taylor1=np.power(impfRatio*E*I/2.407e-3/f_L/q,1.0/3)
    w0m=L0_taylor1*impfRatio
    L_taylor1=np.arange(L0_taylor1,L_max,100)

    for num_iter,element in enumerate(L_taylor1):
        if num_iter<=200:
            us_initial=-0.01
            P0_initial=4e5
            P,P0,deltaT,wm,xiegamaM,u_s=buckleProTaylorImpf1_withVariationAxialFriction(impfRatio,element,
                                                                                    L0_taylor1,us_initial,P0_initial)
        else:
            us_initial=-50
            P0_initial=2e6
            P,P0,deltaT,wm,xiegamaM,u_s=buckleProTaylorImpf1_withVariationAxialFriction(impfRatio,element,L0_taylor1,
                                                                                        us_initial,P0_initial)

        deltaT_all_taylor1_1.setdefault(i,[]).append(deltaT)
        P_all_taylor1_1.setdefault(i,[]).append(P)
        P0_all_taylor1_1.setdefault(i,[]).append(P0)
        wm_all_taylor1_1.setdefault(i,[]).append(wm)
        xiegamaM_all_taylor1_1.setdefault(i,[]).append(xiegamaM)
        u_sAll_taylor1_1.setdefault(i,[]).append(u_s)
        L_all_taylor1_1.setdefault(i,[]).append(element)

deltaT_all_taylor1=[]
P_all_taylor1=[]
P0_all_taylor1=[]
wm_all_taylor1=[]
xiegamaM_all_taylor1=[]
u_sAll_taylor1=[]
L_all_taylor1=[]

for i in range(num):
    temp1_deltaT=Series(deltaT_all_taylor1_1[i])
    temp2_P=Series(P_all_taylor1_1[i])
    temp3_P0=Series(P0_all_taylor1_1[i])
    temp4_wm=Series(wm_all_taylor1_1[i])
    temp5_xiegamaM=Series(xiegamaM_all_taylor1_1[i])
    temp6_us=Series(u_sAll_taylor1_1[i])
    temp7_L=Series(L_all_taylor1_1[i])

    temp1_deltaT=temp1_deltaT.values
    temp2_P=temp2_P.values
    temp3_P0=temp3_P0.values
    temp4_wm=temp4_wm.values
    temp5_xiegamaM=temp5_xiegamaM.values
    temp6_us=temp6_us.values
    temp7_L=temp7_L.values

    temp1_are_nans=np.isnan(temp1_deltaT)
    temp2_are_nans=np.isnan(temp2_P)
    temp3_are_nans=np.isnan(temp3_P0)
    temp4_are_nans=np.isnan(temp4_wm)
    temp5_are_nans=np.isnan(temp5_xiegamaM)
    temp6_are_nans=np.isnan(temp6_us)
    temp7_are_nans=np.isnan(temp7_L)

    temp1_deltaT[temp1_are_nans]=0
    temp2_P[temp2_are_nans]=0
    temp3_P0[temp3_are_nans]=0
    temp4_wm[temp4_are_nans]=0
    temp5_xiegamaM[temp5_are_nans]=0
    temp6_us[temp6_are_nans]=0
    temp7_L[temp7_are_nans]=0

    deltaT_all_taylor1.append(temp1_deltaT)
    P_all_taylor1.append(temp2_P)
    P0_all_taylor1.append(temp3_P0)
    wm_all_taylor1.append(temp4_wm)
    xiegamaM_all_taylor1.append(temp5_xiegamaM)
    u_sAll_taylor1.append(temp6_us)
    L_all_taylor1.append(temp7_L)

#print (L0_taylor1)
print (P0_all_taylor1[0][:10],"\n")
print (P_all_taylor1[0][:10],"\n")
print (L_all_taylor1[0][:10],"\n")
print (u_sAll_taylor1[0][:10],"\n")
#print (deltaT_all_taylor1[0][:100],"\n")

######################################################################################################################
#Calculate the P,P0,wm, and others under specific impfRatios based on Hong's analtyical without axial fricitial variation
######################################################################################################################
deltaT_all_hong=[]
P_all_hong=[]
P0_all_hong=[]
wm_all_hong=[]
xiegamaM_all_hong=[]
u_sAll_hong=[]
L_all_hong=[]
for impfRatio in ImpfRatio:
    L0_hong=np.power((impfRatio)*E*I/1.032e-2/f_L/q,1/3)
    w0m=L0_hong*impfRatio
    #Be careful, for bulking mode 3 in Hong's,2015, this length is only the center part of the pipelines.
    #the entire length of pipe is 2.59*L
    L_hong=np.arange(L0_hong,L_max,100)

    P,P0,deltaT,wm,u_s=buckleProHongImpf1(impfRatio,L_hong,L0_hong)
    deltaT_are_NANs=np.isnan(deltaT)
    P_are_NANs=np.isnan(P)
    P0_are_NANs=np.isnan(P0)
    deltaT[deltaT_are_NANs]=0
    P[P_are_NANs]=0
    P0[P0_are_NANs]=0

    deltaT_all_hong.append(deltaT)
    P_all_hong.append(P)
    P0_all_hong.append(P0)
    wm_all_hong.append(wm)
    u_sAll_hong.append(u_s)
    L_all_hong.append(L_hong)

##############################################################################
##plot the diagrams of different buckling amplitude and pre-buckling force P0.
##############################################################################
fig3,axe3=plt.subplots()
axe3.plot(wm_all_taylor[0],P0_all_taylor[0],'--ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe3.plot(wm_all_taylor[1],P0_all_taylor[1],'--Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe3.plot(wm_all_taylor[2],P0_all_taylor[2],'--8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)
axe3.plot(wm_all_hong[0],P0_all_hong[0],'-ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe3.plot(wm_all_hong[1],P0_all_hong[1],'-Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe3.plot(wm_all_hong[2],P0_all_hong[2],'-8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="Buckling amplitude wm (mm)",fontsize=15)
ax.set_ylabel(ylabel="Pre-buckling force (P0)",fontsize=15)
plt.xlim(0,4000)
plt.ylim(0,1.5e7)
plt.legend()

###########################################################################
##plot the diagrams of different buckling amplitude and temperature increase.
############################################################################

fig,axe1=plt.subplots()
axe1.plot(wm_all_taylor[0],deltaT_all_taylor[0],'--ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe1.plot(wm_all_taylor[1],deltaT_all_taylor[1],'--Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe1.plot(wm_all_taylor[2],deltaT_all_taylor[2],'--8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)
axe1.plot(wm_all_hong[0],deltaT_all_hong[0],'-ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe1.plot(wm_all_hong[1],deltaT_all_hong[1],'-Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe1.plot(wm_all_hong[2],deltaT_all_hong[2],'-8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="wm",fontsize=15)
ax.set_ylabel(ylabel="deltaT",fontsize=15)
plt.xlim(0,4000)
plt.ylim(0,200)
plt.legend()


#####################################################################
##plot the diagrams of different buckling length and buckling force.
#####################################################################
fig2,axe2=plt.subplots()
axe2.plot(L_all_taylor[0],P_all_taylor[0],'--ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe2.plot(L_all_taylor[1],P_all_taylor[1],'--Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe2.plot(L_all_taylor[2],P_all_taylor[2],'--8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)
axe2.plot(L_all_hong[0],P_all_hong[0],'-ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe2.plot(L_all_hong[1],P_all_hong[1],'-Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe2.plot(L_all_hong[2],P_all_hong[2],'-8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="Buckling length L (mm)",fontsize=15)
ax.set_ylabel(ylabel="Buckling force P (N)",fontsize=15)
plt.xlim(0,L_max)
plt.ylim(0,1.0e7)
plt.legend()


#####################################################################
##plot the diagrams of different buckling length and PRE-buckling force.
#####################################################################
fig4,axe4=plt.subplots()
axe4.plot(L_all_taylor[0],P0_all_taylor[0],'--ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe4.plot(L_all_taylor[1],P0_all_taylor[1],'--Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe4.plot(L_all_taylor[2],P0_all_taylor[2],'--8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)
axe4.plot(L_all_hong[0],P0_all_hong[0],'-ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
axe4.plot(L_all_hong[1],P0_all_hong[1],'-Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
axe4.plot(L_all_hong[2],P0_all_hong[2],'-8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="Buckling length L (mm)",fontsize=15)
ax.set_ylabel(ylabel="Pre-buckling force P0 (N)",fontsize=15)
plt.xlim(0,L_max)
plt.ylim(0,1.5e7)
plt.legend()


# ######################################################################################################################
# ##plot the comparison diagrams of pipes with and without axial fricitonal force based on Taylor,1986.
# #######################################################################################################################
#
# fig5,axe5=plt.subplots()
# axe5.plot(L_all_taylor[0],P0_all_taylor[0],'--ok',linewidth=2,label="w0m/L0=0.003",markevery=20,ms=6)
# axe5.plot(L_all_taylor[1],P0_all_taylor[1],'--Dr',linewidth=2,label="w0m/L0=0.007",markevery=20,ms=6)
# axe5.plot(L_all_taylor[2],P0_all_taylor[2],'--8b',linewidth=2,label="w0m/L0=0.01",markevery=20,ms=6)
# axe5.plot(L_all_taylor1[0],P0_all_taylor1[0],'-*b',linewidth=2,label="w0m/L0=0.003",markevery=10,ms=6)
# axe5.plot(L_all_taylor1[1],P0_all_taylor1[1],'-vy',linewidth=2,label="w0m/L0=0.007",markevery=10,ms=6)
# axe5.plot(L_all_taylor1[2],P0_all_taylor1[2],'-^b',linewidth=2,label="w0m/L0=0.01",markevery=10,ms=6)
#
# ax=plt.gca()
# ax.set_xlabel(xlabel="Buckling length L (mm)",fontsize=15)
# ax.set_ylabel(ylabel="Pre-buckling force P0 (N)",fontsize=15)
# plt.xlim(0,L_max)
# plt.ylim(0,1.5e7)
# plt.legend()


fig6,axe6=plt.subplots()
axe6.plot(wm_all_taylor[0],P0_all_taylor[0]/1e6,'--ok',linewidth=2,label="$w_{0m}$/$L_0$=0.003",markevery=20,ms=6)
axe6.plot(wm_all_taylor[1],P0_all_taylor[1]/1e6,'--Dr',linewidth=2,label="$w_{0m}$/$L_0$=0.007",markevery=20,ms=6)
axe6.plot(wm_all_taylor[2],P0_all_taylor[2]/1e6,'--8b',linewidth=2,label="$w_{0m}$/$L_0$=0.01",markevery=20,ms=6)
axe6.plot(wm_all_taylor1[0],P0_all_taylor1[0]/1e6,'-*b',linewidth=2,label="$w_{0m}$/$L_0$=0.003-AxialVar",markevery=10,ms=6)
axe6.plot(wm_all_taylor1[1],P0_all_taylor1[1]/1e6,'-vy',linewidth=2,label="$w_{0m}$/$L_0$=0.007-AxialVar",markevery=10,ms=6)
axe6.plot(wm_all_taylor1[2],P0_all_taylor1[2]/1e6,'-^b',linewidth=2,label="$w_{0m}$/$L_0$=0.01-AxialVar",markevery=10,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="Buckling amplitude $w_m$ (mm)",fontsize=15)
ax.set_ylabel(ylabel="Buckling force $P_0$ (MN)",fontsize=15)
plt.xlim(0,7000)
plt.ylim(0,10)
plt.legend(loc='best',prop = {'size':10})


fig7,axe7=plt.subplots()
plt.plot(wm_all_taylor[0],xiegamaM_all_taylor[0],'--ok',lw=2,label="$w_{0m}$/$L_0$=0.003",markevery=20,ms=6)
plt.plot(wm_all_taylor1[0],xiegamaM_all_taylor1[0],'-*k',lw=2,label="$w_{0m}$/$L_0$=0.003-AxialVar",markevery=10,ms=6)
plt.plot(wm_all_taylor[1],xiegamaM_all_taylor[1],'--Dr',lw=2,label="$w_{0m}$/$L_0$=0.007",markevery=20,ms=6)
plt.plot(wm_all_taylor1[1],xiegamaM_all_taylor1[1],'-^r',lw=2,label="$w_{0m}$/$L_0$=0.007-AxialVar",markevery=10,ms=6)
plt.plot(wm_all_taylor[2],xiegamaM_all_taylor[2],'--8b',lw=2,label="$w_{0m}$/$L_0$=0.01",markevery=20,ms=6)
plt.plot(wm_all_taylor1[2],xiegamaM_all_taylor1[2],'-vb',lw=2,label="$w_{0m}$/$L_0$=0.01-AxialVar",markevery=10,ms=6)
ax=plt.gca()
ax.set_xlabel(xlabel="Buckling amplitude $w_m$ (mm)",fontsize=15)
ax.set_ylabel(ylabel="Max.compression stress $\sigma_m$ (MPa)",fontsize=15)
plt.xlim(0,7000)
plt.ylim(0,800)
plt.legend(loc='best',prop = {'size':10})


fig8,axe8=plt.subplots()
plt.plot(wm_all_taylor[0],P_all_taylor[0]/1e6,'--ok',lw=2,label="$w_{0m}$/$L_0$=0.003",markevery=20,ms=6)
plt.plot(wm_all_taylor1[0],P_all_taylor1[0]/1e6,'-*k',lw=2,label="$w_{0m}$/$L_0$=0.003-AxialVar",markevery=10,ms=6)
plt.plot(wm_all_taylor[1],P_all_taylor[1]/1e6,'--Dr',lw=2,label="$w_{0m}$/$L_0$=0.007",markevery=20,ms=6)
plt.plot(wm_all_taylor1[1],P_all_taylor1[1]/1e6,'-^r',lw=2,label="$w_{0m}$/$L_0$=0.007-AxialVar",markevery=10,ms=6)
plt.plot(wm_all_taylor[2],P_all_taylor[2]/1e6,'--8b',lw=2,label="$w_{0m}$/$L_0$=0.01",markevery=20,ms=6)
plt.plot(wm_all_taylor1[2],P_all_taylor1[2]/1e6,'-vb',lw=2,label="$w_{0m}$/$L_0$=0.01-AxialVar",markevery=10,ms=6)
ax=plt.gca()
ax.set_xlabel(xlabel="Buckling amplitude $w_m$ (mm)",fontsize=15)
ax.set_ylabel(ylabel="Force within a buckle P (MN)",fontsize=15)
plt.xlim(0,7000)
plt.ylim(0,8)
plt.legend(loc='best',prop = {'size':10})


fig9,axe9=plt.subplots()
plt.plot(wm_all_taylor[0],deltaT_all_taylor[0],'--ok',linewidth=2,label="$w_{0m}$/$L_0$=0.003",markevery=20,ms=6)
plt.plot(wm_all_taylor1[0],deltaT_all_taylor1[0],'-*k',linewidth=2,label="$w_{0m}$/$L_0$=0.003-AxialVar",markevery=20,ms=6)
plt.plot(wm_all_taylor[1],deltaT_all_taylor[1],'--Dr',linewidth=2,label="$w_{0m}$/$L_0$=0.007",markevery=20,ms=6)
plt.plot(wm_all_taylor1[1],deltaT_all_taylor1[1],'-^r',linewidth=2,label="$w_{0m}$/$L_0$=0.007-AxialVar",markevery=20,ms=6)
plt.plot(wm_all_taylor[2],deltaT_all_taylor[2],'--8b',linewidth=2,label="$w_{0m}$/$L_0$=0.01",markevery=20,ms=6)
plt.plot(wm_all_taylor1[2],deltaT_all_taylor1[2],'-vb',linewidth=2,label="$w_{0m}$/$L_0$=0.01-AxialVar",markevery=20,ms=6)

ax=plt.gca()
ax.set_xlabel(xlabel="$w_m (mm)$",fontsize=15)
ax.set_ylabel(ylabel="$\Delta T$ ($^oC$)",fontsize=15)
plt.xlim(0,7000)
plt.ylim(0,140)
plt.legend(loc='best',prop = {'size':10})
plt.show()





