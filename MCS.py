
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


######################################################################################
# Using Mento Carlo method to comupte the Casino problem. roll dice !
#######################################################################################
def rollDice():
    roll=random.randint(1,100)
    if roll == 100:
        # print ("roll={}. roll was 100, you lose. (Play again!)".format(roll))
        return False
    elif roll<=50:
        # print("roll={}.roll was 1-50, you lose.".format(roll))
        return False
    elif roll>50 and roll<100:
        # print ("roll={}. roll was 51-99, you win! (play more!)".format(roll))
        return True

#  create a single bettor.
def bettor(fund,initial_wager,wager_count):
    '''
    :param fund: the money you have
    :param initial_wager: the initial bet value
    :param wager_count:  the number of bettor
    :return:
    '''
    value=fund
    wager=initial_wager
    Xwage=[]
    Ywage=[]
    count=1
    while count<=wager_count:
        if rollDice():
            value+=wager
            Xwage.append(count)
            Ywage.append(value)
        else:
            value-=wager
            Xwage.append(count)
            Ywage.append(value)
        count=count+1
        # print ("fund={}".format(value))
    plt.plot(Xwage,Ywage)

# num=0
# plt.figure(figsize=(8,6))
# while num<=100:
#     bettor(10000,100,100000)
#     num+=1
#
# ax=plt.gca()
# ax.set_xlabel("the number of betting time for Casino")
# ax.set_ylabel("the variation of left money after Casino after each betting")
# plt.show()




########################################################################################
# using Mento Carlo method to compute the area of irregular shape. for instance, y=x**2#
########################################################################################
def dotDist(a_x,b_x,a_y,b_y):
    '''
    a_x,b_x,a_y,b_y should follow the specific figure of function.
    '''
    dot_x=random.uniform(a_x,b_x)
    dot_y=random.uniform(a_y,b_y)
    return (dot_x,dot_y)

def fun(x_value,y_value):
    # y=x_value**3
    e=2.71828
    y=(e**(-1*x_value))/(1+(x_value-1)**2)
    sign=True if y_value<=y else False
    return sign

def comp_area_MCS(a_x,b_x,a_y,b_y,num):
    '''
    :param a: the low-boundary of function
    :param b: the upper-boundary of function
    :param num:  the MC simulation number.
    :return:  the area of function between a and b
    '''
    k=0
    for i in range (num):
        dot_x,dot_y=dotDist(a_x,b_x,a_y,b_y)
        sign=fun(dot_x, dot_y)
        # plt.scatter(dot_x,dot_y)
        if sign:
            k+=1
    area=(b_x-a_x)*(b_y-a_y)*k/num
    print ("the area of function between {} and {} is {}".format(a_x,b_x,area))

########################################################################################
# using Mento Carlo method 2 to compute the area of irregular shape. for instance, y=x**2#
########################################################################################

def dotDistribution(a,b):
    return random.uniform(a,b)
def fun(x):
    e=2.71828
    return (e**(-1*x))/(1+(x-1)**2)
def MCS_2(a,b,count):
    all_y=0
    for i in range(count):
        value_x=dotDistribution(a,b)
        all_y +=fun(value_x)
    value_ave=all_y/count
    area=(b-a)*value_ave
    return area

def get_crude_MC_variance(a,b,num_samples):
    """
    This function returns the variance fo the Crude Monte Carlo.
    Note that the inputed number of samples does not neccissarily
    need to correspond to number of samples used in the Monte
    Carlo Simulation.
    Args:
    - num_samples (int)
    Return:
    - Variance for Crude Monte Carlo approximation of f(x) (float)
    """
    # get the average of squares
    ave_all = 0
    squre_all=0
    for i in range(num_samples):
        x = dotDistribution(a,b)
        squre_all += fun(x)**2
        ave_all += fun(x)
    sum_of_sqs = squre_all*(b-a) / num_samples
    sum_of_ave=(ave_all*(b-a)/num_samples)**2
    return sum_of_sqs - sum_of_ave

# print (MCS_2(0,5,1000000))
# print (get_crude_MC_variance(0,5,100000))

########################################################################################
# using Mento Carlo method  to compute the failure probability of limit state
########################################################################################

def paraDist():
    P=random.gauss(10,2)
    L=random.gauss(8,0.1)
    W=random.gauss(100e-6,2e-5)
    T=random.gauss(600e3,1e5)
    return (P,L,W,T)

def limitState(P,L,W,T):
    g=W*T-P*L/4
    return g

def MC_pf(count):
    num_all=count
    sum_I=0
    for i in range(count):
        P,L,W,T=paraDist()
        g=limitState(P,L,W,T)
        I=1 if g<0 else 0
        sum_I+=I
    pf=sum_I/num_all
    return pf

def MC_beta(count):
    num_all=count
    sum_ave=0
    sum_ave2=0
    for i in range(count):
        P,L,W,T=paraDist()
        g=limitState(P,L,W,T)
        sum_ave+=g
        sum_ave2+=g**2
    ave=sum_ave/num_all
    std=np.sqrt(sum_ave2/num_all-(sum_ave/num_all)**2)
    beta=ave/std
    return beta


###############
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def fun(a,b):
    return a**3+b**3-18
def fun_deriv_a(a):
    return 3*a**2
def fun_deriv_b(b):
    return 3*b**2
def gradient2_a(a,sx1):
    grad2=(fun_deriv_a(a)*sx1)**2
    return grad2
def gradient2_b(b,sx2):
    grad2=(fun_deriv_b(b)*sx2)**2
    return grad2

# from scipy.stats import lognorm
# stddev = 0.859455801705594
# mean = 0.418749176686875
# dist=lognorm(s=stddev,loc=0,scale=np.exp(mean))
#
#
# x=np.linspace(0,6,200)
# y=1/(np.sqrt(2*np.pi)*x*stddev)*np.exp(-0.5*((np.log(x)-mean)/stddev)**2)
# cdf=dist.cdf(1)
# r=dist.rvs(1000)
# print (cdf)
# print (r)
# plt.plot(x,dist.pdf(x))
# plt.plot(x,y)
# plt.plot(x,dist.cdf(x))
# plt.hist(r, density=True, histtype='stepfilled', alpha=0.2)
# plt.show()


