__author__ = 'Jie'
# coding=utf8
'''
Calculate the pipe buckling properties based on five different Hobbs’ buckling modes.
The input data file is : “pipeParameters.txt”.
Six different corresponding imperfection shapes of pipes can be output via this code.
Critical buckling force based on elastic foundation method according
to Timoshenko and Gere, 1961 can be calculated.
'''
import numpy as np
from matplotlib import pyplot as plt
import sys
from scipy.optimize import fsolve,root
#sys.path.append('D:\temp\postDoc-Brest\lateralBuckling2019')
from PhilippeEquations import PipeTorsionPhilippe

class PipeData():
    '''
    #a,b,c=[float(el)  for el in input(" input 3 parameters:").split()]
    this script is used for obtaining pipe data for calculations based on a file: "pipeParameters.txt"
    the sequence of data is :  diameter, thickness, length, unit weight, lateral friction...
    '''
    def __init__(self,modeName="bucklingMode1"):
        self.modeName=modeName

    def kValue(self):
        if self.modeName=="bucklingMode1":
            k1,k2,k3,k4,k5=[80.76,6.391e-5,0.5,2.407e-3,0.06938]
        elif self.modeName=="bucklingMode2":
            k1,k2,k3,k4,k5=[4*np.pi**2,1.743e-4,1.0,5.532e-3,0.1088]
        elif self.modeName=="bucklingMode3":
            k1,k2,k3,k4,k5=[34.06,1.668e-4,1.294,1.032e-2,0.1434]
        else:
            k1,k2,k3,k4,k5=[28.20,2.144e-4,1.608,1.047e-2,0.1483]
        return k1,k2,k3,k4,k5

    def pipeData(self):
        with open("pipeParameters.txt",'r') as pipeFile:
            data=[]
            for eachline in pipeFile.readlines():
                ele=eachline.split(':')
                data.append(ele[1].strip())
            od,t,L,q,f_L,f_A,L_b,E,alpha,deltaT,deltaP=list(map(lambda x: float(x),data))
            #pipeFile.close()
        return od,t,L,q,f_L,f_A,L_b,E,alpha,deltaT,deltaP

    def imperfData(self):
        with open("imperfData.txt","r") as imperfFile:
            data=[]
            for eachline in imperfFile.readlines():
                ele=eachline.split(':')
                data.append(ele[1].strip())
            w0m,L0=list(map(lambda x: float(x),data))
            #imperfFile.close()
        return w0m,L0

class PipeHobbsModes():
    '''
    This class is used to calculate the five different Hobbs,1984 modes of pipe buckling
    P: the reduced axial force within a buckle
    P0: the critical buckling force of a perfectly straight pipeline laid on a rigid flat seabed
    w_max:the maximum amplitude of a buckle
    M_max: the maximum bending moment
    k1..k5: coefficients of each analytical solutions.
    f_L: frictional coefficient of the seabed
    L_b:the half-wave length of buckling of pipe.The definition is based on Hobbs' definition.
    I:bending area moment of inertia (2nd moment of area)
    S:steel area of pipe cross-section
    E=2.06e5      # units system is always MPa, mm, N as a convention
    '''
    def __init__(self,E,f_L,q,L,L_b,I,A,k1,k2,k3,k4,k5):
        self.k1,self.k2,self.k3,self.k4,self.k5=k1,k2,k3,k4,k5
        self.f_L=f_L
        self.q=q
        self.L=L
        self.L_b=L_b
        self.I=I
        self.E=E
        self.P=self.k1*self.E*self.I/self.L_b**2
        self.A=A
        s_e=A*self.E*self.f_L*self.q*self.L_b**5/(self.E*self.I)**2
        self.P0=self.P+self.k3*self.f_L*self.q*self.L_b*((1.0+self.k2*s_e)**0.5-1.0)
        self.w_max=self.k4*self.f_L*self.q*self.L_b**4/self.E/self.I
        self.M_max=self.k5*self.f_L*self.q*self.L_b**2

    def buckleResults(self):
        return self.P,self.P0,self.w_max,self.M_max

    def bucklingMode5(self):
        '''
        this mode is only used for a displaying!
        :return:
        '''
        # m=self.f_L*self.q/self.E/self.I
        # n=np.sqrt(self.P/self.E*self.I)
        x5=np.linspace(-1.5*np.pi,1.5*np.pi,num=100)
        # w5=(m/n**4)*(1+n**2*self.L_b**2/8-n**2*x5**2/2-np.cos(n*x5)/np.cos(n*self.L_b/2))
        w5=np.cos(x5)

        plt.plot(x5,w5,'k-',linewidth=4.0,label="Mode 5 buckling")
        plt.rc('legend',fontsize=13)
        plt.xticks([])
        plt.yticks([])
        ax=plt.gca()
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        plt.legend(frameon=False,loc='best')

        plt.show()

    def bucklingMode1(self):
        '''
        the original deflection function of Mode 1 from hobbs is actually the infinite
        mode. Here, the shape function is from Kerr,1979
        In Karampour, 2013:there are still some Mode 1 shapes with different curvature ratio.
        '''
        w_amplitude=(self.f_L*self.q*self.L_b**4)/(self.E*self.I)
        n=4.493*2/self.L_b
        x=np.linspace(-self.L_b/2,self.L_b/2,num=100)
        x1=np.linspace(self.L_b/2,self.L/2,num=100)  # the pipe segment outside of buckle
        x2=np.linspace(-self.L/2,-self.L_b/2,num=100) # the pipe segment outside of buckle
        y1=np.zeros((100,))
        w1=w_amplitude/32/(4.493)**2*(1-(2*x/self.L_b)**2-2*(np.cos(n*x)-np.cos(4.493))/(4.493*np.sin(4.493)))
        plt.plot(x,w1,'k-',linewidth=4,label='Mode 1 buckling')
        plt.plot(x1,y1,'k-',linewidth=8)
        plt.plot(x2,y1,'k-',linewidth=8)
        plt.rc('legend',fontsize=13)
        plt.xlim(-self.L/2,self.L/2)
        plt.ylim(0,)
        #plt.xticks([])
        #plt.yticks([])
        ax=plt.gca()
        ax.set_xlabel(xlabel="Pipe length (mm)",fontsize=13)
        ax.set_ylabel(ylabel="Deflection (mm)", fontsize=13)

        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',-5))
        ax.spines['left'].set_position(('data',0))
        plt.legend(frameon=False)
        #plt.ion()
        #plt.show()
        #plt.pause(2)

    def bucklingMode2(self):
        x1=np.linspace(0,self.L_b,num=100)
        x2=np.linspace(-self.L_b,0,num=100)
        w_amplitude=(self.f_L*self.q*self.L_b**4)/(16*np.pi**4*self.E*self.I)
        w21=-w_amplitude *(1-np.cos(2*np.pi*x1/self.L_b)+\
                                        np.pi*np.sin(2*np.pi*x1/self.L_b)+(2*np.pi**2*x1/self.L_b)*(1-x1/self.L_b))
        w22=w_amplitude *(1-np.cos(2*np.pi*x2/self.L_b)-\
                                        np.pi*np.sin(2*np.pi*x2/self.L_b)-(2*np.pi**2*x2/self.L_b)*(1+x2/self.L_b))

        plt.plot(x1,w21,'k-',label="Mode 2 buckling",linewidth=4)
        plt.plot(x2,w22,'k-',linewidth=4)
        plt.rc('legend',fontsize=13)

        x3=np.linspace(self.L_b,self.L/2,num=100)  # the pipe segment outside of buckle
        x4=np.linspace(-self.L/2,-self.L_b,num=100) # the pipe segment outside of buckle
        y3=np.zeros((100,))
        plt.plot(x3,y3,'k-',linewidth=4)
        plt.plot(x4,y3,'k-',linewidth=4)

        plt.xlim(-self.L/2,self.L/2)
        #plt.ylim(,0)
        # plt.xticks([])
        # plt.yticks([])
        ax=plt.gca()
        ax.set_xlabel(xlabel="Pipe length (mm)",fontsize=13)
        ax.set_ylabel(ylabel="Deflection (mm)",fontsize=13)
        ax.xaxis.set_label_coords(0.25, 0.43)
        ax.yaxis.set_label_coords(0.55,0.87)

        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        plt.legend(frameon=False)
        #plt.ion()
        #plt.show()
        #plt.pause(2)


    def bucklingMode3(self):
        '''
        The expression of equation is from Hong,2015.
        Lb is the buckling length of Mode3 based on Hobbs,which is the centeral lobe width of the buckle.
        tips: For correct display, the entire pipe length L should be larger than 2.59L_b according to definition.

        '''
        w_amplitude=(self.f_L*self.q*self.L_b**4)/(self.E*self.I)
        n3=2.918*2/self.L_b
        x1=np.linspace(-self.L_b/2,self.L_b/2,num=100)
        w31=w_amplitude*(2.59/2)**4*(1.2e-3*np.cos(n3*x1)-(x1/self.L_b)**2*2/(2.59*7.551)**2+2.48e-3)

        x2=np.linspace(self.L_b/2,2.59*self.L_b/2,num=100)
        w32=w_amplitude*(2.59/2)**4*(0.6e-3*np.cos(n3*x2)+0.14e-3*np.sin(n3*x2)-1.36e-2*(x2/self.L_b)*(2/2.59)+\
                         4.48e-3+(x2/self.L_b)**2*(2/(7.551*2.59)**2))

        x3=np.linspace(-self.L_b/2,-2.59*self.L_b/2,num=100)
        w33=w_amplitude*(2.59/2)**4*(0.6e-3*np.cos(n3*x3)-0.14e-3*np.sin(n3*x3)+1.36e-2*(x3/self.L_b)*(2/2.59)+\
                         4.48e-3+(x3/self.L_b)**2*(2/(7.551*2.59)**2))

        plt.plot(x1,w31,'k-',label="Mode 3 buckling",linewidth=4.0)
        plt.plot(x2,w32,'k-',linewidth=4.0)
        plt.plot(x3,w33,'k-',linewidth=4.0)
        plt.rc('legend',fontsize=13)


        x4=np.linspace(2.59*self.L_b/2,self.L/2,num=100)  # the pipe segment outside of buckle
        x5=np.linspace(-self.L/2,-2.59*self.L_b/2,num=100) # the pipe segment outside of buckle
        y4=np.zeros((100,))
        plt.plot(x4,y4,'k-',linewidth=4)
        plt.plot(x5,y4,'k-',linewidth=4)


        #plt.xlim(-2.59*self.L_b/2,2.59*self.L_b/2)
        #plt.ylim(0,)
        # plt.xticks([])
        # plt.yticks(([]))
        ax=plt.gca()
        ax.set_xlabel(xlabel="Pipe length (mm)",fontsize=13)
        ax.set_ylabel(ylabel="Deflection (mm)",fontsize=13)
        ax.xaxis.set_label_coords(0.85,0.33)
        ax.yaxis.set_label_coords(0.45,0.5)

        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        #plt.ion()
        plt.legend(frameon=False)
        plt.show()


    def bucklingMode4(self):
        '''
        the deflection equation comes from Kerr,1979
        '''
        n=5.31/self.L_b
        l=1.61*self.L_b
        w_amplitude=(self.f_L*self.q*l**4)/(self.E*self.I)
        A1,A2,A3,A4=-0.19e-3*w_amplitude,0.57e-3*w_amplitude,4.89e-3*w_amplitude/l,0.19e-3*w_amplitude
        A5,A6,A7,A8=0.02e-3*w_amplitude,0.26e-3*w_amplitude,-12.16e-3*w_amplitude/l,5.12e-3*w_amplitude

        x1=np.linspace(0,self.L_b)
        x2=np.linspace(self.L_b,l)
        w41=A1*np.cos(n*x1)+A2*np.sin(n*x1)+A3*x1+A4-w_amplitude*x1**2/(2*n**2)/l**4
        w42=A5*np.cos(n*x2)+A6*np.sin(n*x2)+A7*x2+A8+w_amplitude*x2**2/(2*n**2)/l**4

        x3=np.linspace(-self.L_b,0)
        x4=np.linspace(-l,-self.L_b)
        w43=-A1*np.cos(n*x3)+A2*np.sin(n*x3)+A3*x3-A4+w_amplitude*x3**2/(2*n**2)/l**4
        w44=-A5*np.cos(n*x4)+A6*np.sin(n*x4)+A7*x4-A8-w_amplitude*x4**2/(2*n**2)/l**4

        plt.plot(x1,w41,'k-',label='Mode 4 buckling',linewidth=4.0)
        plt.plot(x2,w42,'k-',linewidth=4.0)
        plt.plot(x3,w43,'k-',linewidth=4.0)
        plt.plot(x4,w44,'k-',linewidth=4.0)

        x5=np.linspace(l,self.L/2,num=100)  # the pipe segment outside of buckle
        x6=np.linspace(-self.L/2,-l,num=100) # the pipe segment outside of buckle
        y5=np.zeros((100,))
        plt.plot(x5,y5,'k-',linewidth=4)
        plt.plot(x6,y5,'k-',linewidth=4)

        plt.rc('legend',fontsize=13)
        # plt.xticks([])
        # plt.yticks([])
        ax=plt.gca()
        ax.set_xlabel(xlabel="Pipe length (mm)",fontsize=13)
        ax.set_ylabel(ylabel="Deflection (mm)",fontsize=13)
        ax.xaxis.set_label_coords(0.15,0.46)
        ax.yaxis.set_label_coords(0.48,0.80)

        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        plt.legend(frameon=False)
        plt.show()

class BOEF():
    '''
    This is the class of beam-on -elastic foundation method based on Timoshenko and Gere, 1961.
    For a series of given buckling length L_b supported on an elastic foundation. the criticla buckling
    force is calculated.
    the results have been comparable to McCarron,2019, but there is a little difference.
    '''
    def __init__(self,L_b,n,ks,E,I):
        self.L_b=L_b
        self.n=n
        self.ks=ks
        self.E=E
        self.I=I
    def pipe_BOEF(self):
        temp=self.n**2+self.ks*self.L_b**4/self.n**2/np.pi**4/self.E/self.I
        P0=np.pi**2*self.E*self.I*(temp)/self.L_b**2
        return P0

class PipeImp(PipeTorsionPhilippe):
    '''
    This script defines all the five different imperfections modes based on Taylor,1986. The original buckling modes
    are from Hobbs,1984  research.
    '''

    def __init__(self,D,t,L_b,k,E,alpha_t,A,J,I,w0m,L0):
        self.w0m=w0m
        self.L0=L0
        PipeTorsionPhilippe.__init__(self,D,t,L_b,k,E,alpha_t,A,J,I)

    def bucklingMode1_shapes(self):
        '''
        this function is only used for the visulization of
        different imperfections in terms of Mode 1 with different shapes from Zhang Xinhu, 2018.
        '''
        x1=np.linspace(-self.L0/2,0,num=50)
        w011=self.w0m*(8/3 * (2*x1/self.L0)**2-3*2*x1/self.L0+1)*(1+2*x1/self.L0)**3

        x2=np.linspace(0,self.L0/2,num=50)
        w012=self.w0m*(8/3 * (2*x2/self.L0)**2+3*2*x2/self.L0+1)*(1-2*x2/self.L0)**3

        x3=np.linspace(-self.L0/2,self.L0/2,num=50)
        w02=self.w0m/2*(1+np.cos(2*np.pi*x3/self.L0))

        # w031=-self.w0m*(4*2*x1/self.L0-1)*(2*x1/self.L0+1)
        # w032=self.w0m*(4*2*x2/self.L0+1)*(2*x2/self.L0-1)

        w031=self.w0m*(1+2.5*(2*x1/self.L0)**3-1.5*(2*x1/self.L0)**5)
        w032=self.w0m*(1-2.5*(2*x2/self.L0)**3+1.5*(2*x2/self.L0)**5)

        w041=self.w0m*(1+4*(2*x1/self.L0)**3+3*(2*x1/self.L0)**4)
        w042=self.w0m*(1-4*(2*x2/self.L0)**3+3*(2*x2/self.L0)**4)

        plt.plot(x1,w011,'k-',label="N.1",linewidth=2.0)
        plt.plot(x2,w012,'k-',linewidth=2.0)
        plt.plot(x3,w02,'r-',label="N.2",linewidth=2.0)
        plt.plot(x1,w031,'b-',label="N.3",linewidth=2.0)
        plt.plot(x2,w032,'b-',linewidth=2.0)
        plt.plot(x1,w041,'g-',label="N.4",linewidth=2.0)
        plt.plot(x2,w042,'g-',linewidth=2.0)

        plt.xticks([])
        plt.yticks([])
        ax=plt.gca()
        ax.set_xlabel(xlabel="$L_0$",fontsize=20)
        ax.set_ylabel(ylabel="$W_0$",fontsize=20)

        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        #plt.ion()
        plt.legend()
        plt.show()

    def imperfMode2(self):


        x1=np.linspace(0,self.L0,num=100)
        x2=np.linspace(-self.L0,0,num=100)

        w21=-self.w0m/8.62 *(1-np.cos(2*np.pi*x1/self.L0)+\
                                        np.pi*np.sin(2*np.pi*x1/self.L0)+(2*np.pi**2*x1/self.L0)*(1-x1/self.L0))
        w22=self.w0m/8.62 *(1-np.cos(2*np.pi*x2/self.L0)-\
                                        np.pi*np.sin(2*np.pi*x2/self.L0)-(2*np.pi**2*x2/self.L0)*(1+x2/self.L0))

        plt.plot(x1,w21,'k-',linewidth=4,label="Mode 2 imperfection")
        plt.plot(x2,w22,'k-',linewidth=4)

        plt.rc('legend',fontsize=13)
        ax=plt.gca()
        ax.set_xlabel(xlabel="Lo",fontsize=13)
        ax.set_ylabel(ylabel="Wo",fontsize=13)
        ax.xaxis.set_label_coords(0.55,0.48)
        ax.yaxis.set_label_coords(0.48,0.80)
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))
        plt.legend()
        plt.show()

    def imperfMode3(self):
        '''
        Impperfection is based on Kerr, 1978. Mode 3 buckling. It is the simplified version fo bucklingMode3 in
        the former class.
        here, the L_0 is defined as the full length of entire imperfection.
        '''
        l=self.L0/2
        l1=2*l/2.59
        n3=7.551/l

        x1=np.linspace(-l1/2,l1/2,num=100)
        x2=np.linspace(l1/2,2.59*l1/2,num=100)
        x3=np.linspace(-l1/2,-2.59*l1/2,num=100)

        w31=self.w0m/1.484*(1+0.484*np.cos(n3*x1)-2.109*(x1/l1)**2)
        w32=self.w0m/0.821*(1+0.134*np.cos(5.836*x2/l1)+0.031*np.sin(5.836*x2/l1)-2.341*(x2/l1)+1.167*(x2/l1)**2)
        w33=self.w0m/0.821*(1+0.134*np.cos(5.836*x3/l1)-0.031*np.sin(5.836*x3/l1)+2.341*(x3/l1)+1.167*(x3/l1)**2)

        plt.plot(x1,w31,'k-',linewidth=4,label='Mode 3 imperfection')
        plt.plot(x2,w32,'k-',linewidth=4)
        plt.plot(x3,w33,'k-',linewidth=4)
        plt.rc('legend',fontsize=13)

        ax=plt.gca()
        ax.set_xlabel(xlabel="Lo",fontsize=13)
        ax.set_ylabel(ylabel="Wo",fontsize=13)
        ax.xaxis.set_label_coords(0.55,0.18)
        ax.yaxis.set_label_coords(0.48,0.75)
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))

        #plt.xlim(-L_b/2,L_b/2)
        #plt.ylim(0,)
        plt.legend()
        plt.show()

    def imperfMode4(self):
        '''
        Impperfection is based on Kerr, 1978. Mode 4 buckling.It is the simplified version fo bucklingMode4 in
        the former class.
        here, the L_0 is defined as the full length of entire imperfection.
        '''
        l=self.L0/2
        n4=8.54/l  # L_b is the full length of imperfection, Mode 4.
        x1=np.linspace(0,l/1.61,num=100)
        x2=np.linspace(l/1.61,l,num=100)
        x3=np.linspace(-l/1.61,0,num=100)
        x4=np.linspace(-l,-l/1.61,num=100)

        w41=self.w0m/8.2*(1-np.cos(8.54*x1/l)+3*np.sin(8.54*x1/l)+25.74*(x1/l)-36.11*(x1/l)**2)
        w42=self.w0m/0.304*(1+0.00391*np.cos(8.54*x2/l)+0.05078*np.sin(8.54*x2/l)-2.375*(x2/l)+1.3398*(x2/l)**2)
        w43=-self.w0m/8.2*(1-np.cos(8.54*x3/l)-3*np.sin(8.54*x3/l)-25.74*(x3/l)-36.11*(x3/l)**2)
        w44=-self.w0m/0.304*(1+0.00391*np.cos(8.54*x4/l)-0.05078*np.sin(8.54*x4/l)+2.375*(x4/l)+1.3398*(x4/l)**2)

        plt.plot(x1,w41,'k-',linewidth=4,label='Mode 4 imperfection')
        plt.plot(x2,w42,'k-',linewidth=4)
        plt.plot(x3,w43,'k-',linewidth=4)
        plt.plot(x4,w44,'k-',linewidth=4)
        plt.rc('legend',fontsize=13)

        ax=plt.gca()
        ax.set_xlabel(xlabel="Lo",fontsize=13)
        ax.set_ylabel(ylabel="Wo",fontsize=13)
        ax.xaxis.set_label_coords(0.65,0.46)
        ax.yaxis.set_label_coords(0.48,0.80)
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data',0))
        ax.spines['left'].set_position(('data',0))

        #plt.xlim(-L_b/2,L_b/2)
        #plt.ylim(0,)
        plt.legend()
        plt.show()

    def imperfMode6(self):
        Pcn,temp1,temp2=self.find_Pc_freeTorsion(0.00001)
        a=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))+self.xi-2*self.beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(self.gama/self.beta)*np.sqrt(2*np.sqrt(2*self.beta*self.xi*(2*self.alpha*self.beta-Pcn))-self.xi+2*self.beta*Pcn)
        x=np.linspace(-0.5,0.5,num=1000)
        w=1-temp1*np.cosh(a*x)*np.cos(b*x)+temp2*np.sinh(a*x)*np.sin(b*x)
        # plt.plot(x,w,'k-',linewidth=4,label="Imperfection mode with torsion")
        # plt.legend()
        # plt.show()
        return a,b

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

def main():
    pipeData=PipeData("bucklingMode3")
    od,t,L,q,f_L,f_A,L_b,E,alpha,deltaT,deltaP=pipeData.pipeData()
    k1,k2,k3,k4,k5=pipeData.kValue()
    w0m,L0=pipeData.imperfData()
    #od_1,t1,L1,L_b1,k,E,alphaT1=pipeData.philippeData()   #the specific case for Philippe calculation

    print ("##############################################################################################")
    print ("The parameters based on Hobbs'1984 for pipe buckling with different Mode (k1,k2,k3,k4,k5) is:\n",
           k1,k2,k3,k4,k5)
    print ("Inputs pipe data for calculation (od,t,L_max,L_b,q,f_L,f_A,E,alpha): \n",
        od,t,L,L_b,q,f_L,f_A,E,alpha)
    print ("There is no imperfection inputs needed for Hobbs's calculation.")
    print ("############################################################################################\n")

    pipeProp=PipeProp(L,L0,w0m,od,t,E,alpha,deltaT,deltaP)
    I,J,A,P0,L_expand=pipeProp.pipePro()
    print ("##############################################################################################")
    print ("the pipe moment of inertial is:\
     \n{:.2f} mm4,\nThe pipe steel area of cross-section is:\
     \n {:.2f} mm2,\nThe increased pipe length due to Temp increased {} &Pressure increased {} is:\
      \n{:.2f} mm,\nThe corresponding axial force P0 due to thermal$Pressure is:"
           " \n{:.2f}N\n".format(I,A,deltaT,deltaP,L_expand,P0))
    print ("############################################################################################\n")

    pipeHobbsModes=PipeHobbsModes(E,f_L,q,L,L_b,I,A,k1,k2,k3,k4,k5)
    pipeHobbsModes.bucklingMode1()
    P,P0,w_max,M_max=pipeHobbsModes.buckleResults()
    print ("##############################################################################################")
    print ("the corresponding pipe features without imperfection with the buckle length of L_b {}mm"
           " under this buckling mode is (based Hobbs):\
    \nP: {:.1f}N;\nP0:{:.1f}N;\nw_max:{:.5f}mm;\nM_max:{:.1f}Nmm\n ".format(L_b,P,P0,w_max,M_max))
    print ("############################################################################################\n")

    print ("############################################################################################\n")
    print ("Given a series of buckling length, plot the diagrams of deltaT,P0,and P with buckling amplitude via Hobbs")
    L_buckle=np.linspace(0.1,L,num=10000)

    P_all,P0_all,w_maxAll,M_maxAll=[],[],[],[]
    deltaT_all=[]
    for length in L_buckle:
        P_1=k1*E*I/length**2
        s_e=A*E*f_L*q*length**5/(E*I)**2
        P0_1=P_1+k3*f_L*q*length*((1.0+k2*s_e)**0.5-1.0)
        w_max1=k4*f_L*q*length**4/E/I
        M_max1=k5*f_L*q*L_b**2
        deltaT=P0_1/E/A/alpha
        P_all.append(P_1)
        P0_all.append(P0_1)
        w_maxAll.append(w_max1)
        M_maxAll.append(M_max1)
        deltaT_all.append(deltaT)
    fig1,ax1=plt.subplots()
    plt.plot(w_maxAll,P_all,'-k',linewidth=4,label="Perfect pipe_Hobbs")
    ax=plt.gca()
    ax.set_xlabel(xlabel="$w_m (mm)$",fontsize=15)
    ax.set_ylabel(ylabel="Force within a buckle P (N)",fontsize=15)
    plt.ylim(0,6e6)
    plt.xlim(0,6000)
    plt.legend()

    fig2,ax2=plt.subplots()
    plt.plot(w_maxAll,P0_all,'-k',linewidth=4,label="Perfect pipe_Hobbs")
    ax=plt.gca()
    ax.set_xlabel(xlabel="$w_m (mm)$",fontsize=15)
    ax.set_ylabel(ylabel="Buckling force P0 (N)",fontsize=15)
    plt.ylim(0,6e6)
    plt.xlim(0,6000)
    plt.legend()

    fig3,ax3=plt.subplots()
    plt.plot(w_maxAll,deltaT_all,'-k',linewidth=4,label="Perfect pipe_Hobbs")
    ax=plt.gca()
    ax.set_xlabel(xlabel="$w_m (mm)$",fontsize=15)
    ax.set_ylabel(ylabel="Temperature rise $\Delta T$ ($^oC$)",fontsize=15)
    plt.ylim(0,80)
    plt.xlim(0,6000)
    plt.legend()
    #plt.xlim(0,100)

    index=deltaT_all.index(min(deltaT_all))
    print ("the maximum temperature rise delatT is: {}oC".format(min(deltaT_all)))
    print ("the critical buckling force P0 is: {}KN".format(P0_all[index]/1e3))
    print ("the post buckling force P inside the buckle is: {}KN".format(P_all[index]/1e3))
    print ("the corresponding pipe buckled length is: {}m".format(L_buckle[index]/1e3))

    print ("############################################################################################\n")
    print ("Given a series of buckling length, plot the diagram of P0 with buckling amplitude via Timoshenko")

    #n=1
    P0_timoALL={}
    fig4,ax4=plt.subplots()
    for num,n in enumerate(np.arange(1,5)):
        for length1 in L_buckle:
            ks=f_L*q/0.1/od
            boef=BOEF(length1,n,ks,E,I)
            P0_timo=boef.pipe_BOEF()
            P0_timoALL.setdefault(num,[]).append(P0_timo)
        plt.plot(L_buckle,P0_timoALL[num],'-',linewidth=4,label="Perfect pipe_Timoshenko, n="+str(num+1))
        print ("the critical force of n={} via BOEF is: {} KN".format(n,min(P0_timoALL[num])/1e3))
    ax=plt.gca()
    ax.set_xlabel(xlabel="Buckle_length (mm)",fontsize=15)
    ax.set_ylabel(ylabel="Buckling capacity P0 (N)",fontsize=15)
    plt.ylim(0,5e6)
    plt.xlim(0,100000)
    plt.legend()
    plt.show()

    k=f_L/(0.1*od)
    pipeImp=PipeImp(od,t,L_b,k,E,alpha,A,J,I,w0m,L0)
    a,b=pipeImp.imperfMode6()
    print("a of Mode 6 is: {}".format(a))
    print ("b of Mode 6 is: {}".format(b))
    plt.show()

if __name__ == '__main__':
    main()
