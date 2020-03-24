__author__ = 'Jie'

import numpy as np
import matplotlib.pyplot  as plt
import math
import os, sys
import numpy.linalg as lg


class CSproperty:
    '''
    this class is used to calculate all the cross-section properties of structures, including pipe cross-section, circular,
    rectangular, etc.
    '''
    def __init__(self,judge="True"):
        self.judge=judge

    def PipeCs(self,D,d):

        J= math.pi*D**4*(1-(d/D)**4)/32    # the torsion moment of inertia
        Ix=Iy=J/2                          # the bending moment of inertia
        print ("\nthe torsion moment of inertia (J) of this pipeline  is:{:.2f} mm^4\n".format(J))
        print("\nthe bending moment of inertia (Ix=Iy) of this pipeline is:{:.2f}mm^4\n".format(Ix))
        return J, Ix, Iy

def main():
    print ("\nplease input the type of a cross-section (such as: pipe,circular,rectangular,and I ):\n")
    section_Type=input()
    types=section_Type.strip()

    if types== "pipe":
        print("\n please input the dimension of the pipe cross-section (outer diameter and inner diameter):\n")
        inputs=sys.stdin.read()   # using ctrl+ D t terminate the input. using enter to confirm each parameter
        D,d,*others=list(map(float,inputs.split()))
    elif types=="circular":
        print("\n please input the dimension of the circular cross-section (diameter):\n")
        inputs=sys.stdin.read()
        D,*others=list(map(float,inputs.split()))
    elif types=="rectangular":
        print ("\n please input the dimension of the rectangular cross-section (width and height):\n")
        inputs=sys.stdin.read()
        w,h,*others=list(map(float,inputs.split()))
    else:
        print ("\n please input the dimension of the I cross-section (WEB t and height, UPER t and width, Bottom t and width):\n")
        inputs=sys.stdin.read()
        t_w,h_w,t_u,h_u,t_b,h_b,*others=list(map(float,inputs.split()))


    CSPro=CSproperty()
    CSPro.PipeCs(D,d)

if __name__ == '__main__':
        main()