__author__ = 'Jie'

import numpy as np
from abaqus import *
from abaqusConstants import *
from os import environ
environ['ABAQUS_BAT_PATH'] = 'C:\\SIMULIA\\Abaqus\\Commands\\abaqus'
environ['ABAQUS_BAT_SETTING'] = 'noGUI'

def createPerfectPipe(name,od,t,L):
    '''
    create a part for a pipe without initial imperfection
    The part will used for the following calcualtion.
    name: the name of pipe part
    od: the outer diameter of pipe
    t: the thickness of pipe
    L: the length of pipe
    :return:
    '''
    import part
    vp=session.currentViewportName
    m=mdb.models['Model-1']
    s=m.Sketch(name='__profile__',sheetSize=L*1.5)
    s.Line(point1=(0,0),point2=(L,0))
    p=m.Part(name='pipe-part',dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
    p.BaseWire(sketch=s)

def main():
    createPerfectPipe('pipe1',214,14.3,420)

if __name__ == '__main__':
    main()

