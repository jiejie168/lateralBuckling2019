__author__ = 'Jie'
'''
This script is deployed to create a plug-in app for the production of a pipe part.
The module pipeCreate.py should be associated, consisting of a function of pipe named "createPerfectPipe".
07-05-2019
'''
#import numpy as np
#from matplotlib import pyplot as plt
from abaqusGui import *
from abaqusConstants import *
import testUtils
testUtils.setBackwardCompatibility()
import regionToolset
from os import environ
environ['ABAQUS_BAT_PATH'] = 'C:\\SIMULIA\\Abaqus\\Commands\\abaqus'
environ['ABAQUS_BAT_SETTING'] = 'noGUI'

# main GUI , write a form mode, containing only one dialog box
class pipeModel(AFXForm):
    def __init__(self,owner):
        AFXForm.__init__(self,owner)
        #the object "pipeCreate" should be the name of .py file; the method "createPerfectPipe" is one of the fuction in  pipeCreate.py
        self.cmd=AFXGuiCommand(mode=self,method='createPerfectPipe',objectName="pipeCreate")
        self.namekw=AFXStringKeyword(self.cmd,"name",TRUE)
        self.odkw=AFXFloatKeyword(self.cmd,"od",TRUE)
        self.tkw=AFXFloatKeyword(self.cmd,"t",TRUE)
        self.Lkw=AFXFloatKeyword(self.cmd,"L",TRUE)

    def getFirstDialog(self):
        self.cmd.setKeywordValuesToDefaults
        return pipeModelDB(self)

## Form to launch , it is a dialog constructor
class pipeModelDB(AFXDataDialog):
    def __init__(self,form):
        AFXDataDialog.__init__(self,form,"Pipe Part Creator",\
                               self.OK|self.CANCEL,DIALOG_ACTIONS_SEPARATOR)
        v=FXHorizontalFrame(self)
        va=AFXVerticalAligner(v)
        #pic = FXXPMIcon(getAFXApp(), IconData)
        #FXLabel(v, text='', ic=pic, opts=LAYOUT_CENTER_Y)
        AFXTextField(va,15,'Part Name:',form.namekw,0)
        AFXTextField(va,15,'Outer Diameter(od):',form.odkw,0)
        AFXTextField(va,15,'Pipe Thickness(t):',form.tkw,0)
        AFXTextField(va,15,'Pipe Length(L):',form.Lkw,0)


#/# #### Register toolset ###############################################
from symbolicConstants import SymbolicConstant

toolset=getAFXApp().getAFXMainWindow().getPluginToolset()
#A String specifying the string sent to the kernel the first time this plug-in is invoked
# here is by import the pipeCreate module
toolset.registerGuiMenuButton(buttonText='Geometry|pipeCreate',object=pipeModel(toolset),kernelInitString='import pipeCreate',\
                             version='1.0',author='Jie Cai',description='Creates a pipe part in the current model' )



