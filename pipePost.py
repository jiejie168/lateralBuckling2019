__author__ = 'Jie'
'''
This programming is used for the post-processing of .obd files. 20/08/2019
Specifically, it will automatically output the curves data in nodes: displacement-LPF.
tips:
(1): The nodeSets used for history outputs should be prescribed in modelling so that the odb file
     will contain the same nodeSets with inputs in captials such as CENTER, and N_TOP.
(2) The nodeSets file can contain any number of nodes.
(3) For repeatly running, the data under "Tools/XY data/.." should be manually delected in case duplicated
    curves names. Error would be activated.

'''
import numpy as np
#import pandas as pd
#from pandas import Series, DataFrame
from abaqus import *
from abaqusConstants import *
import odbAccess
import testUtils
testUtils.setBackwardCompatibility()
import regionToolset
import displayGroupMdbToolset as dgm
from material_libLateralBuckling import material_lib
from os import environ      # the aim of these import is try to connect pycharm with abaqus liberaries.
environ['ABAQUS_BAT_PATH'] = 'C:\\SIMULIA\\Abaqus\\Commands\\abaqus' # allowing pycharm IDE check types and run abaqus python program
environ['ABAQUS_BAT_SETTING'] = 'noGUI'

class PipeDataProcess():
    '''
    '''
    def dataProcess(self):
        db1=(("file directory",'D:/temp/'),("odb name:",'job1-Riks.odb'),("1st nodeSetName(in Capital):",'N_TOP'),("2ed nodeSet:",'CENTER'))
        flag=True
        while flag:
            try:
                dir,odbName,N_TOP,CENTER=list(map(lambda x: str(x),getInputs(fields=db1,label="please input: file directory,.odb result file (including suffix), 1st nodeSet, 2nd nodeSet")))
            except (ValueError,TypeError):
                raise "Please input string type for odbName"
            if bool(odbName)==False or bool(dir)==False or bool(N_TOP)==False or bool(CENTER)==False:
                flag=True
            else:
                break

        #odb=session.openOdb(name=dir+ odbName)
        odb=odbAccess.openOdb(path=dir+ odbName,readOnly=False)
        print ("\nall the part in odb file are:\n")
        print (odb.parts)
        print ("\nall the steps in odb file are:\n")
        print (odb.steps)
        print ("\nall the instances in odb file are:\n")
        print (odb.rootAssembly.instances)
        #print (odb.__members__)
        #print (odb.rootAssembly.__memebers__)
        print ("\nall the nodeSets in odb file are:\n")
        print (odb.rootAssembly.nodeSets)
        print ("\nall the elementSets in odb file are\n")
        print (odb.rootAssembly.elementSets)

        ##### create the lpf-arc curve for the whole model.
        steps=[]
        for key, _ in odb.steps.items():
            steps.append(key)

        xy_result=session.XYDataFromHistory(name='lpf_arc',odb=odb,
                                            outputVariableName='Load proportionality factor: LPF for Whole Model',
                                            steps=(str(steps[-1]),),)
        # c1=session.Curve(xyData=xy_result)
        # xyp=session.xyPlots['XYPlot-1']
        # chartName=xyp.charts.keys()[0]
        # chart=xyp.charts[chartName]
        # chart.setValues(curvesToPlot=(c1,),)


        ###### create the dis-arc curve for selected nodes. here, the node in pipe center and
        ###### pipe ends are selected.
        # the former scripts are only used for the extraction of node label in a set, so that we can use it
        # for history curve outputs.
        nodeLabel_center=[]
        nodeLabel_top=[]
        step_f=odb.steps[steps[-1]]
        frame=step_f.frames[-1]
        centerSet=odb.rootAssembly.nodeSets[CENTER]  #for nodeSets,there is no need to explore in specific instances.
        topSet=odb.rootAssembly.nodeSets[N_TOP]
        disp_center=frame.fieldOutputs['U']
        disp_top=frame.fieldOutputs['U']
        centerDisp=disp_center.getSubset(region=centerSet)
        topDisp=disp_top.getSubset(region=topSet)
        node_num=0
        node1_num=0
        for v in centerDisp.values:
            #print (v)
            label=v.nodeLabel
            node_num=node_num+1
            nodeLabel_center.append(label)
            xy_result1 = session.XYDataFromHistory(name='dis-arc_center'+str(node_num), odb=odb,
                outputVariableName='Spatial displacement: U2 at Node ' + str(label) + ' in NSET '+ CENTER,
                steps=(str(steps[-1]),),)
            xy1=session.xyDataObjects['dis-arc_center'+str(node_num)]
            xy2 = session.xyDataObjects['lpf_arc']
            xy3 = combine(xy1, xy2)
            xy3.setValues(sourceDescription='combine (-"dis-arc_center"+str(node_num) , "lpf_arc" )')
            tmpName = xy3.name
            session.xyDataObjects.changeKey(tmpName, 'dis-lpf'+CENTER+str(label))

        for v1 in topDisp.values:
            label=v1.nodeLabel
            node1_num=node1_num+1
            nodeLabel_top.append(label)
            xy_result1 = session.XYDataFromHistory(name='dis-arc_top'+str(node_num), odb=odb,
                outputVariableName='Spatial displacement: U1 at Node ' + str(label) + ' in NSET '+ N_TOP,
                steps=(str(steps[-1]),),)
            xy1=session.xyDataObjects['dis-arc_top'+str(node_num)]
            xy2 = session.xyDataObjects['lpf_arc']
            xy3 = combine(-xy1, xy2)
            xy3.setValues(sourceDescription='combine (-"dis-arc_top"+str(node_num) , "lpf_arc" )')
            tmpName = xy3.name
            session.xyDataObjects.changeKey(tmpName, 'dis-lpf'+N_TOP+str(label))
        print("history outputs in the nodeSets are completed.")

        odb.save()
        #odb.close()


def main():
    pipeData=PipeDataProcess()
    pipeData.dataProcess()

if __name__=='__main__':
    main()