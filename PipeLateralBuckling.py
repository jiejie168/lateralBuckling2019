__author__ = 'Jie'
#Filename:pipeLateralBuckling.py
'''
This script is for the overall calculation of lateral buckling of pipelines.
Different functions for pipes are developed.
'''
import numpy as np
from abaqus import *
from abaqusConstants import *
import testUtils
testUtils.setBackwardCompatibility()
import regionToolset
import displayGroupMdbToolset as dgm
from material_libLateralBuckling import material_lib
from os import environ      # the aim of these import is try to connect pycharm with abaqus liberaries.
environ['ABAQUS_BAT_PATH'] = 'C:\\SIMULIA\\Abaqus\\Commands\\abaqus' # allowing pycharm IDE check types and run abaqus python program
environ['ABAQUS_BAT_SETTING'] = 'noGUI'

class DataInput():
    '''
    This class is only used for the input data of simulation.
    default data come from Zeng,2014 (Critical upheaval buckling forces of imperfect pipelines).
    zengData: od=457, t=14.3, L=1200m,seabedWidth=40m, temp_initial=20oC,temp_final=120oC,q=1.5N/mm;u_L=0.6
    zeng: w0=200mm, L_b0=100m;
    Current caseData: McCarron, 2019
    '''
    def dataInput(self):
        db1=(('pipe outer Diameter(mm):','273'),('pipe thickness(mm):','18'),('pipe length(mm):','400000'),
             ('pipe part:','pipe1'),('seabed part:','seabed1'),('seabed width(mm):','40000'),
             ('pipe model:','model1'),('pipeMaterial:','Steel-e'))
        db2=(('pipe mesh size(mm):','1000'),('seabed mesh size(mm):','20000'),('initial temperature(oC):','20'),
             ('final temperature(oC):','120'),('job name:','job1-SD'),('flagEigen','0'))
        db3=(('w0(mm)','273'),('L_b(mm)','76000'))
        db4=(("unit downard force(N/mm):",'0.664'),("frictional coefficient(mu_L)",'0.95'))
        flag1=True
        while flag1:
            od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial=getInputs(fields=db1, label='Please enter the pipe dimensions')
            try:
                od,t,L,w_bed=float(od),float(t),float(L),float(w_bed)
            except (TypeError,ValueError):
                raise ValueError,'please input a float!!'
            pipePart,seabedPart,pipeModel,pipeMaterial=str(pipePart),str(seabedPart),str(pipeModel),str(pipeMaterial)
            if bool(od)==False or bool(t)==False or bool(L)==False or bool(w_bed)==False or\
                            bool(pipePart)==False or bool(seabedPart)==False or bool(pipeModel)==False:
                flag1=True
            else:
                break
        print ("input of pipe dimentsion completed")


        flag2=True
        while flag2:
            sPipe,sSeabed,initialT,finalT,jobName,flagEigen=getInputs(fields=db2,label='Please enter the mesh parameters')
            try:
                sPipe,sSeabed,initialT,finalT,flagEigen=float(sPipe),float(sSeabed),float(initialT),float(finalT),float(flagEigen)
            except (ValueError,TypeError):
                raise ValueError, 'please input a float!!'
            jobName=str(jobName)
            if bool(sPipe)==False or bool(sSeabed)==False or bool(initialT)==False or bool(finalT)==False or bool(jobName)==False :
                flag2=True
            else:
                break
        print ("input of mesh parameters and temperatures completed")


        flag3=True
        while flag3:
            w0,L_b=getInputs(fields=db3,label='Please enter the imperfection:')
            try:
                w0,L_b=float(w0),float(L_b)
            except (ValueError,TypeError):
                raise ValueError,'please input a float for imperfection value'
            if bool(w0)==False or bool(L_b)==False:
                flag3=True
            else:
                break
        print ("input of imperfection and other information completed")

        flag4=True
        while flag4:
            P_u,mu_L=getInputs(fields=db4,label="Please enter other information:")
            try:
                P_u,mu_L=float(P_u),float(mu_L)
            except (ValueError,TypeError):
                raise ValueError,'please input a float for load and friction'

            if bool(P_u)==False or bool(mu_L)==False:
                flag4=True
            else:
                break
        print("input of other information completed")

        return od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,sPipe,sSeabed,initialT,finalT,jobName,flagEigen,w0,L_b,P_u,mu_L

class LateralPipeBuckling(object):
    def __init__(self,od,t,L,pipePart,seabedPart,w_bed,pipeModel,m,sPipe,sSeabed,initialT,finalT,jobName,w0,L_b,P_u,mu_L):
        self.od=od
        self.t=t
        self.L=L
        self.pipePart=pipePart
        self.seabedPart=seabedPart
        self.w_bed=w_bed
        self.pipeModel=pipeModel
        self.m=m
        self.sPipe=sPipe
        self.sSeabed=sSeabed
        self.initialT=initialT
        self.finalT=finalT
        self.jobName=jobName
        self.w0=w0
        self.L_b=L_b
        self.P_u=P_u
        self.mu_L=mu_L

    def createPerfectPipe(self):
        '''
        create a part for a pipe without initial imperfection
        The part will used for the following calcualtion.
        name: the name of pipe part
        :return:
        '''
        import part
        vp=session.currentViewportName
        m=mdb.Model(self.pipeModel)
        s=m.Sketch(name='__profile__',sheetSize=self.L*1.5)
        s.Line(point1=(0,0),point2=(self.L,0))
        p=m.Part(name=self.pipePart,dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p.BaseWire(sketch=s)

    def createPerfectPipe4HT(self):
        '''
        create a pipe without initial imperfection. Extra pipe length is used for a strong boundary condition.
        the extra length in both ends is 5% L, respectively.
        '''
        import part
        vp=session.currentViewportName
        m=mdb.Model(self.pipeModel)
        s=m.Sketch(name='__profile__',sheetSize=self.L*1.5)
        s.Line(point1=(-0.05*self.L,0),point2=(0,0))
        s.Line(point1=(0,0),point2=(self.L,0))
        s.Line(point1=(self.L,0),point2=(1.05*self.L,0))
        p=m.Part(name=self.pipePart,dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p.BaseWire(sketch=s)


    def createSeabed(self):
        '''
        produce a  part for rigid seabed for the seabed simplification.
        the overal length of seabed is 0.94*L
        :return:
        '''
        import part
        m=mdb.models[self.pipeModel] # this part should follow after the pipe part.
        s=m.Sketch(name='__profile__',sheetSize=self.L*1.5)
        s.rectangle(point1=(0,0),point2=(self.L,self.w_bed))
        p=m.Part(name=self.seabedPart,dimensionality=THREE_D,type=DISCRETE_RIGID_SURFACE)
        p.BaseShell(sketch=s)

    def pipePro(self):
        '''
        this function is to define the common properties of pipes with respect to material section, assembly,
        and mesh.
        :return:
        '''
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]

        import section
        material_lib(self.pipeModel)
        m.PipeProfile(name='profile-1', r=self.od/2, t=self.t)
        m.BeamSection(name='Section-1',profile='profile-1', material=self.m,integration=DURING_ANALYSIS)
        pipeEdge=p.edges
        region=(pipeEdge,)
        p.SectionAssignment(region=region, sectionName='Section-1')
        region=(pipeEdge,)
        p.assignBeamSectionOrientation(region=region,method=N1_COSINES,n1=(0,0,-1))
        region=(pipeEdge,)
        p.assignBeamSectionOrientation(region=region,method=N1_COSINES,n1=(0,0,-1))

        import assembly
        a=m.rootAssembly
        a.Instance(name=self.pipePart+str(-1),part=p,dependent=OFF)
        # edges for seeds and mesh
        l_lon=a.instances[self.pipePart+str(-1)].edges.getByBoundingBox(xMin=-1,xMax=self.L+1)
        a.Set(edges=l_lon,name='pipeL')

        import mesh
        a.seedEdgeBySize(edges=l_lon,size=self.sPipe,constraint=FINER)
        regions=a.instances[self.pipePart+str(-1)].edges
        a.generateMesh(regions=(a.instances[self.pipePart+str(-1)],))
        elemType1=mesh.ElemType(elemCode=PIPE32H)
        regions=regionToolset.Region(edges=a.instances[self.pipePart+str(-1)].edges)
        a.setElementType(regions=regions,elemTypes=(elemType1,))

    def seabedPro(self):

        m=mdb.models[self.pipeModel]
        p=m.parts[self.seabedPart]

        import section
        p.ReferencePoint(point=(self.L/2,self.w_bed/2,0))

        import assembly
        a=m.rootAssembly
        a.Instance(name=self.seabedPart+str(-1),part=p,dependent=OFF)

        import mesh
        # for 3D discrete rigid structures,
        # the default element type is R3D4  in ABAQUS
        l_lon=a.instances[self.seabedPart+str(-1)].edges.findAt(((self.sSeabed,0,0),),)
        l_width=a.instances[self.seabedPart+str(-1)].edges.findAt(((0,self.sSeabed,0),),)
        a.Set(name="l_lon",edges=l_lon)
        a.seedEdgeBySize(edges=l_lon,size=self.sSeabed,constraint=FINER)
        a.seedEdgeBySize(edges=l_width,size=self.sSeabed,constraint=FINER)
        a.generateMesh(regions=(a.instances[self.seabedPart+str(-1)],))

    def pipeEigenValue(self):
        '''
        Before the call of this method, the pipe Part should be produced by calling
        createPerfectPipe function. The eigen mode could used for pipe imperfections in the following simulation
        :return:
        '''
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        self.pipePro()
        nodes=a.instances[self.pipePart+str(-1)].nodes
        n_all=nodes.getByBoundingCylinder((-1,0,0),(self.L+1,0,0),0.2)
        n_bottom=nodes.getByBoundingBox(xMin=-1,xMax=0.1)
        n_top=nodes.getByBoundingBox(xMin=self.L-0.1,xMax=self.L+0.1)
        a.Set(nodes=n_all,name='n_all')
        a.Set(nodes=n_bottom,name='n_bottom')
        a.Set(nodes=n_top,name='n_top')

        import step
        m.BuckleStep(name='Step-eigenValue',previous='Initial',description='EigenValue analyses',
                 numEigen=5,eigensolver=SUBSPACE)

        import load
        regions=regionToolset.Region(nodes=n_top)
        m.DisplacementBC(name='top-Simply',createStepName='Initial',region=regions,u2=SET)
        regions=regionToolset.Region(nodes=n_bottom)
        m.DisplacementBC(name='bottom',createStepName='Initial',region=regions,u1=SET,u2=SET)
        regions=regionToolset.Region(nodes=n_all)
        m.DisplacementBC(name='upheaval-fix',createStepName='Initial',region=regions,u3=SET)
        regions=regionToolset.Region(nodes=n_top)
        m.ConcentratedForce(name='unitLoad',createStepName='Step-eigenValue',region=regions,cf1=-1.0)


    def pipeStaticDamping(self):
        '''
        Before the call of this method, the pipe Part and seabed Part should be produced by calling
        functions, and the common modeling function pipePro() should be called
        '''
        m=mdb.models[self.pipeModel]
        p_pipe=m.parts[self.pipePart]
        p_seabed=m.parts[self.seabedPart]
        a=m.rootAssembly

        self.pipePro()
        self.seabedPro()
        a.translate(instanceList=(self.seabedPart+str(-1), ), vector=(0, -self.w_bed/2, 0))
        #a.translate(instanceList=(self.seabedPart+str(-1), ), vector=(0.0, 0.0, -self.od/2))  #-self.od/2-1

        nodes_pipe=a.instances[self.pipePart+str(-1)].nodes
        pipeNode_all=nodes_pipe.getByBoundingCylinder((-0.05*self.L-0.1,0,0),(1.05*self.L+0.1,0,0),0.2)
        pipeNode_bottomBC=nodes_pipe.getByBoundingBox(xMin=-0.05*self.L-0.1,xMax=-0.1)
        pipeNode_bottom=nodes_pipe.getByBoundingBox(xMin=-0.1,xMax=0.1)
        pipeNode_topBC=nodes_pipe.getByBoundingBox(xMin=self.L+0.1,xMax=1.05*self.L+0.1)
        pipeNode_top=nodes_pipe.getByBoundingBox(xMin=self.L-0.1,xMax=self.L+0.1)
        pipeNode_center=nodes_pipe.getByBoundingBox(xMin=self.L/2-0.1,xMax=self.L/2+0.1)
        a.Set(nodes=pipeNode_all,name='pipeNode_all')
        a.Set(nodes=pipeNode_bottomBC,name='pipeNode_bottomBC')
        a.Set(nodes=pipeNode_bottom,name='pipeNode_bottom')
        a.Set(nodes=pipeNode_topBC,name='pipeNode_topBC')
        a.Set(nodes=pipeNode_top,name='pipeNode_top')
        a.Set(nodes=pipeNode_center,name='pipeNode_center')


        ref1=a.instances[self.seabedPart+str(-1)].referencePoints  # only one reference is allowed.
        refPoints1=(ref1[ref1.keys()[0]],)
        a.Set(referencePoints=refPoints1, name=self.seabedPart+'refPoint')

        import step
        #m.StaticStep(name='contact',previous='Initial',description='building contact',initialInc=0.1,nlgeom=ON)
        m.StaticStep(name='pre-load',previous='Initial',\
                         description='this is the step used for pre-loading of pipe',\
                         initialInc=0.1,nlgeom=ON,minInc=1e-9)
        m.StaticStep(name='buckling',previous='pre-load',\
                         description='this is the buckling step used for buckling of pipe',\
                         stabilizationMagnitude=1e-9,stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
                        adaptiveDampingRatio=0.05,nlgeom=ON,
                        initialInc=0.001,minInc=1e-09,maxInc=0.01,maxNumInc=1000)
        m.FieldOutputRequest(name='F_extra',createStepName='pre-load',variables=('NFORC','CFORCE','CDISP',
                                                                           'CSTRESS','NT'))
        regions=a.sets['pipeNode_top']
        m.HistoryOutputRequest(name='n_top',createStepName='pre-load',variables=('U1','CFN','NT'),region=regions)
        regions=a.sets['pipeNode_center']
        m.HistoryOutputRequest(name='pipeNode_center',createStepName='pre-load',variables=('U2','CFN','NT'),region=regions)

        import interaction
        ###### the offset amount "1" should be consistant with the offshore of seabed instance in the former section.
        #####
        s1 = a.instances[self.seabedPart+str(-1)].faces
        side1Faces1 = s1.findAt(((self.L/2, -self.w_bed/2, 0), ))
        a.Surface(side1Faces=side1Faces1, name='seabed-sf')
        c1 = a.instances[self.pipePart+str(-1)].edges
        circumEdges1 = c1.findAt(((self.L/2, 0.0, 0.0), ))
        a.Surface(circumEdges=circumEdges1, name='edge_sf')
        m.ContactProperty('IntProp-1')
        m.interactionProperties['IntProp-1'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
            self.mu_L, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION,
            fraction=0.005, elasticSlipStiffness=None)
        m.interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD)
        region1=a.Surface(side1Faces=side1Faces1, name='seabed-sf')
        region2=a.Surface(circumEdges=circumEdges1, name='edge_sf')
        m.SurfaceToSurfaceContactStd(name='Int-1',createStepName='Initial', master=region1, slave=region2,
            sliding=FINITE,interactionProperty='IntProp-1',enforcement=NODE_TO_SURFACE)

        import load
        regions=regionToolset.Region(nodes=pipeNode_topBC)
        m.DisplacementBC(name='pipeNode_topBC',createStepName='Initial',region=regions,u1=SET,u2=SET,u3=SET)
        # m.boundaryConditions['pipeTop'].setValuesInStep(stepName='contact', u3=-1.01)
        # m.boundaryConditions['pipeTop'].setValuesInStep(stepName='pre-load', u3=SET)
        #m.PinnedBC(name='pipeNode_top',createStepName='Initial',region=regions)

        regions=regionToolset.Region(nodes=pipeNode_bottomBC)
        m.DisplacementBC(name='pipeNode_bottomBC',createStepName='Initial',region=regions,u1=SET,u2=SET,u3=SET)

        regions=regionToolset.Region(referencePoints=refPoints1)
        m.EncastreBC(name='seabed',createStepName='Initial',region=regions)

        #m.boundaryConditions['seabed'].deactivate('pre-load')
        #regions=regionToolset.Region(referencePoints=refPoints1)
        #m.EncastreBC(name='seabed-1',createStepName='pre-load',region=regions)

        # regions=a.sets['pipeNode_all']
        # m.Temperature(name='initialTem',createStepName='Initial',region=regions,magnitudes=(self.initialT,))
        pipe_nodes_WOends=a.instances[self.pipePart+str(-1)].nodes.getByBoundingBox(xMin=-0.1,xMax=self.L+0.1)
        a.Set(nodes=pipe_nodes_WOends,name='pipe_nodes_WOends')
        regions=a.sets['pipe_nodes_WOends']
        m.Temperature(name='Temp',createStepName='Initial',region=regions,magnitudes=(self.initialT,))
        m.predefinedFields['Temp'].setValuesInStep(stepName='buckling',magnitudes=(self.finalT,))

        pipeEdge_lon=a.instances[self.pipePart+str(-1)].edges.getByBoundingBox(xMin=-0.1,xMax=self.L+0.1)
        region = regionToolset.Region(edges=pipeEdge_lon)
        m.LineLoad(name='lineLoad', createStepName='pre-load',region=region, comp3=-self.P_u)
        #m.loads['lineLoad'].deactivate('buckling')

    def pipeRiks_seabed(self):

        m=mdb.models[self.pipeModel]
        p_pipe=m.parts[self.pipePart]
        p_seabed=m.parts[self.seabedPart]
        a=m.rootAssembly

        m.setValues(absoluteZero=-273.15)
        self.pipeStaticDamping()

        m.StaticRiksStep(name='buckling',previous='pre-load',maintainAttributes=True,maxLPF=2.0,maxNumInc=1000,
                         initialArcInc=0.001)
        m.historyOutputRequests['n_top'].suppress()
        m.historyOutputRequests['pipeNode_center'].suppress()
        regions=a.sets['pipeNode_top']
        m.HistoryOutputRequest(name='n_top1',createStepName='buckling',variables=('U1','NT'),region=regions)
        regions=a.sets['pipeNode_center']
        m.HistoryOutputRequest(name='n_center1',createStepName='buckling',variables=('U2','NT'),region=regions)


    def pipeRiks(self):
        '''
        this function can be simply used for the buckling analysis based on modified Riks mehtod.
        Pipe is under uniaxial compression, There is no introduction of the seabed.
        Before the call of this method, the pipe Part should be produced by calling
        createPerfectPipe function, the common modeling function pipePro() should be called as well
        :return: no return, but produce a model for riks simulation.
        '''
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        self.pipePro()
        nodes=a.instances[self.pipePart+str(-1)].nodes
        n_all=nodes.getByBoundingCylinder((-1,0,0),(self.L+1,0,0),0.2)
        n_bottom=nodes.getByBoundingBox(xMin=-1,xMax=0.1)
        n_top=nodes.getByBoundingBox(xMin=self.L-0.1,xMax=self.L+0.1)
        n_center=nodes.getByBoundingBox(xMin=self.L/2-0.1,xMax=self.L/2+0.1)
        a.Set(nodes=n_all,name='n_all')
        a.Set(nodes=n_bottom,name='n_bottom')
        a.Set(nodes=n_top,name='n_top')
        a.Set(nodes=n_center,name='n_center')

        import interaction

        import step
        m.StaticRiksStep(name='Step-Riks',previous='Initial',\
                         description='this is the step used for stability analysis',\
                         maxLPF=1,maxNumInc=100,initialArcInc=0.01,minArcInc=1e-10,maxArcInc=0.1,nlgeom=ON)
        regions=a.sets['n_top']
        m.HistoryOutputRequest(name='H-Output-2',createStepName='Step-Riks',variables=('UT','MISES','S11','E11'),region=regions)
        regions=a.sets['n_center']
        m.HistoryOutputRequest(name='centerHistory',createStepName='Step-Riks',variables=('UT',),region=regions)

        import load
        regions=regionToolset.Region(nodes=n_top)
        m.DisplacementBC(name='top-Simply',createStepName='Initial',region=regions,u2=SET)
        regions=regionToolset.Region(nodes=n_bottom)
        m.DisplacementBC(name='bottom',createStepName='Initial',region=regions,u1=SET,u2=SET)
        regions=regionToolset.Region(nodes=n_all)
        m.DisplacementBC(name='upheaval-fix',createStepName='Initial',region=regions,u3=SET)
        # regions=regionToolset.Region(nodes=n_top)
        # m.PinnedBC(name='BC-top',createStepName='Initial',region=regions)
        # regions=regionToolset.Region(nodes=n_bottom)
        # m.PinnedBC(name='BC-bottom',createStepName='Initial',region=regions)
        # regions=a.sets['pipeL']
        # m.Temperature(name='initialTem',createStepName='Initial',region=regions,magnitudes=(self.initialT,))
        # m.predefinedFields['initialTem'].setValuesInStep(stepName='Step-Riks',magnitudes=(self.finalT,))


    def jobEigen(self):
        import job
        mdb.Job(name='pipeEigen',model=self.pipeModel,description='Eigenvalue analyese starts')
        #mdb.jobs['pipeEigen'].submit()

    def jobRiks(self):
        import job
        mdb.Job(name=self.jobName, model=self.pipeModel, description='Stability analyses based on Riks method starts',\
                numCpus=6, numDomains=6)
        #mdb.jobs[self.jobName].submit()

    def jobStaticDamping(self):
        import job
        mdb.Job(name=self.jobName,model=self.pipeModel,description="Buckling based on static damping",numCpus=6,
                numDomains=6)

    def jobRiks_seabed(self):
        import job
        mdb.Job(name=self.jobName,model=self.pipeModel,description="Buckling based on Riks with seabed contact",
                numCpus=6,numDomains=6)

    def initImperfectionMode1(self):
        '''
        this method is used to produce a real initial imperfection  mode 1 via offset of the pipe nodes
        the shape function is based on buckling modes, as illustrated in script: pipeCalculation.py
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        #regions=regionToolset.Region(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]
        num=4.493*2/self.L_b   # nL=8.986

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.
            dy2=self.w0/15.7*(1+8.9868**2/8-(num*x0)**2/2-np.cos(num*x0)/np.cos(8.9868/2))  # equation from Taylor, 1986
            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

    def initImperfectionMode2(self):
        '''
        the shape function is based on buckling modes, as illustrated in script: pipeCalculation.py.
        The function derivation is from Taylor, 1986
        the definition of buckling length refers to Paper: pipeReview2019
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-self.L_b-0.2,0,0),(0.5*self.L+self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.

            L_b_half=self.L_b/2
            if x0>=0:
                dy2=-self.w0/8.62 *(1-np.cos(2*np.pi*x0/L_b_half)+\
                                np.pi*np.sin(2*np.pi*x0/L_b_half)+(2*np.pi**2*x0/L_b_half)*(1-x0/L_b_half))
            else:
                dy2=self.w0/8.62 *(1-np.cos(2*np.pi*x0/L_b_half)-\
                                    np.pi*np.sin(2*np.pi*x0/L_b_half)-(2*np.pi**2*x0/L_b_half)*(1+x0/L_b_half))
            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

    def initImperfectionMode3(self):
        '''
        the shape function of imperfection is based on buckling modes, as illustrated in script: pipeCalculation.py.
        The function derivation is based on Kerr,1979
        the definition of buckling length refers to Paper: pipeReview2019
        L_b here is changed to the entire length of buckling, instead of
        the central lobe length . Noted that for mode 3, extra lobes are included.
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        # select the region for L_b. L_b here is changed to the entire length of buckling, instead of
        # the central lobe length . Noted that for mode 3, extra lobes are included.
        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]
        n3=2*7.551/self.L_b
        l1=self.L_b/2.59

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.

            if x0>=-self.L_b/2/2.59 and x0<=self.L_b/2/2.59:
                dy2=self.w0/1.484*(1+0.484*np.cos(n3*x0)-2.109*(x0/l1)**2)
            elif x0>self.L_b/2/2.59 and x0<=self.L_b/2:
                dy2=self.w0/0.821*(1+0.134*np.cos(5.836*x0/l1)+0.031*np.sin(5.836*x0/l1)-\
                                   2.341*(x0/l1)+1.167*(x0/l1)**2)
            else:
                dy2=self.w0/0.821*(1+0.134*np.cos(5.836*x0/l1)-0.031*np.sin(5.836*x0/l1)+\
                                   2.341*(x0/l1)+1.167*(x0/l1)**2)
            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

    def initImperfectionMode4(self):
        '''
        Directly introduce structural imperfections.
        the shape function of imperfection is based on buckling modes, as illustrated in script: pipeCalculation.py.
        The function derivation is based on Kerr,1979
        the definition of buckling length refers to Paper: pipeReview2019
        L_b here is changed to the entire length of buckling. Noted that for mode 4, extra
        lobes are included.
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        # select the region for L_b. L_b here is changed to the entire length of buckling. Noted that for mode 34, extra
        # arc is included.
        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]
        l=self.L_b/2
        n4=8.54/l

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.

            if x0>=0 and x0<=l/1.61:
                dy2=self.w0/8.2*(1-np.cos(8.54*x0/l)+3*np.sin(8.54*x0/l)+25.74*(x0/l)-36.11*(x0/l)**2)
            elif x0>l/1.61 and x0<=l:
                dy2=self.w0/0.304*(1+0.00391*np.cos(8.54*x0/l)+0.05078*np.sin(8.54*x0/l)-2.375*(x0/l)+1.3398*(x0/l)**2)
            elif x0<0 and x0>=-l/1.61:
                dy2=-self.w0/8.2*(1-np.cos(8.54*x0/l)-3*np.sin(8.54*x0/l)-25.74*(x0/l)-36.11*(x0/l)**2)
            else:
                dy2=-self.w0/0.304*(1+0.00391*np.cos(8.54*x0/l)-0.05078*np.sin(8.54*x0/l)+2.375*(x0/l)+1.3398*(x0/l)**2)

            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

def main():
    dataInput=DataInput()
    od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,sPipe,sSeabed,initialT,finalT,\
    jobName,flagEigen,w0,L_b,P_u,mu_L=dataInput.dataInput()
    ## creat an input data based for GUI inputs
    lateralPipeBucklingInst=LateralPipeBuckling(od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,sPipe,sSeabed,\
                                                initialT,finalT,jobName,w0,L_b,P_u,mu_L)
    # lateralPipeBucklingInst.createPerfectPipe()
    # lateralPipeBucklingInst.createSeabed()
    if flagEigen==1:
        lateralPipeBucklingInst.createPerfectPipe()
        lateralPipeBucklingInst.pipeEigenValue()
        lateralPipeBucklingInst.jobEigen()
    else:
        lateralPipeBucklingInst.createPerfectPipe4HT()
        lateralPipeBucklingInst.createSeabed()
        #lateralPipeBucklingInst.pipeRiks_seabed()
        lateralPipeBucklingInst.pipeStaticDamping()
        lateralPipeBucklingInst.initImperfectionMode2()
        # lateralPipeBucklingInst.jobStaticDamping()

if __name__ == '__main__':
    main()

