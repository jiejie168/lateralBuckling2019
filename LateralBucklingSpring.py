__author__ = 'Jie'
#filename: LateralBucklingSpring.py
'''
This script is used for the produce of pipe model with a fixed buckling length
the pipe-seabed interaction is modeled by a spring connector.
No seabed is needed.

# for case with imperfection 6, the L & Lb should be set to half for inputting.
# it should use combining with the "pipeSpringRiks_FreeRotation_Bottom_symm.py" model.
'''
import numpy as np
from abaqus import *
from abaqusConstants import *
import testUtils
testUtils.setBackwardCompatibility()
import regionToolset
import displayGroupMdbToolset as dgm
from material_libLateralBuckling import material_lib
import connectorBehavior
# from scipy.optimize import root, fsolve
from os import environ      # the aim of these import is try to connect pycharm with abaqus libraries.
environ['ABAQUS_BAT_PATH'] = 'C:\\SIMULIA\\Abaqus\\Commands\\abaqus' # allowing pycharm IDE check types and run abaqus python program
environ['ABAQUS_BAT_SETTING'] = 'noGUI'

class DataInput():
    '''
    This class is only used for the input data of simulation.
    default data come from Zeng,2014 (Critical upheaval buckling forces of imperfect pipelines).
    zengData: od=457, t=14.3, L=1200m,seabedWidth=40m, temp_initial=20oC,temp_final=120oC,q=1.5N/mm;u_L=0.6
    zeng: w0=200mm, L_b0=100m;
    Data from McCarron, 2019:
    default data from Philippe: od=324mm,t=17mm,L_b=L=100m,stiffness density (k, N/mm2) : 5248e-6....

    in this scripts, the L should equal to L_b
    K=131.2e-6; N=64853.4 N ;a=2.235 and b=3.600
    K=13120e-6; N=728719.2 N; a=6.709 and b=11.598

    a and b are the coefficients of imperf6. It should be calculated in advance, via "PhilippeEquations.py".
    The default values of a=2.235 and b=3.600 are for the case with K=131.2 N/mm2.

    '''
    def dataInput(self):
        db1=(('pipe outer Diameter(mm):','324'),('pipe thickness(mm):','17'),('pipe length(mm):','100000'),
             ('pipe part:','pipe1'),('seabed part:','seabed1'),('seabed width(mm):','40000'),
             ('pipe model:','model1'),('pipeMaterial:','Steel-e'),('refLoad(N,only valid for Riks)','64853.4'))
        db2=(('pipe mesh size(mm):','1000'),('seabed mesh size(mm):','20000'),('initial temperature(oC):','20'),
             ('final temperature(oC):','120'),('job name:','job1-SD'),('flagEigen','0'))
        db3=(('w0(mm)','1000'),('L_b(mm)','100000'),('a','6.709'),('b','11.598'))   # s=0.01. a and b are the coefficients of imperf6. It should be calculated in advance.
        db4=(("unit downard force(N/mm):",'0.664'),("Stiffness density K(N/mm2)",'131.2e-6'),
             ("frictional coefficient(mu_L)",'0.95'),("Young's Modulus(E,MPa)",'2.1e5'))
        flag1=True
        while flag1:
            od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,Pref=getInputs(fields=db1, label='Please enter the pipe dimensions')
            try:
                od,t,L,w_bed=float(od),float(t),float(L),float(w_bed)
            except (TypeError,ValueError):
                raise ValueError,'please input a float!!'
            Pref=float(Pref)
            pipePart,seabedPart,pipeModel,pipeMaterial=str(pipePart),str(seabedPart),str(pipeModel),str(pipeMaterial)
            if bool(od)==False or bool(t)==False or bool(L)==False or bool(w_bed)==False or bool(Pref)==False or\
                            bool(pipePart)==False or bool(seabedPart)==False or bool(pipeModel)==False:
                flag1=True
            else:
                break
        print ("input of pipe dimension completed")

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
            w0,L_b,a,b=getInputs(fields=db3,label='Please enter the imperfection:')
            try:
                w0,L_b,a,b=float(w0),float(L_b),float(a),float(b)
            except (ValueError,TypeError):
                raise ValueError,'please input a float for imperfection value'
            if  bool(L_b)==False or bool(a)==False or bool(b)==False:
                flag3=True
            else:
                break
        print ("input of imperfection and other information completed")

        flag4=True
        while flag4:
            P_u,K,mu_L,E=getInputs(fields=db4,label="Please enter other information:")
            try:
                P_u,K,mu_L,E=float(P_u),float(K),float(mu_L),float(E)
            except (ValueError,TypeError):
                raise ValueError,'please input a float for effective load and friction'

            if bool(P_u)==False or bool(K)==False or bool(mu_L)==False or bool(E)==False:
                flag4=True
            else:
                break
        print("input of other information completed")

        return od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,Pref,sPipe,sSeabed,initialT,finalT,jobName,flagEigen,w0,L_b,P_u,K,mu_L,E,a,b

class LateralBucklingSpring():

    def __init__(self,od,t,L,pipePart,seabedPart,w_bed,pipeModel,m,Pref,sPipe,sSeabed,initialT,finalT,jobName,w0,L_b,P_u,K,mu_L,E,a,b):
        self.od=od
        self.t=t
        self.L=L
        self.pipePart=pipePart
        self.seabedPart=seabedPart
        self.w_bed=w_bed
        self.pipeModel=pipeModel
        self.m=m
        self.Pref=Pref
        self.sPipe=sPipe
        self.sSeabed=sSeabed
        self.initialT=initialT
        self.finalT=finalT
        self.jobName=jobName
        self.w0=w0
        self.L_b=L_b
        self.P_u=P_u
        self.K=K
        self.mu_L=mu_L
        self.E=E
        self.a=a
        self.b=b

    def pipeCal(self):
        d=self.od-2*self.t
        alpha1=d/self.od
        I=np.pi*self.od**4*(1-alpha1**4)/64
        J=2*I
        A=(1-alpha1**2)*np.pi*self.od**2/4
        # ibxilong_temp=self.alpha*self.deltaT
        # ibxilong_pressure=self.deltaP*self.D*(1-2*0.3)/(4*self.t*self.E)
        # P0=self.E*A*(ibxilong_temp+ibxilong_pressure)
        # L_expand=self.L*(ibxilong_temp+ibxilong_pressure)

        G=self.E/2.0/(1+0.3)
        alpha=I/A/(self.od**2/4)
        beta=G/self.E
        gama=A*self.L_b**2/I
        xi=self.K*(self.od**2/4)/self.E/A

        return I,J,A,G,alpha,beta,gama,xi

    def resultsTorsion_fun(self,Pcn):

        I,J,A,G,alpha,beta,gama,xi=self.pipeCal()
        # a1=np.sqrt(4*gama*np.sqrt(2*xi*(2*alpha*beta-Pcn)/beta)+gama*(xi-2*beta*Pcn)/beta)/4
        # b1=np.sqrt(4*self.gama*np.sqrt(2*self.xi*(2*self.alpha*self.beta-Pcn)/self.beta)-self.gama*(self.xi-2*self.beta*Pcn)/self.beta)/4

        a=0.5/np.sqrt(2)*np.sqrt(gama/beta)*np.sqrt(2*np.sqrt(2*beta*xi*(2*alpha*beta-Pcn))+xi-2*beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(gama/beta)*np.sqrt(2*np.sqrt(2*beta*xi*(2*alpha*beta-Pcn))-xi+2*beta*Pcn)
        func=b*(-alpha*xi*gama**2+Pcn*gama*(a**2+b**2)+3*a**4+2*(a**2)*(b**2)-b**4)*np.sinh(a)+\
             a*(alpha*xi*gama**2+Pcn*gama*(a**2+b**2)+a**4-2*(a**2)*(b**2)-3*b**4)*np.sin(b)
        return np.array(func)

    def find_Pc_freeTorsion(self,guess_value):

        I,J,A,G,alpha,beta,gama,xi=self.pipeCal()
        Pcn=fsolve(self.resultsTorsion_fun,guess_value)
        r=root(self.resultsTorsion_fun,guess_value)
        print ("The solving of P0 is: {}".format(r.success))
        a=0.5/np.sqrt(2)*np.sqrt(gama/beta)*np.sqrt(2*np.sqrt(2*beta*xi*(2*alpha*beta-Pcn))+xi-2*beta*Pcn)
        b=0.5/np.sqrt(2)*np.sqrt(gama/beta)*np.sqrt(2*np.sqrt(2*beta*xi*(2*alpha*beta-Pcn))-xi+2*beta*Pcn)

        temp1=(2*a*b*np.cosh(a/2)*np.cos(b/2)+(a**2-b**2)*np.sinh(a/2)*np.sin(b/2))/(2*a*b*(np.cosh(a/2)**2*np.cos(b/2)**2+np.sinh(a/2)**2*np.sin(b/2)**2))
        temp2=((a**2-b**2)*np.cosh(a/2)*np.cos(b/2)-2*a*b*np.sinh(a/2)*np.sin(b/2))/(2*a*b*(np.cosh(a/2)**2*np.cos(b/2)**2+np.sinh(a/2)**2*np.sin(b/2)**2))
        x=np.linspace(-0.5,0.5,num=1000)
        w=1-temp1*np.cosh(a*x)*np.cos(b*x)+temp2*np.sinh(a*x)*np.sin(b*x)

        return Pcn, a,b,temp1,temp2

    def createPerfectPipe(self):
        '''
        create a part for a pipe without initial imperfection
        The part will used for the following calculation.
        name: the name of pipe part.
        :return:
        '''
        import part
        vp=session.currentViewportName
        m=mdb.Model(self.pipeModel)
        s=m.Sketch(name='__profile__',sheetSize=self.L*1.5)
        s.Line(point1=(0,0),point2=(self.L,0))
        p=m.Part(name=self.pipePart,dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p.BaseWire(sketch=s)

    def createPerfectPipeWithPartition(self):
        '''
        this part is only used for the spring connector purporse.
        :return:
        '''
        import part
        m=mdb.Model(self.pipeModel)
        s=m.Sketch(name='__profile__',sheetSize=self.L*1.5)
        num_elem=np.round(self.L/self.sPipe)

        for i in range(num_elem):
            s.Line(point1=(self.sPipe*i,0),point2=(self.sPipe*(i+1),0))
        p=m.Part(name=self.pipePart,dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p.BaseWire(sketch=s)

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
        a.Instance(name=self.pipePart+str(-2),part=p,dependent=OFF)
        ins_seabed=a.instances[self.pipePart+str(-2)]
        ins_seabed.translate(vector=(0,0,-100))
        # edges for seeds and mesh
        l_lon=a.instances[self.pipePart+str(-1)].edges.getByBoundingBox(xMin=-1,xMax=self.L+1)
        a.Set(edges=l_lon,name='pipeL')
        l_lon_seabed=a.instances[self.pipePart+str(-2)].edges.getByBoundingBox(xMin=-1,xMax=self.L+1,zMin=-101,zMax=-90)
        a.Set(edges=l_lon_seabed,name='pipeL_seabed')

        import mesh
        a.seedEdgeBySize(edges=l_lon+l_lon_seabed,size=self.sPipe,constraint=FINER)
        #regions=a.instances[self.pipePart+str(-1)].edges
        a.generateMesh(regions=(a.instances[self.pipePart+str(-1)],a.instances[self.pipePart+str(-2)]))
        elemType1=mesh.ElemType(elemCode=PIPE32H)
        regions=regionToolset.Region(edges=a.instances[self.pipePart+str(-1)].edges+a.instances[self.pipePart+str(-2)].edges)
        a.setElementType(regions=regions,elemTypes=(elemType1,))

    def SpringConnector(self):
        '''
        This one is used to create the property of connector, used for the mimic of seabed effect. It should be noted that
        this simplification is only suitable for a small lateral displacement. The typical breakout distance for pipe is
        0.1D.
        Adoped the linear coupled connector. Using CARTISIAN coordiante system and Euler system for rotation.
        '''
        import interaction
        import section
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly

        ###################################################
        #Create Wires for the definition of connector.
        ###################################################
        v_all=a.instances[self.pipePart+str(-1)].vertices
        v_center=v_all.getByBoundingBox(xMin=self.sPipe-1,xMax=self.L+1-self.sPipe,zMin=-1,zMax=1)
        v_ends=v_all.findAt(((0,0,0),),((self.L,0,0),))
        v_all_seabed=a.instances[self.pipePart+str(-2)].vertices
        v_center_seabed=v_all_seabed.getByBoundingBox(xMin=self.sPipe-1,xMax=self.L+1-self.sPipe,zMin=-101,zMax=-99)
        v_ends_seabed=v_all_seabed.findAt(((0,0,-100),),((self.L,0,-100),))
        #length_ends=len(v_ends)
        #None_value=(None for i in range(length_ends))
        points_ends=zip(v_ends_seabed,v_ends)
        wire=a.WirePolyLine(points=points_ends,mergeWire=OFF,meshable=OFF)
        points_center=zip(v_center_seabed,v_center)
        wire2=a.WirePolyLine(points=points_center,mergeWire=OFF,meshable=OFF)
        e1=a.edges
        edges_center=e1.getByBoundingBox(xMin=0.1,xMax=self.L-0.1)
        a.Set(edges=edges_center,name='wire_center')
        edges_ends=e1.findAt(((0,0,-1),),((self.L,0,-1),))
        a.Set(edges=edges_ends,name='wire_ends')

        ###################################################
        #Define the connector section properties.
        ###################################################
        m.ConnectorSection(name='connectorSec_center',translationalType=CARTESIAN,rotationalType=EULER)
        elastic_center=connectorBehavior.ConnectorElasticity(components=(2,5),coupling=COUPLED,
                table=((self.K*self.sPipe,-self.K*self.sPipe*self.od/2,self.K*self.sPipe*self.od**2/4),))
        elastic_center.ConnectorOptions()
        m.sections["connectorSec_center"].setValues(behaviorOptions=(elastic_center,))

        m.ConnectorSection(name='connectorSec_end',translationalType=CARTESIAN,rotationalType=EULER)
        # in order to define a uniform distribution of connector. At ends, the spring stiffness is half.
        elastic_ends=connectorBehavior.ConnectorElasticity(components=(2,5),coupling=COUPLED,
            table=((self.K*self.sPipe/2,-self.K*self.sPipe*self.od/2/2,self.K*self.sPipe*self.od**2/4/2),))
        elastic_ends.ConnectorOptions()
        m.sections["connectorSec_end"].setValues(behaviorOptions=(elastic_ends,))

        ##########################################################################################
        #Assign the connectorSection Properties to the Wires.
        ##########################################################################################
        a.DatumCsysByThreePoints(name='d1',coordSysType=CARTESIAN,origin=(0,0,0),point1=(1,0,0),
                                 point2=(0,1,0))
        datum=a.datums
        key1=datum.keys()
        key1=key1[0]
        datum1=datum[key1]
        region=a.sets['wire_center']
        csa=a.SectionAssignment(sectionName='connectorSec_center',region=region)
        a.ConnectorOrientation(region=csa.getSet(),localCsys1=datum1)
        region=a.sets['wire_ends']
        csa1=a.SectionAssignment(sectionName='connectorSec_end',region=region)
        a.ConnectorOrientation(region=csa1.getSet(),localCsys1=datum1)


    def pipeSpringEigen_FreeRotation_Bottom(self):
        '''
        This function is used for the linear analysis of pipe with free rotation at pipeEnds, The connection is
        equal to the connection with bottom line. X direction is the pipe longitudinal direction; Y direction is the
        pipe lateral direction.  The displacement in Z direction (upheaval) is fixed. Simplify-supported Bc.
        :return:
        '''
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        self.pipePro()
        self.SpringConnector()
        nodes=a.instances[self.pipePart+str(-1)].nodes
        n_all=nodes.getByBoundingCylinder((-1,0,0),(self.L+1,0,0),0.2)
        n_bottom=nodes.getByBoundingBox(xMin=-1,xMax=0.1)
        n_top=nodes.getByBoundingBox(xMin=self.L-0.1,xMax=self.L+0.1)
        n_center=nodes.getByBoundingBox(xMin=self.L/2-0.1,xMax=self.L/2+0.1)
        n_side=nodes.getByBoundingBox(xMin=self.L/5-0.1,xMax=self.L/5+0.1)
        a.Set(nodes=n_all,name='n_all')
        a.Set(nodes=n_bottom,name='n_bottom')
        a.Set(nodes=n_top,name='n_top')
        a.Set(nodes=n_center,name='n_center')
        a.Set(nodes=n_side,name='n_side')

        nodes_seabed=a.instances[self.pipePart+str(-2)].nodes
        a.Set(nodes=nodes_seabed,name="nodes_seabed")

        import step
        m.BuckleStep(name='Step-eigenValue',previous='Initial',description='EigenValue analyses',
                 numEigen=5,eigensolver=SUBSPACE,maxIterations=100)

        import load
        regions=regionToolset.Region(nodes=n_top)
        m.DisplacementBC(name='top-Simply',createStepName='Initial',region=regions,u2=SET)
        regions=regionToolset.Region(nodes=n_bottom)
        m.DisplacementBC(name='bottom',createStepName='Initial',region=regions,u1=SET,u2=SET)
        regions=regionToolset.Region(nodes=n_all)
        m.DisplacementBC(name='upheaval-fix',createStepName='Initial',region=regions,u3=SET)
        regions=regionToolset.Region(nodes=n_top)
        m.ConcentratedForce(name='unitLoad',createStepName='Step-eigenValue',region=regions,cf1=-1.0)
        region=regionToolset.Region(nodes=nodes_seabed)
        m.EncastreBC(name="seabed",createStepName='Initial',region=region)


    def pipeSpringEigen_FixedRotation_neutral(self):
        '''
        This function is used for the linear analysis of pipe with fixed rotation at pipeEnds, The connection is
        equal to the connection with neutral line. X direction is the pipe longitudinal direction; Y direction is the
        pipe lateral direction.  The displacement in Z direction (upheaval) is fixed. Simplify-supported Bc, with UR1=0
        '''
        self.pipeSpringEigen_FreeRotation_Bottom()
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        m.boundaryConditions['bottom'].setValues(ur1=SET)
        m.boundaryConditions['top-Simply'].setValues(ur1=SET)

        elastic_center = connectorBehavior.ConnectorElasticity(components=(2, ), coupling=COUPLED,
                                            table=((self.K*self.sPipe, ), ))
        elastic_center.ConnectorOptions()
        m.sections['connectorSec_center'].setValues(behaviorOptions = (elastic_center, ) )

        elastic_ends = connectorBehavior.ConnectorElasticity(components=(2, ), coupling=COUPLED,
                                            table=((self.K*self.sPipe/2, ), ))
        elastic_center.ConnectorOptions()
        m.sections['connectorSec_end'].setValues(behaviorOptions = (elastic_ends, ) )

    def pipeSpringEigen_FixedRotation_Bottom(self):
        '''
        This function is used for the linear analysis of pipe with fixed rotation at pipeEnds, The connection is
        equal to the connection with bottom line. X direction is the pipe longitudinal direction; Y direction is the
        pipe lateral direction.  The displacement in Z direction (upheaval) is fixed. Simplify-supported Bc, with UR1=0
        '''
        self.pipeSpringEigen_FreeRotation_Bottom()
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        m.boundaryConditions['bottom'].setValues(ur1=SET)
        m.boundaryConditions['top-Simply'].setValues(ur1=SET)

    def pipeSpringRiks_FreeRotation_Bottom(self):
        """
        This function is directly descended from the former function: pipeSpringEigen_FreeRotation_Bottom()
        It is used for the Riks simulation of a pipe with fixed buckling length.
        """
        self.pipeSpringEigen_FreeRotation_Bottom()
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        m.StaticRiksStep(name='Step-eigenValue', previous='Initial',
        maxLPF=1.0, maxNumInc=1000, minArcInc=1e-9,initialArcInc=0.001,nlgeom=ON)
        m.steps.changeKey(fromName='Step-eigenValue', toName='Step-riks')

        region=a.sets['wire_ends']
        m.HistoryOutputRequest(name='connector-op',createStepName='Step-riks', variables=('CU1',
         'CU2', 'CU3', 'CUR1', 'CUR2', 'CUR3'), region=region)

        regions=a.sets['n_center']
        m.HistoryOutputRequest(name='n_center',createStepName='Step-riks',variables=('U',),region=regions)
        regions=a.sets['n_side']
        m.HistoryOutputRequest(name='n_side',createStepName='Step-riks',variables=('U',),region=regions)
        regionDef=a.sets['wire_ends']
        m.FieldOutputRequest(name='connector-op', createStepName='Step-riks', variables=('CTF', 'CEF', 'CU', 'CUE',
                'CUP'), region=regionDef)
        m.FieldOutputRequest(name='sectionForce',createStepName='Step-riks',variables=('SF',))
        regions=a.sets['n_top']
        m.ConcentratedForce(name='refLoad',createStepName='Step-riks',region=regions,cf1=-self.Pref)


    def pipeSpringStatic_FixRotation_Bottom(self):
        """
        This function is directly descended from the former function: pipeSpringEigen_FreeRotation_Bottom()
        It is used for the Quasi-static simulation of a pipe with fixed buckling length.
        """
        self.pipeSpringEigen_FreeRotation_Bottom()
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly

        m.StaticStep(name='Step-eigenValue', previous='Initial',maxNumInc=1000,timeIncrementationMethod=FIXED,
        initialInc=0.001,noStop=OFF,nlgeom=ON)
        m.steps.changeKey(fromName='Step-eigenValue', toName='static')

        region=a.sets['wire_ends']
        m.HistoryOutputRequest(name='connector-op',createStepName='static', variables=('CU1',
         'CU2', 'CU3', 'CUR1', 'CUR2', 'CUR3'), region=region)

        regions=a.sets['n_center']
        m.HistoryOutputRequest(name='n_center',createStepName='static',variables=('U',),region=regions)
        regions=a.sets['n_side']
        m.HistoryOutputRequest(name='n_side',createStepName='static',variables=('U',),region=regions)
        regionDef=a.sets['wire_ends']
        m.FieldOutputRequest(name='connector-op', createStepName='static', variables=('CTF', 'CEF', 'CU', 'CUE',
                'CUP'), region=regionDef)
        m.FieldOutputRequest(name='sectionForce',createStepName='static',variables=('SF',))

        # regions=a.sets['n_top']
        # m.loads['unitLoad'].suppress()
        m.boundaryConditions['top-Simply'].move('Initial', 'static')
        m.boundaryConditions['top-Simply'].setValues(u1=-self.L/10)

        m.boundaryConditions['bottom'].setValues(ur1=SET)
        m.boundaryConditions['top-Simply'].setValues(ur1=SET)

    def pipeSpringRiks_FreeRotation_Bottom_symm(self):
        """
        This function is directly descended from the former function: pipeSpringEigen_FreeRotation_Bottom()
        It is used for the Riks simulation of a pipe with fixed buckling length.
        The symmetry is accounted for based on analytical solutions. Hence, only half pipe is considered.
        Tips: The input length should be set to half under this circumstance.
        """
        self.pipeSpringEigen_FreeRotation_Bottom()
        m=mdb.models[self.pipeModel]
        p=m.parts[self.pipePart]
        a=m.rootAssembly
        m.StaticRiksStep(name='Step-eigenValue', previous='Initial',
        maxLPF=1.0, maxNumInc=100, minArcInc=1e-9,initialArcInc=0.001,nlgeom=ON)
        m.steps.changeKey(fromName='Step-eigenValue', toName='Step-riks')

        region=a.sets['wire_ends']
        m.HistoryOutputRequest(name='connector-op',createStepName='Step-riks', variables=('CU1',
         'CU2', 'CU3', 'CUR1', 'CUR2', 'CUR3'), region=region)

        regions=a.sets['n_bottom']
        m.HistoryOutputRequest(name='n_center',createStepName='Step-riks',variables=('U',),region=regions) # due to half-model
        regions=a.sets['n_side']
        m.HistoryOutputRequest(name='n_side',createStepName='Step-riks',variables=('U',),region=regions)
        regionDef=a.sets['wire_ends']
        m.FieldOutputRequest(name='connector-op', createStepName='Step-riks', variables=('CTF', 'CEF', 'CU', 'CUE',
                'CUP'), region=regionDef)
        regions=a.sets['n_top']
        m.ConcentratedForce(name='refLoad',createStepName='Step-riks',region=regions,cf1=-self.Pref)

        m.boundaryConditions['bottom'].suppress()
        regions=a.sets['n_bottom']
        m.XsymmBC(name='center',createStepName='Initial',region=regions)


    def initImperfectionMode1(self):
        '''
        this method is used to produce a real initial imperfection  mode 1 via offset of the pipe nodes
        the shape function is based on buckling modes, as illustrated in script: pipeCalculation.py
        L_b: buckling length; L: overall length of pipe
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]
        num=4.493*2/self.L_b   # nL=8.986

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2    # x0: [-L/2,L/2] for full model.
            # x0=node_coord[0]   # x0: [0,L]   # only for symmetry model
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.
            dy2=-self.w0/15.7*(1+8.9868**2/8-(num*x0)**2/2-np.cos(num*x0)/np.cos(8.9868/2))  # equation from Taylor, 1986
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

    def initImperfectionMode6_sym(self):
        '''
        this method is used to produce a real initial imperfection  mode 6 via offset of the pipe nodes
        Due to symm, only half of the model is modified.
        the shape function is based on buckling modes, as illustrated in script: pipeCalculation.py
        L_b: buckling length; L: overall length of pipe. The imperfection shape is the buckling shape accounting for torsion
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]

        # Pcn,a,b,temp1,temp2=self.find_Pc_freeTorsion(0.00001)
        temp1=(2*self.a*self.b*np.cosh(self.a/2)*np.cos(self.b/2)+(self.a**2-
            self.b**2)*np.sinh(self.a/2)*np.sin(self.b/2))/(2*self.a*self.b*(np.cosh(self.a/2)**2*np.cos(self.b/2)**2+np.sinh(self.a/2)**2*np.sin(self.b/2)**2))
        temp2=((self.a**2-self.b**2)*np.cosh(self.a/2)*np.cos(self.b/2)-
               2*self.a*self.b*np.sinh(self.a/2)*np.sin(self.b/2))/(2*self.a*self.b*(np.cosh(self.a/2)**2*np.cos(self.b/2)**2+np.sinh(self.a/2)**2*np.sin(self.b/2)**2))

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]  # x0: [0,L]
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.
            dy2=self.w0*(1-temp1*np.cosh(self.a*x0/self.L_b/2)*np.cos(self.b*x0/self.L_b/2)+temp2*np.sinh(self.a*x0/self.L_b/2)*np.sin(self.b*x0/self.L_b/2))
            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

    def initImperfectionMode6(self):
        '''
        this method is used to produce a real initial imperfection  mode 6 via offset of the pipe nodes
        this one is used for non_symetry imperfection
        the shape function is based on buckling modes, as illustrated in script: pipeCalculation.py
        L_b: buckling length; L: overall length of pipe. The imperfection shape is the buckling shape accounting for torsion
        '''
        m=mdb.models[self.pipeModel]
        a=m.rootAssembly
        nodes=a.instances[self.pipePart+str(-1)].nodes

        pipeDef=nodes.getByBoundingCylinder((0.5*self.L-0.5*self.L_b-0.2,0,0),(0.5*self.L+0.5*self.L_b+0.2,0,0),0.2) # tolerance 0.2
        defLen=len(pipeDef)
        a.Set(nodes=pipeDef,name='pipeDef')
        nodeDeflib=[]

        # Pcn,a,b,temp1,temp2=self.find_Pc_freeTorsion(0.00001)
        temp1=(2*self.a*self.b*np.cosh(self.a/2)*np.cos(self.b/2)+(self.a**2-
            self.b**2)*np.sinh(self.a/2)*np.sin(self.b/2))/(2*self.a*self.b*(np.cosh(self.a/2)**2*np.cos(self.b/2)**2+np.sinh(self.a/2)**2*np.sin(self.b/2)**2))
        temp2=((self.a**2-self.b**2)*np.cosh(self.a/2)*np.cos(self.b/2)-
               2*self.a*self.b*np.sinh(self.a/2)*np.sin(self.b/2))/(2*self.a*self.b*(np.cosh(self.a/2)**2*np.cos(self.b/2)**2+np.sinh(self.a/2)**2*np.sin(self.b/2)**2))

        for node in pipeDef:
            node_coord=node.coordinates
            label=node.label
            x0=node_coord[0]-self.L/2  # x0: [-L,L]
            y0=node_coord[1]
            z0=node_coord[2]
            #dy2=self.w0*sin(pi*x0/self.L_b) # the default value of L_b is equal to the pipe length L.
            dy2=self.w0*(1-temp1*np.cosh(self.a*x0/self.L_b)*np.cos(self.b*x0/self.L_b)+temp2*np.sinh(self.a*x0/self.L_b)*np.sin(self.b*x0/self.L_b))
            offset_2=dy2
            n_select=nodes[label-1:label]
            a.editNode(nodes=n_select,offset2=offset_2)
            nodeDeflib.append(label)
            defLen=defLen-1
            a.Set(name='deformedNode',nodes=nodes.sequenceFromLabels(nodeDeflib))
            print (defLen,label,dy2)

    def jobEigen_freeRotation_bottom(self):
        import job
        mdb.Job(name='pipeEigen_freeRotation_bottom',model=self.pipeModel,description='Eigenvalue analyses starts')
        #mdb.jobs['pipeEigen'].submit()

    def jobEigen_fixedRotation_bottom(self):
        import job
        mdb.Job(name='pipeEigen_fixedRotation_bottom_'+str(self.pipeModel),model=self.pipeModel,description='Eigenvalue analyses starts')
        #mdb.jobs['pipeEigen'].submit()

    def jobEigen_fixedRotation_neutral(self):
        import job
        mdb.Job(name='pipeEigen_fixedRotation_neutral_'+str(self.pipeModel),model=self.pipeModel,description='Eigenvalue analyses starts')
        #mdb.jobs['pipeEigen'].submit()

    def jobRiks_freeRotation_bottom(self):
        import job
        mdb.Job(name='Riks_'+str(self.pipeModel),model=self.pipeModel,description='Riks starts',numCpus=8,numDomains=8)
        #mdb.jobs['pipeEigen'].submit()

    def jobStatic_fixRotation_bottom(self):
        import job
        mdb.Job(name='QS_'+str(self.pipeModel),model=self.pipeModel,description='Static starts',numCpus=8, numDomains=8)
        #mdb.jobs['pipeEigen'].submit()

def main():
    dataInput=DataInput()
    od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,Pref,sPipe,sSeabed,initialT,finalT,\
    jobName,flagEigen,w0,L_b,P_u,K,mu_L,E,a,b=dataInput.dataInput()
    ## creat an input data based for GUI inputs
    lateralPipeBucklingSpring=LateralBucklingSpring(od,t,L,pipePart,seabedPart,w_bed,pipeModel,pipeMaterial,Pref,sPipe,sSeabed,\
                                                initialT,finalT,jobName,w0,L_b,P_u,K,mu_L,E,a,b)
    lateralPipeBucklingSpring.createPerfectPipeWithPartition()

    # ##############################################################
    # #Combination 1: free rotation ends, bottom, eigenvalue
    # lateralPipeBucklingSpring.pipeSpringEigen_FreeRotation_Bottom()
    # lateralPipeBucklingSpring.jobEigen_freeRotation_bottom()
    # ##############################################################

    # ##############################################################
    # #Combination 2: fixed rotation ends, bottom, eigenvalue
    # lateralPipeBucklingSpring.pipeSpringEigen_FixedRotation_Bottom()
    # lateralPipeBucklingSpring.jobEigen_fixedRotation_bottom()
    # ##############################################################

    # ##############################################################
    # #Combination 3: fixed rotation ends, neural, eigenvalue
    # lateralPipeBucklingSpring.pipeSpringEigen_FixedRotation_neutral()
    # lateralPipeBucklingSpring.jobEigen_fixedRotation_neutral()
    # ##############################################################


    ##############################################################
    #Combination 1: free rotation ends, with linear torsion, bottom, Riks,Mode1
    # tips:  when there is a symmetry, only half model is accounted for.
    # lateralPipeBucklingSpring.pipeSpringStatic_FixRotation_Bottom()

    lateralPipeBucklingSpring.pipeSpringRiks_FreeRotation_Bottom()
    lateralPipeBucklingSpring.initImperfectionMode1()
    # lateralPipeBucklingSpring.initImperfectionMode6()

    # lateralPipeBucklingSpring.pipeSpringRiks_FreeRotation_Bottom_symm()
    # lateralPipeBucklingSpring.initImperfectionMode6_sym()

    lateralPipeBucklingSpring.jobRiks_freeRotation_bottom()
    ##############################################################


if __name__ == '__main__':
        main()