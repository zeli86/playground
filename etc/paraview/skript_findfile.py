try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


import os,sys,glob,struct,math,subprocess
import numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from xml.dom import minidom

cwd_folder=os.getcwd()

#-------------------------------------------------------------
def window(LT, title):
    RenderView1 = CreateRenderView()
    RenderView1.ViewSize = [700, 600]
    RenderView1.CacheKey = 0.0
    #RenderView1.StereoType = 0.0
    #RenderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
    RenderView1.Background = [1.0, 1.0, 1.0]  
    RenderView1.UseLight = 1
    #RenderView1.StereoRender = 0
    RenderView1.LightSwitch = 0
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.InteractionMode = '2D'
    #RenderView1.StereoCapableWindow = 0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
    
    RenderView1.CameraFocalPoint = [8.951139508746778, -0.5593922620017225, 0.0]
    RenderView1.CenterAxesVisibility = 0
    RenderView1.CameraParallelScale = 14.142135623730947
    RenderView1.CenterOfRotation = [10.0, 0.0, 0.0]
    RenderView1.CameraClippingRange = [54.09460598986377, 55.460631393648214]
    RenderView1.CameraPosition = [8.951139508746778, -0.5593922620017225, 54.64101615137755]
    
    #ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='Psi_ref', TitleFontSize=12, Enabled=1, LabelFontSize=10, LabelColor=[0.0, 0.0, 0.0], TitleColor=[0.0, 0.0, 0.0], Position=[0.9, 0.25], LookupTable=LT )
    #GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

    # get color legend/bar for psirefLUT in view renderView1
    psirefLUTColorBar = GetScalarBar(LT, RenderView1)
    psirefLUTColorBar.Position = [0.85, 0.25]
    psirefLUTColorBar.Position2 = [0.12, 0.42999999999999994]
    psirefLUTColorBar.Title = title
    psirefLUTColorBar.ComponentTitle = ' ' 
    psirefLUTColorBar.TitleColor=[0.0, 0.0, 0.0]
    psirefLUTColorBar.TitleFontSize=12
    psirefLUTColorBar.LabelColor=[0.0, 0.0, 0.0]
    psirefLUTColorBar.LabelFontSize=10    
    GetRenderView().Representations.append(psirefLUTColorBar)

def layout(representation, array, LT, PF):
    DataRepresentation1 = Show()
    
    #DataRepresentation1.CubeAxisTitleFont =  ['Arial', '20', '0', '0', 'Arial', '20', '0', '0', 'Arial', '20', '0', '0', 'Arial', '20', '0', '0']  
    DataRepresentation1.CubeAxesXTitle = 'X'
    DataRepresentation1.CubeAxesYTitle = 'Y'
  
    
    DataRepresentation1.AxesOrigin = [0.0, 0.0, 0.0]
    DataRepresentation1.Origin = [-0.0, 0.0, 0.0]
    DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
    #DataRepresentation1.SelectionPointFieldDataArrayName = 'Potential'
    DataRepresentation1.ScalarOpacityFunction = PF
    DataRepresentation1.ScalarOpacityUnitDistance = 0.5802894042140944
    DataRepresentation1.LookupTable = LT
    DataRepresentation1.ScaleFactor = 2.0
    DataRepresentation1.CubeAxesVisibility = 1
    DataRepresentation1.CubeAxesColor = [0.0, 0.0, 0.0]
    
    #DataRepresentation1.CustomRangeActive = [1, 1, 1]
    #DataRepresentation1.CustomBoundsActive = [1, 1, 1]
    #DataRepresentation1.CustomRange = [xmin, xmax, ymin, ymax, 0.0, 0.0]
    #DataRepresentation1.CustomBounds = [xmin, xmax, ymin, ymax, 0.0, 0.0]
    
    DataRepresentation1.ColorArrayName = ('POINT_DATA', array)
    DataRepresentation1.Representation = representation

#-------------------------------------------------------------------------
#Durchsucht alles unter "pathname" und macht liste aus gefundenen Datein:  
file_list=[]
sp = subprocess.Popen( 'find . -type f -name "final-*.pvtu" | sed \'s:^\./::\'', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
file_list = sp.communicate()[0].split()

#print file_list

for file in file_list:

    print file 
    
    current_folder = cwd_folder + '/' + os.path.dirname( file ) + '/';
    os.chdir( current_folder );

    filename = os.path.basename( file )
    filename_noext = os.path.splitext(filename)[0]

    #-------------------------------------------------------------
    #XML-File-PATH:
    xmldoc = minidom.parse('info.xml')

    mu=xmldoc.getElementsByTagName('MU')[0].firstChild.data
    gs=xmldoc.getElementsByTagName('GS')[0].firstChild.data
    N=xmldoc.getElementsByTagName('N')[0].firstChild.data
    xmin=xmldoc.getElementsByTagName('XMIN')[0].firstChild.data
    xmax=xmldoc.getElementsByTagName('XMAX')[0].firstChild.data
    ymin=xmldoc.getElementsByTagName('YMIN')[0].firstChild.data
    ymax=xmldoc.getElementsByTagName('YMAX')[0].firstChild.data
    zmin=xmldoc.getElementsByTagName('ZMIN')[0].firstChild.data
    zmax=xmldoc.getElementsByTagName('ZMAX')[0].firstChild.data
    npm1=float(xmldoc.getElementsByTagName('NPM1')[0].firstChild.data)
    Smax=float(xmldoc.getElementsByTagName('PSI_REF_MAX')[0].firstChild.data)
    Smin=float(xmldoc.getElementsByTagName('PSI_REF_MIN')[0].firstChild.data)

    #print("mu=\t{}\ngs=\t{}\nN=\t{}" .format(mu, gs, N))
    #print("x=\t{}\t{}\ny=\t{}\t{}\nz=\t{}\t{}" .format(xmin, xmax, ymin, ymax, zmin, zmax))
    #print("Psi_ref=\t{}\t{}" .format(Smin, Smax))

    #-------------------------------------------------------------
    psirefLUT = MakeBlueToRedLT( Smin, Smax )
    psirefPWF = GetOpacityTransferFunction('Psiref')

    #-----------------------------------
    subdomainLUT = MakeBlueToRedLT( 0, npm1 ) 
    subdomainPWF = GetOpacityTransferFunction('subdomain')

    final_pvtu = XMLPartitionedUnstructuredGridReader(PointArrayStatus=['Psi_ref', 'Potential', 'system_rhs', 'subdomain'], FileName=[filename] )

    #-------------------------------------------------------------------
    window(psirefLUT, '       Psi_ref')
    layout('Wireframe', 'Psi_ref', psirefLUT, psirefPWF)
    Render()
    WriteImage( filename_noext +"_psi_ref_WF.png")
    RenderView = GetRenderView()
    Delete(RenderView)

    window(psirefLUT, '       Psi_ref')
    layout('Surface', 'Psi_ref', psirefLUT, psirefPWF)
    Render()
    WriteImage( filename_noext + "_psi_ref_SF.png")
    RenderView = GetRenderView()
    Delete(RenderView)

    #layout('Wireframe', 'subdomain', subdomainLUT, subdomainPWF)
    #Render()
    #WriteImage( filename_noext + "_subdomain_WF.png")
    
    Delete(final_pvtu)

    #for f in GetSources().values():
    #  print GetSources().values()
    #  if f.GetProperty("Input") is not None:
    #    Delete(f)    
