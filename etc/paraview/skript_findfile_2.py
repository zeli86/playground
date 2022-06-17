try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


import os,sys,glob,struct,math,subprocess
import numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk

cwd_folder=os.getcwd()

#-------------------------------------------------------------
def window(LT, title):
    RenderView1 = CreateRenderView()
    RenderView1.ViewSize = [800, 600]
    RenderView1.CacheKey = 0.0
    RenderView1.Background = [1.0, 1.0, 1.0]  
    RenderView1.UseLight = 1
    RenderView1.LightSwitch = 0
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
    
    RenderView1.CameraFocalPoint = [15, 0, 0.0] # for GOST
    RenderView1.CameraPosition = [15, 0, 80] # for GOST
    RenderView1.CenterAxesVisibility = 0
    #RenderView1.CameraParallelScale = 14.142135623730947
    #RenderView1.CenterOfRotation = [0.0, 0.0, 0.0]
    #RenderView1.CameraClippingRange = [55, 55]
    
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
    final_pvtu = XMLPartitionedUnstructuredGridReader(PointArrayStatus=['Psi_sol','error per cell'], FileName=[filename] )
    pd = final_pvtu.PointData

    #for n in range(pd.GetNumberOfArrays()):
    #  print "Range for array ", pd.GetArray(n).GetName(), " ", pd.GetArray(n).GetRange()
    
    my_range = pd.GetArray("Psi_sol").GetRange()
    my_range_error = pd.GetArray("error per cell").GetRange()
    my_range_subdomains = pd.GetArray("subdomain").GetRange()
    
    psirefLUT = MakeBlueToRedLT(my_range[0],my_range[1])
    errorLUT = MakeBlueToRedLT(my_range_error[0],my_range_error[1])
    subdomainLUT = MakeBlueToRedLT(my_range_subdomains[0],my_range_subdomains[1])
    
    psirefPWF = GetOpacityTransferFunction('Psi_sol')
    errorPWF = GetOpacityTransferFunction('error per cell')
    subdomainPWF = GetOpacityTransferFunction('subdomain')
    
    #-------------------------------------------------------------------
    window(errorLUT, '       error per cell')
    layout('Wireframe', 'error per cell', errorLUT, errorPWF)
    Render()
    WriteImage( current_folder +"error.png")
    Delete(GetRenderView())

    window(psirefLUT, '       Psi_sol')
    layout('Surface', 'Psi_sol', psirefLUT, psirefPWF)
    Render()
    WriteImage( current_folder + "sol.png")
    Delete(GetRenderView())

    window(subdomainLUT, '       subdomain')
    layout('Wireframe', 'subdomain', subdomainLUT, subdomainPWF)
    Render()
    WriteImage( current_folder + "subdomain.png")
    Delete(GetRenderView())
    
    Delete(final_pvtu)

    #for f in GetSources().values():
    #  print GetSources().values()
    #  if f.GetProperty("Input") is not None:
    #    Delete(f)    
