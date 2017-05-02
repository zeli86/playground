try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os,sys,glob,struct,math,subprocess
import vtk

#os.chdir("/home/zeli/tmp/delme1/")
cwd_folder=os.getcwd()

#-------------------------------------------------------------
def window(LT, title):
    RenderView1 = GetRenderView()
    RenderView1.ViewSize = [450, 800]
    RenderView1.CacheKey = 0.0
    RenderView1.Background = [1.0, 1.0, 1.0]  
    RenderView1.UseLight = 1
    RenderView1.LightSwitch = 0
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
    
    RenderView1.CameraFocalPoint = [15, 7.5, 0.0] # for GOST
    RenderView1.CameraPosition = [15, 7.5, -57] # for GOST
    RenderView1.CameraViewAngle = 30 # for GOST
    RenderView1.CenterAxesVisibility = 0


def layout(representation, array, LT, PF):

    DataRepresentation1 = Show()
    
    DataRepresentation1.CubeAxesXTitle = ''
    DataRepresentation1.CubeAxesYTitle = ''
      
    DataRepresentation1.AxesOrigin = [0.0, 0.0, 0.0]
    DataRepresentation1.Origin = [-0.0, 0.0, 0.0]
    DataRepresentation1.EdgeColor = [0.0, 0.0, 0.0]
    #DataRepresentation1.SelectionPointFieldDataArrayName = 'Potential'
    DataRepresentation1.ScalarOpacityFunction = PF
    DataRepresentation1.ScalarOpacityUnitDistance = 0.5802894042140944
    DataRepresentation1.LookupTable = LT
    DataRepresentation1.CubeAxesVisibility = 1
    DataRepresentation1.ScaleFactor = 1.0
    DataRepresentation1.CubeAxesColor = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesXAxisVisibility = 1
    DataRepresentation1.CubeAxesYAxisVisibility = 1
    DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 0
    DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 0
    DataRepresentation1.CubeAxesTickLocation = 'Outside'
    
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
cam1 = GetActiveCamera()
cam1.Roll(90)
#cam1.Azimuth(180)

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
    psirefLUT = MakeBlueToRedLT(my_range[0],my_range[1])
    psirefPWF = GetOpacityTransferFunction('Psi_sol')
    
    clip1 = Clip(Input=final_pvtu)
    clip1.ClipType = 'Plane'
    clip1.Scalars = ['POINTS', 'Psi_sol']
    clip1.Value = -0.7669909596442928
    clip1.InsideOut = 0
    clip1.Crinkleclip = 0
    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [20.0, 15.0, 0.0]
    clip1.ClipType.Normal = [0.0, -1.0, 0.0]
    clip1.ClipType.Offset = 0.0

    Hide(final_pvtu,GetRenderView())

    clip2 = Clip(Input=clip1)
    clip2.ClipType = 'Plane'
    clip2.Scalars = ['POINTS', 'Psi_sol']
    clip2.ClipType.Origin = [30.0, 7.5, 0.0]
    clip2.ClipType.Normal = [-1.0, 0.0, 0.0]
    clip2.ClipType.Offset = 0.0

    Hide(clip1,GetRenderView())

    # show data in view
    SetActiveSource(clip2)
    clip2Display = Show()
    clip2Display.Representation = 'Surface'
    clip2Display.ColorArrayName = ['POINTS', 'Psi_sol']

    window(psirefLUT, '       Psi_sol')
    layout('Surface', 'Psi_sol', psirefLUT, psirefPWF)
    Render()
    WriteImage( current_folder + "sol.png")
    Delete(psirefLUT)
    del psirefLUT
    
    Delete(psirefPWF)
    del psirefPWF
    
    Delete(clip2)
    del clip2
    
    Delete(clip1)
    del clip1

    Delete(final_pvtu)
    del final_pvtu

    #for f in GetSources().values():
    #  print GetSources().values()
    #  if f.GetProperty("Input") is not None:
    #    Delete(f)    
