try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os,sys,glob,struct,math,subprocess
import vtk

#os.chdir("/home/zeli/tmp/delme1/")
cwd_folder=os.getcwd()

#-------------------------------------------------------------
def window(LT, title):
    #RenderView1 = CreateRenderView()
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
    
    #text1 = Text()
    #text1.Text = 'r'
    #SetActiveSource(text1)
    
    #text1Display = Show(text1, RenderView1)
    #text1Display.Interactivity = 1
    #text1Display.Color = [0, 0, 0]
    #text1Display.Opacity = 1.0
    #text1Display.FontFamily = 'Arial'
    #text1Display.Bold = 0
    #text1Display.Italic = 0
    #text1Display.Shadow = 0
    #text1Display.FontSize = 20
    #text1Display.Justification = 'Left'
    #text1Display.WindowLocation = 'AnyLocation'
    #text1Display.Position = [0.5, 0.02]    

    #text2 = Text()
    #text2.Text = 'z'
    #SetActiveSource(text2)
    
    #text2Display = Show(text2, RenderView1)
    #text2Display.Interactivity = 1
    #text2Display.Color = [0, 0, 0]
    #text2Display.Opacity = 1.0
    #text2Display.FontFamily = 'Arial'
    #text2Display.Bold = 0
    #text2Display.Italic = 0
    #text2Display.Shadow = 0
    #text2Display.FontSize = 20
    #text2Display.Justification = 'Left'
    #text2Display.WindowLocation = 'AnyLocation'
    #text2Display.Position = [0.04, 0.5]    

    #ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='Psi_ref', TitleFontSize=12, Enabled=1, LabelFontSize=10, LabelColor=[0.0, 0.0, 0.0], TitleColor=[0.0, 0.0, 0.0], Position=[0.9, 0.25], LookupTable=LT )
    #GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

    # get color legend/bar for psirefLUT in view renderView1
    ##psirefLUTColorBar = GetScalarBar(LT, RenderView1)
    ##psirefLUTColorBar.Position = [0.85, 0.25]
    ##psirefLUTColorBar.Position2 = [0.12, 0.42999999999999994]
    ##psirefLUTColorBar.Title = title
    ##psirefLUTColorBar.ComponentTitle = ' ' 
    ##psirefLUTColorBar.TitleColor=[0.0, 0.0, 0.0]
    ##psirefLUTColorBar.TitleFontSize=12
    ##psirefLUTColorBar.LabelColor=[0.0, 0.0, 0.0]
    ##psirefLUTColorBar.LabelFontSize=10    
    ##GetRenderView().Representations.append(psirefLUTColorBar)

def layout(representation, array, LT, PF):

    SetActiveSource(final_pvtu)
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
sp = subprocess.Popen( 'find . -type f -name "final-*.vtu" | sed \'s:^\./::\'', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
    final_pvtu = XMLUnstructuredGridReader(PointArrayStatus=['Psi_sol'], FileName=[filename] )
    pd = final_pvtu.PointData

    #for n in range(pd.GetNumberOfArrays()):
    #  print "Range for array ", pd.GetArray(n).GetName(), " ", pd.GetArray(n).GetRange()
    
    my_range = pd.GetArray("Psi_sol").GetRange()
    psirefLUT = MakeBlueToRedLT(my_range[0],my_range[1])
    psirefPWF = GetOpacityTransferFunction('Psi_sol')

    window(psirefLUT, '       Psi_sol')
    layout('Surface', 'Psi_sol', psirefLUT, psirefPWF)
    Render()
    WriteImage( current_folder + "sol.png")
    #Delete(FindSource("Text1"))
    #Delete(FindSource("Text2"))
    
    Delete(psirefLUT)
    del psirefLUT
    
    Delete(psirefPWF)
    del psirefPWF

    Delete(final_pvtu)
    del final_pvtu

    #for f in GetSources().values():
    #  print GetSources().values()
    #  if f.GetProperty("Input") is not None:
    #    Delete(f)    
