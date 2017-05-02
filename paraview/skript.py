try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


import os,sys,glob,struct,math
import numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from xml.dom import minidom

#file = '/home/merle/PV-z/solution-00024.pvtu'
files = '/home/merle/PV-z/test/final-00022.pvtu'

pfad_fig="/home/merle/PV-z/test/"


#-------------------------------------------------------------
#XML-File-PATH:
xmldoc = minidom.parse('/home/merle/PV-z/test/info.xml')


mu=xmldoc.getElementsByTagName('MU')[0].firstChild.data
gs=xmldoc.getElementsByTagName('GS')[0].firstChild.data
N=xmldoc.getElementsByTagName('N')[0].firstChild.data
xmin=xmldoc.getElementsByTagName('XMIN')[0].firstChild.data
xmax=xmldoc.getElementsByTagName('XMAX')[0].firstChild.data
ymin=xmldoc.getElementsByTagName('YMIN')[0].firstChild.data
ymax=xmldoc.getElementsByTagName('YMAX')[0].firstChild.data
zmin=xmldoc.getElementsByTagName('ZMIN')[0].firstChild.data
zmax=xmldoc.getElementsByTagName('ZMAX')[0].firstChild.data
Smax=float(xmldoc.getElementsByTagName('PSI_REF_MAX')[0].firstChild.data)
Smin=float(xmldoc.getElementsByTagName('PSI_REF_MIN')[0].firstChild.data)

print("mu=\t{}\ngs=\t{}\nN=\t{}" .format(mu, gs, N))
print("x=\t{}\t{}\ny=\t{}\t{}\nz=\t{}\t{}" .format(xmin, xmax, ymin, ymax, zmin, zmax))
print("Psi_ref=\t{}\t{}" .format(Smin, Smax))


#-------------------------------------------------------------


a1_Psi_ref_PiecewiseFunction = CreatePiecewiseFunction( Points=[Smin, 0.0, 0.5, 0.0, Smax, 1.0, 0.5, 0.0] )

a=0.11111100000000002*Smax*2.0-abs(Smin)
b=0.3650795*Smax*2.0-abs(Smin)
c=0.4920634999999999*Smax*2.0-abs(Smin)
d=0.6190474999999999*Smax*2.0-abs(Smin)
e=0.873016*Smax*2.0-abs(Smin)

a1_Psi_ref_PVLookupTable = GetLookupTableForArray( "Psi_ref", 1, NanColor=[0.498039, 0.0, 0.0], RGBPoints=[Smin, 0.0, 0.0, 0.562493, a, 0.0, 0.0, 1.0, b, 0.0, 1.0, 1.0, c, 0.500008, 1.0, 0.500008, d, 1.0, 1.0, 0.0, e, 1.0, 0.0, 0.0, Smax, 0.500008, 0.0, 0.0], ScalarOpacityFunction=a1_Psi_ref_PiecewiseFunction, VectorMode='Magnitude', ScalarRangeInitialized=1.0 )


#a1_system_rhs_PiecewiseFunction = CreatePiecewiseFunction( Points=[-3.6654039377026493e-06, 0.0, 0.5, 0.0, 3.6654039377026493e-06, 1.0, 0.5, 0.0] )

#a1_Potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 70.0, 1.0, 0.5, 0.0] )
#a1_Potential_PVLookupTable = GetLookupTableForArray( "Potential", 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 35.0, 0.865, 0.865, 0.865, 70.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ScalarOpacityFunction=a1_Potential_PiecewiseFunction, ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

a1_subdomain_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0] )
a1_subdomain_PVLookupTable = GetLookupTableForArray( "subdomain", 1, RGBPoints=[0.0, 0.231373, 0.298039, 0.752941, 1.5, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902], VectorMode='Magnitude', NanColor=[0.247059, 0.0, 0.0], ScalarOpacityFunction=a1_subdomain_PiecewiseFunction, ColorSpace='Diverging', ScalarRangeInitialized=1.0 )


#a1_system_rhs_PiecewiseFunction = CreatePiecewiseFunction( Points=[-3.6654039377026493e-06, 0.0, 0.5, 0.0, 3.6654039377026493e-06, 1.0, 0.5, 0.0] )
#a1_system_rhs_PVLookupTable = GetLookupTableForArray( "system_rhs", 1, RGBPoints=[-3.6654039377026493e-06, 0.23, 0.299, 0.754, 0.0, 0.865, 0.865, 0.865, 3.6654039377026493e-06, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ScalarOpacityFunction=a1_system_rhs_PiecewiseFunction, ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

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


def window(LT):
    RenderView1 = CreateRenderView()
    RenderView1.ViewSize = [700, 600]
    RenderView1.CacheKey = 0.0
    RenderView1.StereoType = 0
    #RenderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
    RenderView1.Background = [1.0, 1.0, 1.0]  
    RenderView1.UseLight = 1
    RenderView1.StereoRender = 0
    RenderView1.LightSwitch = 0
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.InteractionMode = '2D'
    RenderView1.StereoCapableWindow = 0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
    
    RenderView1.CameraFocalPoint = [8.951139508746778, -0.5593922620017225, 0.0]
    RenderView1.CenterAxesVisibility = 0
    RenderView1.CameraParallelScale = 14.142135623730947
    RenderView1.CenterOfRotation = [10.0, 0.0, 0.0]
    RenderView1.CameraClippingRange = [54.09460598986377, 55.460631393648214]
    RenderView1.CameraPosition = [8.951139508746778, -0.5593922620017225, 54.64101615137755]
    
    ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='Psi_ref', TitleFontSize=12, Enabled=1, LabelFontSize=10, LabelColor=[0.0, 0.0, 0.0], TitleColor=[0.0, 0.0, 0.0], Position=[0.9, 0.25], LookupTable=LT )
    GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)


final00022_pvtu = XMLPartitionedUnstructuredGridReader(PointArrayStatus=['Psi_ref', 'Potential', 'system_rhs', 'subdomain'], FileName=[files] )
#-------------------------------------------------------------------
#window(a1_subdomain_PVLookupTable)
window(a1_Psi_ref_PVLookupTable)
layout('Wireframe', 'Psi_ref', a1_Psi_ref_PVLookupTable, a1_Psi_ref_PiecewiseFunction)
Render()
WriteImage(pfad_fig + "psi_ref_WF.png")
#pic1='pic_WF.png'

layout('Surface', 'Psi_ref', a1_Psi_ref_PVLookupTable, a1_Psi_ref_PiecewiseFunction)
Render()
WriteImage(pfad_fig + "psi_ref_SF.png")
#pic2='pic_SF.png'


window(a1_subdomain_PVLookupTable)

layout('Surface', 'subdomain', a1_subdomain_PVLookupTable, a1_subdomain_PiecewiseFunction)
Render()
WriteImage(pfad_fig + "subdomain_SF.png")
#pic3='pic_SF.png'

layout('Wireframe', 'subdomain', a1_subdomain_PVLookupTable, a1_subdomain_PiecewiseFunction)
Render()
WriteImage(pfad_fig + "subdomain_WF.png")
#pic4='pic_SF.png'

#---Bilder-zusammen----------------------------------------------------------------
#f, axarr = plt.subplots(2, 2)
#f.text(.5,.95, "title", ha='center')

#image_3d=(pfad_fig + pic1)
#image3d=plt.imread(image_3d)
#axarr[0, 0].imshow(image3d)
#axarr[0, 0].set_title('...')
#axarr[0, 0].axis('off')

#image_x=(pfad_fig + pic2)
#imagex=plt.imread(image_x)
#axarr[0, 1].imshow(imagex)
#axarr[0, 1].set_title('...')
#axarr[0, 1].axis('off')

#image_y=(pfad_fig +pic3)
#imagey=plt.imread(image_y)
#axarr[1, 0].imshow(imagey)
#axarr[1, 0].set_title('...')
#axarr[1, 0].axis('off')

#image_z=(pfad_fig + pic4)
#imagez=plt.imread(image_z)
#axarr[1, 1].imshow(imagez)
#axarr[1, 1].set_title('...')
#axarr[1, 1].axis('off')


#plt.savefig(pfad_fig + "all_in_one.png", dpi=400)
