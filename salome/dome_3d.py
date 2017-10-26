# -*- coding: utf-8 -*-

# ATUS-PRO - The ATUS-PRO package is atom interferometer Toolbox developed at ZARM
# (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
# founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
# 50WM0942, 50WM1042, 50WM1342.
# Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
#
# This file is part of ATUS-PRO.
#
# ATUS-PRO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ATUS-PRO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ATUS-PRO.  If not, see <http://www.gnu.org/licenses/>.

### IMPORTANT NOTE ###
# 1) export mesh to med file format 
# 2) open med file with gmsh (http://gmsh.info/)
# 3) save mesh in gmsh format
# This makes sure that the physical volumes and surfaces are read in correctly by deal.ii

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

Rmax = 5
dz = 0.2
h = 10

NO_LAYER=int(Rmax/dz)-1

print 'no of layer = %d' % (NO_LAYER) 

arr_disks = []
arr_trans = []
arr_faces = []
arr_solids = []
arr_faces_2 = []
surface_parts = []

Divided_Cylinder = geompy.MakeDividedCylinder(Rmax, h, GEOM.SQUARE)
geompy.addToStudy( Divided_Cylinder, 'Divided Cylinder' )

for i in range(0, NO_LAYER+1):
    arr_disks.append(geompy.MakeDividedDisk(math.sqrt(Rmax**2-(i*dz)**2), 1, GEOM.SQUARE))
    arr_trans.append(geompy.MakeTranslation(arr_disks[-1], 0, 0, -i*dz))
    arr_faces.append(geompy.ExtractShapes(arr_trans[-1], geompy.ShapeType["FACE"], True))
    geompy.addToStudy( arr_disks[-1], 'Divided_Disk_%d' % (i) )
    geompy.addToStudy( arr_trans[-1], 'Translation_Disk_%d' % (i) )
    for j in range(0, 5):
        geompy.addToStudyInFather( arr_disks[-1], arr_faces[i-1][j], 'Face_%d_%d' % (i,j) )  
        
for i in range(0,NO_LAYER):
    for j in range(0,5):
        arr_solids.append(geompy.MakeHexa2Faces(arr_faces[i][j], arr_faces[i+1][j]))
        geompy.addToStudy( arr_solids[-1], 'Hexahedral_Solid_%d_%d' % (i,j) )
                
arr_solids.append( Divided_Cylinder )

domain = geompy.MakeCompound(arr_solids)
geompy.addToStudy( domain, 'domain' )

arr_faces_2.append(geompy.ExtractShapes(Divided_Cylinder, geompy.ShapeType["FACE"], True))
surface_parts.append(arr_faces_2[-1][1])
surface_parts.append(arr_faces_2[-1][2])
surface_parts.append(arr_faces_2[-1][4])
surface_parts.append(arr_faces_2[-1][6])
surface_parts.append(arr_faces_2[-1][11])
surface_parts.append(arr_faces_2[-1][16])
surface_parts.append(arr_faces_2[-1][18])
surface_parts.append(arr_faces_2[-1][19])
surface_parts.append(arr_faces_2[-1][20])

for i in range(0,NO_LAYER):
    for j in range(0,5):
        arr_faces_2.append(geompy.ExtractShapes(arr_solids[j+5*i], geompy.ShapeType["FACE"], True))
        for k in range(0,6):
            geompy.addToStudyInFather( arr_solids[j+5*i], arr_faces_2[-1][k], 'Hexahedral_Face_%d_%d_%d' % (i,j,k) )
            #path = "/Geometry/Hexahedral_Solid_%d_%d/Hexahedral_Face_%d_%d_%d" % (i,j,i,j,k)
            #SubFaceList = geompy.SubShapeAllSortedCentres(theStudy.FindObjectByPath( path ), geompy.ShapeType["FACE"])
            if j == 0 and k == 1 : 
                surface_parts.append(arr_faces_2[-1][k])
            elif j == 1 and k == 1 : 
                surface_parts.append(arr_faces_2[-1][k])
            elif j == 3 and k == 4 : 
                surface_parts.append(arr_faces_2[-1][k])
            elif j == 4 and k == 4 : 
                surface_parts.append(arr_faces_2[-1][k])
            elif i == NO_LAYER-1 and j == 0 and k == 3 :
                surface_parts.append(arr_faces_2[-1][k])
            elif i == NO_LAYER-1 and j == 1 and k == 3 :
                surface_parts.append(arr_faces_2[-1][k])
            elif i == NO_LAYER-1 and j == 2 and k == 2 :
                surface_parts.append(arr_faces_2[-1][k])
            elif i == NO_LAYER-1 and j == 3 and k == 2 :
                surface_parts.append(arr_faces_2[-1][k])
            elif i == NO_LAYER-1 and j == 4 and k == 2 :
                surface_parts.append(arr_faces_2[-1][k])
                
Surface_Group_1 = geompy.CreateGroup(domain, geompy.ShapeType["FACE"])
geompy.UnionList(Surface_Group_1, surface_parts)
geompy.addToStudyInFather( domain, Surface_Group_1, 'Surface_Group_1' )

Vol_Group_1 = geompy.CreateGroup(domain, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Vol_Group_1, geompy.SubShapeAllIDs(domain, geompy.ShapeType["SOLID"]))
geompy.addToStudyInFather( domain, Vol_Group_1, 'Vol_Group_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

from salome.StdMeshers import StdMeshersBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(domain)
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(10)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Quadrangle_Parameters_1 = Quadrangle_2D.QuadrangleParameters(StdMeshersBuilder.QUAD_TRIANGLE_PREF,-1,[],[])
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
isDone = Mesh_1.Compute()
Surface_Group_1_1 = Mesh_1.GroupOnGeom(Surface_Group_1,'Phys_Surface',SMESH.FACE)
phys_vol = Mesh_1.CreateEmptyGroup( SMESH.VOLUME, 'Phys_Volume' )
phys_vol.AddFrom( Mesh_1.GetMesh() )


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Surface_Group_1_1, 'Phys_Surface')
smesh.SetName(phys_vol, 'Phys_Volume')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)

