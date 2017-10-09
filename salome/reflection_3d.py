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
h = 100

NO_LAYER=int(Rmax/dz)-1

arr_disks = []
arr_trans = []
arr_faces = []
arr_solids = []

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
status = Mesh_1.RemoveHypothesis(Hexa_3D)
isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
