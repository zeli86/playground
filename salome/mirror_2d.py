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

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/zeli')

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
geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)

h = 20
wh = 5
bp = -2
NumberOfSegments = 100

Vertex_1 = geompy.MakeVertex(-wh, 0, 0)
Vertex_2 = geompy.MakeVertex(wh, 0, 0)
Vertex_3 = geompy.MakeVertex(0, bp, 0)
Vertex_4 = geompy.MakeVertex(-wh, h, 0)
Vertex_5 = geompy.MakeVertex(wh, h, 0)

Arc_1 = geompy.MakeArc(Vertex_1, Vertex_3, Vertex_2)
Arc_1_vertex_2 = geompy.GetSubShape(Arc_1, [2])
Line_1 = geompy.MakeLineTwoPnt(Arc_1_vertex_2, Vertex_4)
Line_1_vertex_3 = geompy.GetSubShape(Line_1, [3])
Line_2 = geompy.MakeLineTwoPnt(Line_1_vertex_3, Vertex_5)
Arc_1_vertex_3 = geompy.GetSubShape(Arc_1, [3])
Line_3 = geompy.MakeLineTwoPnt(Vertex_5, Arc_1_vertex_3)

Face_1 = geompy.MakeFaceWires([Arc_1, Line_1, Line_2, Line_3], 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Arc_1, 'Arc_1' )
geompy.addToStudyInFather( Arc_1, Arc_1_vertex_2, 'Arc_1:vertex_2' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_3, 'Line_1:vertex_3' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudyInFather( Arc_1, Arc_1_vertex_3, 'Arc_1:vertex_3' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Face_1, 'Face_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Face_1)
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(NumberOfSegments)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
