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

#import salome_notebook
#notebook = salome_notebook.NoteBook(theStudy)
#sys.path.insert( 0, r'/home/zeli/github/playground/salome')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

Delta = 0.2

slit_height = 6*Delta # even number
slit_width = 5*Delta
slit_depth = 10*Delta
slit_gap = 2*Delta # even number

box_height = 40*Delta
box_width = 40*Delta
box_depth = 40*Delta

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

slit_one = geompy.MakeBoxDXDYDZ(slit_depth, slit_width, slit_height)
translation_slit_one = geompy.MakeTranslation(slit_one, -0.5*slit_depth, -(slit_gap+slit_width), -0.5*slit_height)

slit_two = geompy.MakeBoxDXDYDZ( slit_depth, slit_width, slit_height)
translation_slit_two = geompy.MakeTranslation(slit_two, -0.5*slit_depth, slit_gap, -0.5*slit_height)

box_one = geompy.MakeBoxDXDYDZ(box_width, box_depth, box_height)
translation_box_one = geompy.MakeTranslation(box_one, 0.5*slit_depth, -0.5*box_depth, -0.5*box_height)

box_two = geompy.MakeBoxDXDYDZ(box_width, box_depth, box_height)
translation_box_two = geompy.MakeTranslation(box_two, -0.5*slit_depth-box_width, -0.5*box_depth, -0.5*box_height)

domain = geompy.MakeCompound([translation_slit_one, translation_slit_two, translation_box_one, translation_box_two])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( slit_one, 'slit_one' )
geompy.addToStudy( translation_slit_one, 'translation_slit_one' )
geompy.addToStudy( slit_two, 'slit_two' )
geompy.addToStudy( translation_slit_two, 'translation_slit_two' )
geompy.addToStudy( box_one, 'box_one' )
geompy.addToStudy( translation_box_one, 'translation_box_one' )
geompy.addToStudy( box_two, 'box_two' )
geompy.addToStudy( translation_box_two, 'translation_box_two' )
geompy.addToStudy( domain, 'domain' )

#for k in range(0,6):
#    f_ind_1 = geompy.GetSubShapeID(translation_box_one, box_one_faces[k])
#    print f_ind_1

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

from salome.StdMeshers import StdMeshersBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(domain)
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(Delta,None,1e-11)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Quadrangle_Parameters_1 = Quadrangle_2D.QuadrangleParameters(StdMeshersBuilder.QUAD_QUADRANGLE_PREF,-1,[],[])
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
isDone = Mesh_1.Compute()

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_FreeFaces,SMESH.FT_Undefined,0,SMESH.FT_Undefined,SMESH.FT_Undefined,6.9528e-310)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
Group_phys_surface = Mesh_1.MakeGroupByFilter( 'phys_surface', aFilter )

phys_vol = Mesh_1.CreateEmptyGroup( SMESH.VOLUME, 'Phys_Volume' )
phys_vol.AddFrom( Mesh_1.GetMesh() )

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(phys_vol, 'Phys_Volume')
#smesh.SetName(Group_1, 'Group_1')

isDone = Mesh_1.Compute()


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
