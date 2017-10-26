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

Delta = 0.25
HalfThickness = 2*Delta

NZ1 = 200
NZ2 = 0.5*NZ1-HalfThickness
NY1 = 100    
NX1 = 50
NX2 = 50

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

box_one = geompy.MakeBoxDXDYDZ( NX1*Delta, NY1*Delta, NZ1*Delta )
box_two = geompy.MakeBoxDXDYDZ( NX2*Delta, NY1*Delta, NZ2*Delta )
box_three = geompy.MakeBoxDXDYDZ( NX2*Delta, NY1*Delta, NZ2*Delta )

translation_box_two = geompy.MakeTranslation(box_two, NX1*Delta, 0, 0 )
translation_box_three = geompy.MakeTranslation(box_three, NX1*Delta, 0, (NZ2+2*HalfThickness)*Delta )

domain = geompy.MakeCompound([box_one, translation_box_two, translation_box_three])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( box_one, 'box_one' )
geompy.addToStudy( box_two, 'box_two' )
geompy.addToStudy( box_three, 'box_three' )
geompy.addToStudy( translation_box_two, 'translation_box_two' )
geompy.addToStudy( translation_box_three, 'translation_box_three' )
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


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
