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

Delta = 0.1

slit_Lx = 5*Delta
slit_Ly = 7*Delta
slit_distance = 2*Delta # even number
dom_Lx = 100*Delta
dom_Ly = 100*Delta

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

slit_1 = geompy.MakeFaceHW(slit_Lx, slit_Ly, 1)
slit_1_translation = geompy.MakeTranslation(slit_1, 0.5*(slit_Lx+slit_distance), 0, 0)
slit_2 = geompy.MakeFaceHW(slit_Lx, slit_Ly, 1)
slit_2_translation = geompy.MakeTranslation(slit_2, -0.5*(slit_Lx+slit_distance), 0, 0)

dom_1 = geompy.MakeFaceHW(dom_Lx, dom_Ly, 1)
dom_1_translation = geompy.MakeTranslation(dom_1, 0, 0.5*(slit_Ly+dom_Ly), 0)
dom_2 = geompy.MakeFaceHW(dom_Lx, dom_Ly, 1)
dom_2_translation = geompy.MakeTranslation(dom_2, 0, -0.5*(slit_Ly+dom_Ly), 0)

domain = geompy.MakeCompound([slit_1_translation, slit_2_translation, dom_1_translation, dom_2_translation])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( slit_1, 'slit_1' )
geompy.addToStudy( slit_1_translation, 'slit_1_translation' )
geompy.addToStudy( slit_2, 'slit_2' )
geompy.addToStudy( slit_2_translation, 'slit_2_translation' )
geompy.addToStudy( dom_1, 'dom_1' )
geompy.addToStudy( dom_1_translation, 'dom_1_translation' )
geompy.addToStudy( dom_2, 'dom_2' )
geompy.addToStudy( dom_2_translation, 'dom_2_translation' )
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
isDone = Mesh_1.Compute()

Mesh_1.MergeNodes(Mesh_1.FindCoincidentNodesOnPart( Mesh_1, 1e-5, [], 0 ),[])

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.EDGE,SMESH.FT_FreeBorders,SMESH.FT_Undefined,0)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
Group_1 = Mesh_1.MakeGroupByFilter( 'Phys_Boundary', aFilter )
Group_2 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'Phys_Domain' )
Group_2.AddFrom( Mesh_1.GetMesh() )

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Group_1, 'Phys_Boundary')
smesh.SetName(Group_2, 'Phys_Domain')

isDone = Mesh_1.Compute()

#phys_boundary = Mesh_1.CreateEmptyGroup( SMESH.EDGE, 'Phys_Boundary' )
#phys_boundary.AddFrom( Mesh_1.GetFreeBorders() )

## Criterion : Free Edges
#aBorders = Mesh_1.GetFreeBorders() 
## create groups
#aGroupF = Mesh_1.CreateEmptyGroup(SMESH.FACE, "Faces with free edges")
#aGroupN = Mesh_1.CreateEmptyGroup(SMESH.NODE, "Nodes on free edges")
 
## fill groups with elements, corresponding to the criterion
#print ""
#print "Criterion: Free edges Nb = ", len(aBorders)
#for i in range(len(aBorders)):
   #aBorder = aBorders[i]
   #print "Face # ", aBorder.myElemId, " : Edge between nodes (",
   #print aBorder.myPnt1, ", ", aBorder.myPnt2, ")"
   #aGroupF.Add([aBorder.myElemId])
   #aGroupN.Add([aBorder.myPnt1, aBorder.myPnt2])
   
#smesh.SetName(aGroupF, 'aGroupF')
#smesh.SetName(aGroupN, 'aGroupN')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
