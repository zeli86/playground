# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

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

gap = 2*Delta # even number
dom1_Lx = 50*Delta
dom1_Ly = 200*Delta
dom2_Lx = dom1_Lx
dom2_Ly = 0.5*(dom1_Ly-gap)
dom3_Lx = dom2_Lx
dom3_Ly = dom2_Ly

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

dom_1 = geompy.MakeFaceHW(dom1_Lx, dom1_Ly, 1)
dom_1_translation = geompy.MakeTranslation(dom_1, -0.5*dom1_Lx, 0, 0)
dom_2 = geompy.MakeFaceHW(dom2_Lx, dom2_Ly, 1)
dom_2_translation = geompy.MakeTranslation(dom_2, 0.5*dom2_Lx, 0.5*(dom2_Ly+gap), 0)
dom_3 = geompy.MakeFaceHW(dom2_Lx, dom2_Ly, 1)
dom_3_translation = geompy.MakeTranslation(dom_2, 0.5*dom3_Lx, -0.5*(dom2_Ly+gap), 0)

domain = geompy.MakeCompound([dom_1_translation, dom_2_translation, dom_3_translation])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( dom_1, 'dom_1' )
geompy.addToStudy( dom_1_translation, 'dom_1_translation' )
geompy.addToStudy( dom_2, 'dom_2' )
geompy.addToStudy( dom_2_translation, 'dom_2_translation' )
geompy.addToStudy( dom_3, 'dom_3' )
geompy.addToStudy( dom_3_translation, 'dom_3_translation' )
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

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')

isDone = Mesh_1.Compute()

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
