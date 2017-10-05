# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

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

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')

isDone = Mesh_1.Compute()


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
