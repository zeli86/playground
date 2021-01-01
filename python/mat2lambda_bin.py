#!/usr/bin/env python3

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

# DANGER will Robinson
# This script is supposed to extract the final control parameter from the matlab file 
# and store it in binary and txt format.

import os,sys,struct 
import numpy as np
import scipy.io as sio

noargs = len(sys.argv)

if ( noargs != 2 ):
  print( "No filename specified." )
  exit()

matlabmat = sio.loadmat(sys.argv[1])

lambdas = matlabmat['finlam']
dims = np.shape(lambdas)
NT = 201.0
T = 2.0
dt = T/(NT-1)

print("dims = (%ld,%ld)\n" % (dims[0],dims[1]))
print("dt = (%g)\n" % (dt))

with open("lambda_fin.bin", "wb") as binout:
  binout.write(struct.pack("ii", *bytearray([dims[1],dims[0]])))

  for s in range(0, dims[1]):
    for i in range(0, dims[0]):
        binout.write(struct.pack("d", lambdas[i][s] ))

with open("lambda_fin.txt", "w") as txtout:
  for i in range(0, dims[0]):
    txtout.write( "%g\t" % (i*dt) )
    for s in range(0, dims[1]):
      txtout.write( "%g\t" % (lambdas[i][s]) )
    txtout.write( "\n" )
