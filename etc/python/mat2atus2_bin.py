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

# DANGER Will Robinson
# This script is supposed to extract the final oct wavefunction and store it in the atus2 bin format.

import os,sys,struct,math

import numpy as np
import scipy.io as sio

noargs = len(sys.argv)

if ( noargs != 2 ):
  print( "No filename specified." )
  exit()

matlabmat = sio.loadmat(sys.argv[1])

#newfilename = os.path.splitext( sys.argv[1] )[0] + ".bin"
newfilename = "octfin.bin"

#data for the binary header
nself = 1380
nDatatyp = 16
nself_and_data = 1380+16*int(matlabmat['nDimX']*matlabmat['nDimY'])
nDims = 2
nDimX = int(matlabmat['nDimX'])
nDimY = int(matlabmat['nDimY'])
nDimZ = 1
bAtom = 1
bComplex = 1
t = 0
xMin = float(matlabmat['xMin'])
xMax = float(matlabmat['xMax'])
yMin = float(matlabmat['yMin'])
yMax = float(matlabmat['yMax'])
zMin = 0
zMax = 0
dx = float(matlabmat['dx'])
dy = float(matlabmat['dy'])
dz = 1
dkx = float(math.pi/(xMin-xMax))
dky = float(math.pi/(yMin-yMax))
dkz = 0
dt = 0.001

#  int       ks; // Koordinaten-System 
#  int       nFuture[100];
#  double    dFuture[100];  

#header_raw = fh.read(struct.calcsize("llllllliidddddddddddddd"))
#header = struct.unpack( "llllllliidddddddddddddd", header_raw )

print("dims = (%ld,%ld,%ld)\n" % (nDimX,nDimY,nDimZ))
print("xrange = (%g,%g)\n" % (xMin,xMax))
print("yrange = (%g,%g)\n" % (yMin,yMax))
print("zrange = (%g,%g)\n" % (zMin,zMax))

try:
    fh = open(newfilename, 'wb+')
    fh.write( struct.pack( "llllllliidddddddddddddd", nself, nDatatyp, nself_and_data, nDims, nDimX, nDimY, nDimZ, bAtom, bComplex, t, xMin, xMax, yMin, yMax, zMin, zMax, dx, dy, dz, dkx, dky, dkz, dt )) 
    
    for i in range(0, 301):
        fh.write( struct.pack("i", 0 ) );
        
    for i in range(0, nDimX*nDimY):
        fh.write( struct.pack("dd", matlabmat['finpsi'][i].real, matlabmat['finpsi'][i].imag ) );        
        

finally:
    fh.close()
