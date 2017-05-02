/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

typedef double (*TEF)(const double, const double);

#pragma pack(push)
#pragma pack(4)
struct generic_header
{
  long long nself;    // Grösse dieser Struktur  
  long long nDatatyp; // Grösse des zu speichernden Datentyps
  long long nself_and_data; 
  long long nDims;
  long long nDimX;
  long long nDimY;
  long long nDimZ;
  int       bAtom;
  int       bComplex;
  double    t;
  double    xMin;
  double    xMax;
  double    yMin;
  double    yMax;
  double    zMin;
  double    zMax;
  double    dx;
  double    dy;
  double    dz;
  double    dkx;
  double    dky;
  double    dkz;
  double    dt;
  int       ks; // Koordinaten-System 
  int       nFuture[100];
  double    dFuture[100];  
};
#pragma pack(pop)
#endif

#if __APPLE__&&__MACH__
  #ifndef APPLE_COMPAT
  #define APPLE_COMPAT
    extern void sincos( double, double*, double* );
    extern void sincosf( float, float*, float* );
  #endif
#endif

