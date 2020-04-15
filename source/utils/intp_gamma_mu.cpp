//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
//
// This file is part of atus-pro testing.
//
// atus-pro testing is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// atus-pro testing is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
//

/*
Želimir Marojević

aufruf z.B. intp_gamma_mu gamma_von_mu.txt 10 20 30 40 50 100 200

Interpoliert aus den gamma von mu Kurven zu den gegebenen gammas 10 20 usw den zugehörigen Eigenwert.
*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "strtk.hpp"
 
using namespace std; 

typedef vector<double> Row;
vector<Row> table;

int main(int argc, char *argv[])
{
  Row gammas;
  
  if( argc < 3 ) return EXIT_SUCCESS;
  
  for( int i=2; i<argc; ++i )
  {
    double data = atof(argv[i]);
    gammas.push_back( data );
    //cout << data << endl;
  }  

  // Load table from file
  ifstream file( argv[1] );
  if( !file.is_open() )
  {
    cout << "Could not open file: " << argv[1] << "." << endl;
    return EXIT_FAILURE;
  }
  
  vector<std::string> vec;
  
  const char * raute = "#";
  char first;
  
  while(file)
  {
    vec.clear();
    
    string line;
    getline(file, line);
    
    first = line[0];
    
    if( strcmp(raute,&first) == 0 ) continue;

    strtk::parse(line,";",vec);
    
    if( vec.size() == 0 ) continue;

    Row row;
    for( string str : vec )
    {
      double data = stod(str);
      row.push_back(data);
    }
    table.push_back(row);
  }

  const unsigned dim = table.size();
  
  if( dim == 0 )
  {
    cout << "File " << argv[1] << " is empty." << endl;
    return EXIT_FAILURE;
  }
 
  double *x = new double[dim];
  double *y = new double[dim];
  double mu, gamma;  
  
  for( unsigned i=0; i<dim; ++i )
  {
    Row row = table[i];
    
    x[i] = row[1]*row[2];
    y[i] = row[0];   
    
//     cout << x[i] << ", " << y[i] << endl;   
  }

  const double max_gamma = x[dim-1];
  const double min_gamma = x[0];
  const double L1 = 0.5;
  const double L2 = 0.6;
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, dim);
  gsl_spline_init(spline, x, y, dim);
  
  double pos = 0.1;
  for( vector<double>::iterator it = gammas.begin(); it != gammas.end(); it++ )
  {
    gamma = *it;
    if( gamma < max_gamma && gamma > min_gamma )
    {
      mu = gsl_spline_eval(spline, gamma, acc);
      cout << pos << "\t" << mu << "\n";
      cout << 0.5*(2.0*pos+L1) << "\t" << mu << "\n";
      cout << pos+L1 << "\t" << mu << "\n";
      cout << "\n";
    }
    pos+= L2;
  }
  
  // clean up
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  delete [] x;
  delete [] y;
return EXIT_SUCCESS;
}