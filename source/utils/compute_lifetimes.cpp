//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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
*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <locale>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
//#include "strtk.hpp"
#include "pugixml.hpp"
 
using namespace std; 

typedef vector<double> Row;
typedef vector<Row> Table;

bool operator<( const Row& rhs, const Row& lhs )
{
return (rhs[0] < lhs[0]);
}

bool operator==( const Row& rhs, const Row& lhs )
{
return( fabs(rhs[0]-lhs[0]) < 1e-5 );
}

bool ReadFile( const std::string& filename, Table& table )
{
  table.clear();

  ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    cout << "Could not open file: " << filename << "." << endl;
    return false;
  }

  string line;
  vector<string> vec;

  while(file)
  {
    line.clear();
    vec.clear();
    getline(file, line);
    strtk::parse(line,"\t",vec);
    if( line.size() == 0 ) continue;

    Row row;
    double data;
    for( auto i : vec )
    {
      try
      { 
        data = stod(i);
      }
      catch( const std::invalid_argument& ia )
      {
        continue;
      }      
      row.push_back(data);
    }

    if( fabs(row[1]-1) < 1e-4 )
    {
      continue;
    } 

    table.push_back(row);
    if( row[2] > 10 ) break;
  }

  // removes duplicate rows
  std::sort( table.begin(), table.end() );
  table.erase( std::unique(table.begin(), table.end()), table.end() );  

  for( auto & i : table )
  {
//    i[2] = log(2*M_PI*i[2]);
    i[2] = log(i[2]);
  }

  if( table.size() < 500 )
    cout << "WARNING: " << filename << endl;
  return true;
}

void do_fit( const Table& table, Row& retval )
{
  const int N = 500;
  const int shift = 0;
  const int n = table.size();

  if( n < N ) return;

  double * t = new double[N];
  double * y = new double[N];

  for( int i=0; i<N; ++i )
  {
    t[i] = table[i+shift][0];
    y[i] = table[i+shift][2];
  }  

  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_linear( t, 1, y, 1, N, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

  // fit ist ax+b
  retval.push_back( log(2.0)/c1 );
  retval.push_back( c1 ); // a 
  retval.push_back( c0 ); // b
  retval.push_back( cov00 );
  retval.push_back( cov01 );
  retval.push_back( cov11 );
  retval.push_back( chisq );

/*
  printf ("# best fit: Y = %g + %g X\n", c0, c1);
  printf ("# covariance matrix:\n");
  printf ("# [ %g, %g\n#   %g, %g]\n", cov00, cov01, cov01, cov11);
  printf ("# chisq = %g\n", chisq);
*/

  delete [] t;
  delete [] y;
}

void output_table( const string& filename, const Table& table )
{
  ofstream out( filename );
  //out.setf(ios::fixed, ios::floatfield);
  //out.precision(10);

  vector<ios::fmtflags> fmt;
  fmt.push_back(ios::fixed);
  fmt.push_back(ios::fixed);
  fmt.push_back(ios::fixed);
  fmt.push_back(ios::scientific);
  fmt.push_back(ios::scientific);
  fmt.push_back(ios::scientific);
  fmt.push_back(ios::scientific);
  fmt.push_back(ios::scientific);

  //out << "# N\ln2/c0\tc0\tc1\tcov00\tcov01\tcov11\tchisq" << endl;

  for( auto row : table )
  {
    out.flags(fmt[0]); 
    out << row[0]; 
    for( int i=1; i<row.size(); ++i )
    {
      //out << "\t" << row[i];
      out.flags(fmt[i]); 
      out << "," << row[i]; 
    }
    out << endl;
  }
}

double my_stod (std::string const& s) {
    std::istringstream iss (s);
    iss.imbue (std::locale("en_US.UTF-8"));
    double d;
    iss >> d;
    // insert error checking.
    return d;
}

int main(int argc, char *argv[])
{
  char tmp[10];
  double part_no;

  vector<string> qn_id;
  qn_id.push_back("0_1_0");
  qn_id.push_back("0_1_1");
  qn_id.push_back("1_0_0");
  qn_id.push_back("1_0_1");
  qn_id.push_back("2_0_0");
  qn_id.push_back("2_0_1");

  Table data, results; 
  for( auto qn_str : qn_id )
  {
    results.clear();
    for( int j=10; j<95; j+=5 )
    {
      sprintf( tmp, "%04d", j  );
      string xmlfile = qn_str + "/" + tmp + "/" + "info.xml"; 
      string errfile = qn_str + "/" + tmp + "/" + "err_prop.txt"; 
            
      pugi::xml_document doc;
      if (!doc.load_file(xmlfile.c_str())) throw;
        
      std::string tmpstr = doc.child( "INFO" ).child( "N" ).child_value();
      part_no = my_stod( tmpstr );

      Row new_row;
      new_row.push_back(part_no);
      ReadFile( errfile, data );
      do_fit( data, new_row );
      if( new_row.size() <= 1 ) continue;
      results.push_back( new_row );
    }  
    output_table( qn_str + "_fits.csv", results );
  }
return EXIT_SUCCESS;
}