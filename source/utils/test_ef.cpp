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

/** Želimir Marojević
 */

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdio>

#include "functions_cs.h"
#include "gsl/gsl_sf.h"

class CEigenfunctionsAi : public Function<1>
{
  public:
    CEigenfunctionsAi ( unsigned  QN, const double a ) : Function<1>() 
    { 
      m_QN=QN; 
      m_fak=pow(a,1.0/3.0);
    }
    virtual double value ( const Point<1> &p, const unsigned int  component = 0) const;
    virtual Tensor<1,1> gradient (const Point<1> &p, const unsigned int component=0) const;
  private:
    unsigned m_QN;
    double m_fak;
};

double CEigenfunctionsAi::value( const Point<1> &p, const unsigned int component ) const
{
  double retval = (*EF_AIRY[m_QN])(m_fak,p(0));
return retval;
}

Tensor<1,1> CEigenfunctionsAi::gradient( const Point<1> &p, const unsigned int component ) const
{
  Point<1> retval;
  retval[0] = (*EF_AIRY_DERIV[m_QN])(m_fak,p(0));
return retval;
}

class CEigenfunctionsPHO : public Function<1>
{
  public:
    CEigenfunctionsPHO ( unsigned QN, unsigned s, const double a ) : Function<1>() 
    { 
      m_QN=QN; 
      m_s=s;
      m_fak=a;
    }
    virtual double value ( const Point<1> &p, const unsigned int  component = 0) const;
    virtual Tensor<1,1> gradient (const Point<1> &p, const unsigned int component=0) const;
  private:
    unsigned m_QN;
    unsigned m_s;
    double m_fak;
};

double CEigenfunctionsPHO::value( const Point<1> &p, const unsigned int component ) const
{
  double retval = (*EF_PHO[m_QN+51*m_s])(m_fak,p(0));
return retval;
}

Tensor<1,1> CEigenfunctionsPHO::gradient( const Point<1> &p, const unsigned int component ) const
{
  Point<1> retval;
  retval[0] = (*EF_PHO_DERIV[m_QN+51*m_s])(m_fak,p(0));
return retval;
}

class CPotentialGrav : public Function<1>
{
  public:
    CPotentialGrav ( const double a  ) : Function<1>() 
    { 
      m_fak=a; 
    }
    virtual double value ( const Point<1> &p, const unsigned int  component = 0) const;

    double m_fak;
};
  
/*************************************************************************************************/
double CPotentialGrav::value( const Point<1> &p, const unsigned int component ) const
{
return m_fak*p(0);
}

class CPotentialPHO : public Function<1>
{
  public:
    CPotentialPHO ( const double a, const long s  ) : Function<1>() { m_fak=a*a; m_m = double(s*s); }
    virtual double value ( const Point<1> &p, const unsigned int  component = 0) const;

    double m_fak;
    double m_m;
};

/*************************************************************************************************/
double CPotentialPHO::value( const Point<1> &p, const unsigned int component ) const
{
  double retval=0;
  double rq = p(0)*p(0);
  if( m_m == 0 )
  {
    retval = m_fak*rq;
  }
  else
  {
    if( rq > 0 ) 
      retval = m_fak*rq + m_m/rq;
    else
      retval = 0;
  }
return retval;
}

int main ( int argc, char *argv[] )
{
  const double s = 0;
  const double a = 0;
  const double b = 20;
  const unsigned N = 2000;
  const double dx = (b-a)/double(N-1);

  std::cout << "dx == " << dx << std::endl;

  unsigned A = std::stoi(argv[1]);
  unsigned B = std::stoi(argv[2]);
  unsigned C = std::stoi(argv[3]);
  unsigned D = std::stoi(argv[4]);
  unsigned E = std::stoi(argv[5]);
  
  unsigned QN1[] = {A,B,E};
  unsigned QN2[] = {C,D,E};  
  vector<double> omega = {0.5,0.5,0};
  //vector<double> omega = {1,1,0};


  std::ofstream out( "vals.txt" );

  double retval = 0;
  double retval2 = 0;
  double retval3 = 0;
  double retval4 = 0;
  double retval5 = 0;
  double retval6 = 0;

 #pragma omp parallel 
 {
    CEigenfunctionsAi eigai( A, omega[0]);
    CEigenfunctionsAi eigai2( C, omega[0]);
    CPotentialGrav potg( omega[0]  );

    CEigenfunctionsPHO eigpho( B, E, omega[1] );
    CEigenfunctionsPHO eigpho2( D, E, omega[1] );
    CPotentialPHO potho( omega[1], E );  
    Point<1> ptr;
    Point<1> ptz;

    CEigenfunctions Ef1( QN1, omega );
    CEigenfunctions Ef2( QN2, omega );    
    CPotential pot( omega, C );

    Point<2> pt;

    #pragma omp for reduction(+:retval,retval2)
    for( unsigned i=0; i<N; ++i )
    {
      pt[0] = i*dx;
      ptz[0] = i*dx;
      for( unsigned j=0; j<N; ++j )
      {
        pt[1] = j*dx;
        ptr[0] = j*dx;
 
        double val = Ef1.value(pt) * Ef2.value(pt);
        double grad_val = Ef1.gradient(pt) * Ef2.gradient(pt);
      
        //out << r << "\t" << Ef1.value(pt) << "\t" << Ef2.value(pt) << std::endl;

        //retval += pt[1]*(grad_val + pot.value(pt) * val);
        
        retval += pt[1] * ( eigai.gradient(ptz)*eigai2.gradient(ptz) + eigpho.gradient(ptr) *  eigpho2.gradient(ptr) + (potg.value(ptz) + potho.value(ptr)) * eigai.value(ptz)* eigai2.value(ptz) * eigpho.value(ptr) * eigpho2.value(ptr) );
        retval2 += pt[1]*(eigai.value(ptz)* eigai2.value(ptz) * eigpho.value(ptr) * eigpho2.value(ptr));
      }
    }
  }
  

  CEigenfunctionsAi eigai( A, omega[0]);
  CEigenfunctionsAi eigai2( C, omega[0]);
  CPotentialGrav potg( omega[0]  );

  CEigenfunctionsPHO eigpho( B, E, omega[1] );
  CEigenfunctionsPHO eigpho2( D, E, omega[1] );
  CPotentialPHO potho( omega[1], E );  

  Point<1> pt1;
  for( unsigned i=0; i<N; ++i )
  {
    pt1[0] = i*dx;

    double val1 = eigai.value(pt1) * eigai2.value(pt1);
    double grad_val1 = eigai.gradient(pt1) * eigai2.gradient(pt1);
    double val2 = eigpho.value(pt1) * eigpho2.value(pt1);
    double grad_val2 = eigpho.gradient(pt1) * eigpho2.gradient(pt1);
      
    out << pt1[0] << "\t" << eigai.value(pt1) << "\t" << eigpho.value(pt1) << std::endl;

    retval3 += (grad_val1 + potg.value(pt1) * val1);
    retval4 += val1;
    retval5 += pt1[0]*(grad_val2 + potho.value(pt1) * val2);
    retval6 += pt1[0]*val2;
  }   

  //retval *= dx*dx;
  retval *= dx*dx;
  retval2 *= dx*dx;
  retval3 *= dx;
  retval4 *= dx;
  retval5 *= dx;
  retval6 *= dx;

  //double Energy = pow(omega[0],2.0/3.0) * fabs(gsl_sf_airy_zero_Ai(A+1)) + omega[1] * ( E + 2*B + 1 );

  std::cout << "QN1 == " << QN1[0] << ", " << QN1[1] << ", " << QN1[2] << std::endl;
  std::cout << "QN2 == " << QN2[0] << ", " << QN2[1] << ", " << QN2[2] << std::endl;

  std::cout << "<Psi|H|Psi> = " << retval << std::endl;
  std::cout << "<Psi|Psi> = " << retval2 << std::endl;
  std::cout << "E == " << retval3 << ", N ==" << retval4 << ", " << pow(omega[0],2.0/3.0) * fabs(gsl_sf_airy_zero_Ai(A+1)) << std::endl;
  std::cout << "E == " << retval5 << ", N ==" << retval6 << ", " << 2 * omega[1] * ( E + 2*B + 1 ) << std::endl;
  
return EXIT_SUCCESS;
}
