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


#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS

#include <deal.II/base/function.h>
#include "eigenfunctions_AiryAi.h"
#include "eigenfunctions_Polar_HO.h"

using namespace std;
using namespace dealii;
   
class CEigenfunctions : public Function<2>
{
  public:
    CEigenfunctions ( unsigned int QN[3], std::vector<double> a ) : Function<2>() 
    { 
      m_QNx=QN[0]; m_QNy=QN[1]; m_QNz=QN[2]; 
      m_fakx=pow(a[0],1.0/3.0); m_faky=a[1]*a[1]; 
    }
    virtual double value ( const Point<2> &p, const unsigned int  component = 0) const;
    virtual Tensor<1,2> gradient (const Point<2> &p, const unsigned int component=0) const;

  private:
    unsigned int m_QNx;
    unsigned int m_QNy;
    unsigned int m_QNz;
    double m_fakx;
    double m_faky;
};

double CEigenfunctions::value( const Point<2> &p, const unsigned int ) const
{
  double retval = (*EF_AIRY[m_QNx])(m_fakx,p(0)) * (*EF_PHO[m_QNy+50*m_QNz])(m_faky,p(1));
return retval;
}

Tensor<1,2> CEigenfunctions::gradient (const Point<2> &p, const unsigned int ) const
{
  Point<2> retval;
  retval[0] = (*EF_AIRY_DERIV[m_QNx])(m_fakx,p(0));
  retval[1] = (*EF_PHO_DERIV[m_QNy+50*m_QNz])(m_faky,p(1));
return retval;
}

/*************************************************************************************************/

class CPotential : public Function<2>
{
  public:
    CPotential ( const std::vector<double> a, const long l  ) : Function<2>() { m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2]; m_m = double(l*l); }
    virtual double value ( const Point<2> &p, const unsigned int  component = 0) const;

    double m_fakx;
    double m_faky;
    double m_fakz;
    double m_m;
};
  
/*************************************************************************************************/
double CPotential::value( const Point<2> &p, const unsigned int ) const
{
  double retval=0;
  double rq = p(1)*p(1);
  if( m_m == 0 )
  {
    retval = m_fakx*p(0) + m_faky*rq;
  }
  else
  {
    if( rq > 0 ) 
      retval = m_fakx*p(0) + m_faky*rq + m_m/rq;
    else
      retval = 0;
  }
return retval;
}
#endif