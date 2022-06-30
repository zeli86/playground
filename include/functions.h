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


#pragma once

#include <deal.II/base/function.h>
#include "eigenfunctions_AiryAi.h"
#include "eigenfunctions_HO.h"
#include "eigenfunctions_Polar_HO.h"

using namespace std;
using namespace dealii;

template <int dim>
class CEigenfunctions : public Function<dim>
{
public:
  CEigenfunctions(unsigned QN[3], std::vector<double> a) : Function<dim>()
  {
    m_QNx = QN[0];
    m_QNy = QN[1];
    m_QNz = QN[2];

#if POTENTIAL==1
    m_fakx = a[0];
#endif
#if POTENTIAL==2
    m_fakx = pow(a[0], 2.0 / 3.0);
#endif

    m_faky = a[1];
    m_fakz = a[2];
  }

  CEigenfunctions(const std::vector<int>& QN, std::vector<double> a) : Function<dim>()
  {
    m_QNx = QN[0];
    m_QNy = QN[1];
    m_QNz = QN[2];

#if POTENTIAL==1
    m_fakx = a[0];
#endif
#if POTENTIAL==2
    m_fakx = pow(a[0], 2.0 / 3.0);
#endif

    m_faky = a[1];
    m_fakz = a[2];
  }
  virtual double value(const Point<dim>& p, const unsigned component = 0) const;
  virtual Tensor<1, dim> gradient(const Point<2>& p, const unsigned component = 0) const;

private:
  unsigned int m_QNx;
  unsigned int m_QNy;
  unsigned int m_QNz;
  double m_fakx;
  double m_faky;
  double m_fakz;
};

template <int dim>
double CEigenfunctions<dim>::value(const Point<dim>& p, const unsigned) const
{
  double retval;

  switch (dim)
  {
  case 1:
#if POTENTIAL==1
    retval = (*EF_HO[m_QNx])(m_fakx, p(0));
#endif
#if POTENTIAL==2
    retval = (*EF_AIRY[m_QNx])(m_fakx, p(0));
#endif
    break;
  case 2:
#if POTENTIAL==1
    retval = (*EF_HO[m_QNx])(m_fakx, p(0)) * (*EF_HO[m_QNy])(m_faky, p(1));
#endif
#if POTENTIAL==2
    retval = (*EF_AIRY[m_QNx])(m_fakx, p(0)) * (*EF_HO[m_QNy])(m_faky, p(1));
#endif
    break;
  case 3:
#if POTENTIAL==1
    retval = (*EF_HO[m_QNx])(m_fakx, p(0)) * (*EF_HO[m_QNy])(m_faky, p(1)) * (*EF_HO[m_QNz])(m_fakz, p(2));
#endif
#if POTENTIAL==2
    retval = (*EF_AIRY[m_QNx])(m_fakx, p(0)) * (*EF_HO[m_QNy])(m_faky, p(1)) * (*EF_HO[m_QNz])(m_fakz, p(2));
#endif
    break;
  }
  return retval;
}

template <int dim>
Tensor<1, dim> CEigenfunctions<dim>::gradient(const Point<2>& p, const unsigned) const
{
  Point<dim> retval;
  switch (dim)
  {
  case 1:
#if POTENTIAL==1
    retval[0] = (*EF_HO_DERIV[m_QNx])(m_fakx, p(0));
#endif
#if POTENTIAL==2
    retval[0] = (*EF_AIRY_DERIV[m_QNx])(m_fakx, p(0));
#endif
    break;
  case 2:
#if POTENTIAL==1
    retval[0] = (*EF_HO_DERIV[m_QNx])(m_fakx, p(0));
    retval[1] = (*EF_HO_DERIV[m_QNy])(m_faky, p(1));
#endif
#if POTENTIAL==2
    retval[0] = (*EF_AIRY_DERIV[m_QNx])(m_fakx, p(0));
    retval[1] = (*EF_HO_DERIV[m_QNy])(m_faky, p(1));
#endif
    break;
  case 3:
#if POTENTIAL==1
    retval[0] = (*EF_HO_DERIV[m_QNx])(m_fakx, p(0));
    retval[1] = (*EF_HO_DERIV[m_QNy])(m_faky, p(1));
    retval[2] = (*EF_HO_DERIV[m_QNz])(m_fakz, p(2));
#endif
#if POTENTIAL==2
    retval[0] = (*EF_AIRY_DERIV[m_QNx])(m_fakx, p(0));
    retval[1] = (*EF_HO_DERIV[m_QNy])(m_faky, p(1));
    retval[2] = (*EF_HO_DERIV[m_QNz])(m_fakz, p(2));
#endif
    break;
  }
  return retval;
}

/*************************************************************************************************/

template <int dim>
class CPotential : public Function<dim>
{
public:
  explicit CPotential(std::vector<double> a) : Function<dim>()
  {
    m_fakx = a[0];
    m_faky = a[1];
    m_fakz = a[2];
  }
  virtual double value(const Point<dim>& p, const unsigned component = 0) const;

  double m_fakx;
  double m_faky;
  double m_fakz;
};

/*************************************************************************************************/
template <int dim>
double CPotential<dim>::value(const Point<dim>& p, const unsigned) const
{
  double retval;

  switch (dim)
  {
  case 1:
#if POTENTIAL==1
    retval = m_fakx * m_fakx * p(0) * p(0);
#endif
#if POTENTIAL==2
    retval = m_fakx * p(0);
#endif
    break;
  case 2:
#if POTENTIAL==1
    retval = m_fakx * m_fakx * p(0) * p(0) + m_faky * m_faky * p(1) * p(1);
#endif
#if POTENTIAL==2
    retval = m_fakx * p(0) + m_faky * m_faky * p(1) * p(1);
#endif
    break;
  case 3:
#if POTENTIAL==1
    retval = m_fakx * m_fakx * p(0) * p(0) + m_faky * m_faky * p(1) * p(1) + m_fakz * m_fakz * p(2) * p(2);
#endif
#if POTENTIAL==2
    retval = m_fakx * p(0) + m_faky * m_faky * p(1) * p(1) + m_fakz * m_fakz * p(2) * p(2);
#endif
    break;
  }
  return retval;
}
