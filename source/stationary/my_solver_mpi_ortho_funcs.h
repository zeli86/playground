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

/** Želimir Marojević
 */

namespace variant_1
{
  template<int dim>
  double fun( const gsl_vector* x, void *params )
  {
    MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(params); 
    double f1, f2;
    double t[2];
    t[0] = gsl_vector_get (x, 0);
    t[1] = gsl_vector_get (x, 1);

    f1 = t[0]*sol->m_T[0] + t[1]*sol->m_I12 + t[0]*t[0]*t[0]*sol->m_W[0] + 3.0*t[0]*t[1]*t[1]*sol->m_W[2] + 3.0*t[0]*t[0]*t[1]*sol->m_W[3] + t[1]*t[1]*t[1]*sol->m_W[4];
    f2 = t[1]*sol->m_T[1] + t[0]*sol->m_I12 + t[1]*t[1]*t[1]*sol->m_W[1] + 3.0*t[0]*t[1]*t[1]*sol->m_W[4] + 3.0*t[0]*t[0]*t[1]*sol->m_W[2] + t[0]*t[0]*t[0]*sol->m_W[3];

  return f1*f1+f2*f2;
  }

  template<int dim>
  void fun_df( const gsl_vector* x, void *params, gsl_vector* df )
  {
    MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(params); 
    double f1, f2;
    double t[2];
    t[0] = gsl_vector_get (x, 0);
    t[1] = gsl_vector_get (x, 1);

    f1 = t[0]*sol->m_T[0] + t[1]*sol->m_I12 + t[0]*t[0]*t[0]*sol->m_W[0] + 3.0*t[0]*t[1]*t[1]*sol->m_W[2] + 3.0*t[0]*t[0]*t[1]*sol->m_W[3] + t[1]*t[1]*t[1]*sol->m_W[4];
    f2 = t[1]*sol->m_T[1] + t[0]*sol->m_I12 + t[1]*t[1]*t[1]*sol->m_W[1] + 3.0*t[0]*t[1]*t[1]*sol->m_W[4] + 3.0*t[0]*t[0]*t[1]*sol->m_W[2] + t[0]*t[0]*t[0]*sol->m_W[3];

    gsl_vector_set(df, 0, 2.0*f1*(sol->m_T[0]+3*t[0]*t[0]*sol->m_W[0]+3.0*t[1]*t[1]*sol->m_W[2]+6.0*t[0]*t[1]*sol->m_W[3]) + 2.0*f2*(sol->m_I12+6.0*t[0]*t[1]*sol->m_W[2]+3.0*t[0]*t[0]*sol->m_W[3]+3*t[1]*t[1]*sol->m_W[4]) );
    gsl_vector_set(df, 1, 2.0*f1*(sol->m_I12+3.0*t[1]*t[1]*sol->m_W[4]+6.0*t[0]*t[1]*sol->m_W[2]+3*t[0]*t[0]*sol->m_W[3]) + 2.0*f2*(sol->m_T[1]+3*t[1]*t[1]*sol->m_W[1]+6.0*t[0]*t[1]*sol->m_W[4]+3.0*t[0]*t[0]*sol->m_W[2]) );
  }

  template<int dim>
  void fun_fdf(const gsl_vector *x, void *params, double *f, gsl_vector* df) 
  {
    *f = fun<dim>(x, params); 
    fun_df<dim>(x, params, df);
  }
}

namespace variant_2
{
  template<int dim>
  double fun( const gsl_vector* x, void *params )
  {
    MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(params); 
    double t[2];
    t[0] = gsl_vector_get (x, 0);
    t[1] = gsl_vector_get (x, 1);

  return (0.25*sol->m_W[0]*pow(t[0],4) + 0.5*sol->m_T[0]*pow(t[0],2) + 0.25*sol->m_W[1]*pow(t[1],4) + t[0]*sol->m_W[4]*pow(t[1],3) +1.5*sol->m_W[2]*pow(t[0],2)*pow(t[1],2) +0.5*sol->m_T[1]*pow(t[1],2) +pow(t[0],3)*sol->m_W[3]*t[1] +t[0]*sol->m_I12*t[1]);
  }

  template<int dim>
  void fun_df( const gsl_vector* x, void *params, gsl_vector* df )
  {
    MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(params); 
    double f1, f2;
    double t[2];
    t[0] = gsl_vector_get (x, 0);
    t[1] = gsl_vector_get (x, 1);

    f1 = t[0]*sol->m_T[0] + t[1]*sol->m_I12 + t[0]*t[0]*t[0]*sol->m_W[0] + 3.0*t[0]*t[1]*t[1]*sol->m_W[2] + 3.0*t[0]*t[0]*t[1]*sol->m_W[3] + t[1]*t[1]*t[1]*sol->m_W[4];
    f2 = t[1]*sol->m_T[1] + t[0]*sol->m_I12 + t[1]*t[1]*t[1]*sol->m_W[1] + 3.0*t[0]*t[1]*t[1]*sol->m_W[4] + 3.0*t[0]*t[0]*t[1]*sol->m_W[2] + t[0]*t[0]*t[0]*sol->m_W[3];

    gsl_vector_set(df, 0, f1 );
    gsl_vector_set(df, 1, f2 );
  }

  template<int dim>
  void fun_fdf(const gsl_vector *x, void *params, double *f, gsl_vector* df) 
  {
    *f = fun<dim>(x, params); 
    fun_df<dim>(x, params, df);
  }
}
