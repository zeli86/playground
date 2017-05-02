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

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::Expectation_value_position( Vector<double>& vec, double* retval )
  {
    double JxWxn;
    
    *retval=0;
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<1> spacept;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
	 JxWxn = fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_vals[q_point][0]+vec_vals[q_point][1]*vec_vals[q_point][1]);
	 spacept = fe_values.quadrature_point(q_point);
         *retval += spacept[0]*JxWxn;
      }
    }
  }  

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::Expectation_value_momentum( Vector<double>& vec, double* retval )
  {
    *retval=0;
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> vec_grads(n_q_points, vector<Tensor<1,1>>(2));

    typename DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      fe_values.get_function_gradients( vec, vec_grads );
	
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        *retval += fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_grads[q_point][1][0] - vec_vals[q_point][1]*vec_grads[q_point][0][0]);
      }
    }
  }  

  template <int no_time_steps, int no_lam>
  double MySolver<no_time_steps,no_lam>::Particle_Number( Vector<double>& vec )
  {
    double retval=0;
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        retval += fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_vals[q_point][0]+vec_vals[q_point][1]*vec_vals[q_point][1]);
    }

  return retval;
  }    