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

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::Expectation_value_position( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxWxn = fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
          spacept = fe_values.quadrature_point(qp);
          tmp[0] += spacept[0]*JxWxn;
          tmp[1] += spacept[1]*JxWxn;
#if dim==3
          tmp[2] += spacept[2]*JxWxn;
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  
  
  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::Expectation_value_momentum( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxW;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> vec_grads(n_q_points, vector<Tensor<1,dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        fe_values.get_function_gradients( vec, vec_grads );

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          tmp[0] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][0] - vec_vals[qp][1]*vec_grads[qp][0][0]);
          tmp[1] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][1] - vec_vals[qp][1]*vec_grads[qp][0][1]);
#if dim == 3
          tmp[2] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][2] - vec_vals[qp][1]*vec_grads[qp][0][2]); 
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template <int dim, int no_time_steps, int no_lam>
  double MySolver<dim,no_time_steps,no_lam>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0.0;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for( unsigned qp=0; qp<n_q_points; qp++ )
         tmp1 += fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
      }
    }

    double retval;
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }    