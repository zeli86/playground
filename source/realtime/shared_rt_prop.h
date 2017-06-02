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

  #include "utils_complex.h"

  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0.0;
    
     const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
         tmp1 += fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_vals[q_point][0]+vec_vals[q_point][1]*vec_vals[q_point][1]);
      }
    }

    double retval;
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }  
  
  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
#if DIMENSION==2
    Point<dim,double> pt1( m_xmin, m_ymin );
    Point<dim,double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim,double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim,double> pt2( m_xmax, m_ymax, m_zmax );
#endif

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);

//     unsigned int tmp1[2], tmp2[2];
//     tmp1[0] = triangulation.n_cells();
//     tmp1[1] = triangulation.n_active_cells();
//     MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::setup_system( const bool initial_step )
  {
    m_computing_timer.enter_section(__func__);
    if( initial_step )
    {
      dof_handler.distribute_dofs (fe);

      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

      m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    }

    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_t.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());

    int myrank;
    MPI_Comm_rank( mpi_communicator, &myrank );    
    cout << "(" << myrank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;
    
    system_rhs = 0;

    vector<bool> mask (dof_handler.get_fe().n_components(), true );    
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints, ComponentMask(mask));
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    m_computing_timer.exit_section();
  }  
  
  template<int dim>
  void MySolver<dim>::Expectation_value_position( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	 JxWxn = fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_vals[q_point][0]+vec_vals[q_point][1]*vec_vals[q_point][1]);
	 spacept = fe_values.quadrature_point(q_point);
         tmp[0] += spacept[0]*JxWxn;
         tmp[1] += spacept[1]*JxWxn;
#if dim == 3
         tmp[2] += spacept[2]*JxWxn;
#endif
	}
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template<int dim>
  void MySolver<dim>::Expectation_value_variance( LA::MPI::Vector& vec, double* retval, double* x1 )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	 JxWxn = fe_values.JxW(q_point)*(vec_vals[q_point][0]*vec_vals[q_point][0]+vec_vals[q_point][1]*vec_vals[q_point][1]);
	 spacept = fe_values.quadrature_point(q_point);
         tmp[0] += (spacept[0]-x1[0])*(spacept[0]-x1[0])*JxWxn;
         tmp[1] += (spacept[1]-x1[1])*(spacept[1]-x1[1])*JxWxn;
#if dim == 3
         tmp[2] += (spacept[2]-x1[2])*(spacept[2]-x1[2])*JxWxn;
#endif
	}
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  
  
  template<int dim>
  void MySolver<dim>::Expectation_value_momentum( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxW;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> vec_grads(n_q_points, vector<Tensor<1,dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        fe_values.get_function_gradients( vec, vec_grads );
	
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
	  JxW = fe_values.JxW(q_point);
          tmp[0] += JxW*(vec_vals[q_point][0]*vec_grads[q_point][1][0] - vec_vals[q_point][1]*vec_grads[q_point][0][0]);
	  tmp[1] += JxW*(vec_vals[q_point][0]*vec_grads[q_point][1][1] - vec_vals[q_point][1]*vec_grads[q_point][0][1]);
#if dim == 3
	  tmp[2] += JxW*(vec_vals[q_point][0]*vec_grads[q_point][1][2] - vec_vals[q_point][1]*vec_grads[q_point][0][2]); 
#endif
	}
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    m_computing_timer.enter_section(__func__);
    string filename;
    ComputeIntensity<dim> intensities;
    ComputePhase<dim> phase;
    
    vector<std::string> solution_names;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi , intensities);
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();
    
    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    m_computing_timer.exit_section();
  }   

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi);

    triangulation.save( filename.c_str() );
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system(true);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(system_rhs);

    m_Psi = system_rhs;
  }    

  