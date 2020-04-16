/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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
    
  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);

    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);

    m_computing_timer.exit_section();
  }
  
  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

    constraints.distribute(vec);
    m_workspace_1=vec;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
         tmp1 += fe_values.JxW(q_point)*fabs(fe_values.quadrature_point(q_point)[1])*(vec_vals[q_point][0]*vec_vals[q_point][0] + vec_vals[q_point][1]*vec_vals[q_point][1]);
      }
    }

    tmp1 *= 2*M_PI;
    double retval;
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }

  template <int dim>
  void MySolver<dim>::setup_system( const bool initial_step )
  {
    m_computing_timer.enter_section(__func__);
    if( initial_step )
    {
      dof_handler.distribute_dofs (fe);
      //DoFRenumbering::component_wise (dof_handler);
      
      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

      m_Psi.reinit (locally_owned_dofs, mpi_communicator);
    }
    
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    newton_update.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_t.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_0.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());


    //cout << "(" << m_rank << ") n_active_cells = " << triangulation.n_active_cells()  << endl;
    pcout << "triangulation.n_global_active_cells() = " << triangulation.n_global_active_cells()  << endl;
    cout << "(" << m_rank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(2), constraints);
    
    if( m_QN1[2] > 0 )
       VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    m_computing_timer.exit_section();
  }
  
  template<int dim>
  void MySolver<dim>::setup_boundary_ids()
  {
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for( unsigned f=0; f<GeometryInfo<dim>::faces_per_cell; f++ )
      {
        const Point<dim> face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary() && !(face_center[1]==0) )    
        {
          cell->face(f)->set_all_boundary_ids(1);
        }
      }
    }
  }
  
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    std::vector<const LA::MPI::Vector*> x_system (2);
    
    m_workspace_1=m_Psi;    
    m_workspace_2=m_Psi_0;    
    x_system[0] = &m_workspace_1;
    x_system[1] = &m_workspace_2;
    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(x_system);

    triangulation.save( filename.c_str() );
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_boundary_ids();

    setup_system(true);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi);

    constraints.distribute(m_Psi);

//     ofstream out ("grid.svg");
//     GridOut grid_out;
//     grid_out.write_svg (triangulation, out);
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    m_computing_timer.enter_section(__func__);
    string filename;
    ComputeIntensity<dim> intensities( "Intensity Psi" ) ;
    ComputeIntensity<dim> intensities2( "Intensity Psi_0" ) ;
    //ComputePhase<dim> phase;
    
    vector<std::string> solution_names;

//     Vector<float> subdomain (triangulation.n_active_cells());
//     for (unsigned int i=0; i<subdomain.size(); ++i)
//       subdomain(i) = triangulation.locally_owned_subdomain();

//     m_Psi.update_ghost_values();
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi , intensities);
    solution_names.clear();
    solution_names.push_back ("Re Psi_0");
    solution_names.push_back ("Im Psi_0");    
    data_out.add_data_vector (m_Psi_0, solution_names);
    data_out.add_data_vector (m_Psi_0, intensities2);
    
    //data_out.add_data_vector (subdomain, "subdomain");
    //data_out.add_data_vector (m_perturbation, "perturbation");
    data_out.build_patches ();
    
    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );
    
    static vector<pair<double,string> > times_and_names;
    times_and_names.push_back (pair<double,string> (m_t, filename));
    ofstream pvd_output ("solution.pvd");
    data_out.write_pvd_record (pvd_output, times_and_names);    

    m_computing_timer.exit_section();
  } 
  
  template<int dim>
  void MySolver<dim>::Expectation_value_position( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    constraints.distribute(vec);
    m_workspace_1=vec;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );
        for( unsigned qp=0; qp<n_q_points; qp++ )
      	{
      	  JxWxn = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1])*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
	        spacept = fe_values.quadrature_point(qp);
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
  void MySolver<dim>::Expectation_value_momentum( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxW;
    
    constraints.distribute(vec);
    m_workspace_1=vec;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> vec_grads(n_q_points, vector<Tensor<1,dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );
        fe_values.get_function_gradients( m_workspace_1, vec_grads );
	
        for( unsigned qp=0; qp<n_q_points; qp++ )
	      {
	        JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
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
