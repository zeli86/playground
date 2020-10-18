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
    TimerOutput::Scope timing_section(m_computing_timer, "");
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

//     unsigned int tmp1[2], tmp2[2];
//     tmp1[0] = triangulation.n_cells();
//     tmp1[1] = triangulation.n_active_cells();
//     MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    
  }
  
  template <int dim>
  void MySolver<dim>::setup_system( const bool initial_step )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
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
    
  }  

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
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

    
  }   

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_Psi);

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

  