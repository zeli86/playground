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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>

#include "mpi.h"
#include "functions.h"
//#include "SinCos.h"
#include "ParameterReader.h"
#include "my_table.h"
#include "ref_pt_list.h"
#include "MyRealTools.h"
#include "MyComplexTools.h"

namespace BreedSolver
{
  template <int dim>
  class MySolver;

  using namespace std;
  using namespace dealii;

  #include "my_solver_mpi_ortho_funcs.h"

#ifdef __variant_1__
  using namespace variant_1;
#endif
#ifdef __variant_2__
  using namespace variant_2;
#endif

  #include "CBase.h"

  template <int dim>
  class MySolver : public CBase<2>
  {
  public:
    MySolver( ParameterHandler & );
    virtual ~MySolver();

    void run ();

    double m_T[2];
    double m_W[5];
    double m_I12; 
  protected:
    int DoIter_std_newton( string="" );
    int DoIter( string="" );

    void save( string );
    void make_grid_custom();
    void setup_system( const bool );
    void do_superposition();
    void Expectation_value_position( LA::MPI::Vector&, double* );

    void solve();
    void compute_contributions();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );

    void output_results ( string, string = "step" );
    void output_guess ();

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector newton_update;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_ng;
    LA::MPI::Vector m_Psi_i;
    LA::MPI::Vector m_Psi_f;    
    Vector<double> m_error_per_cell;

    string m_Potential_str;
    
    MyTable m_table;
    MyTable m_results;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  template <int dim>
  MySolver<dim>::MySolver ( ParameterHandler &ph ) 
    : 
    CBase<2>(ph),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (1),
    dof_handler (triangulation)
  {
    m_Potential_str = "0.5*x^2+0.5*y^2";
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    m_computing_timer.enter_section(__func__);

    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );
    
#if DIMENSION==2
    Point<dim,double> pt1( m_xmin, m_ymin );
    Point<dim,double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim,double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim,double> pt2( m_xmax, m_ymax, m_zmax );
#endif
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    
#if DIMENSION==2
    triangulation.refine_global(5);    
    double isovalues[] = {34,32,30,28,26};
#endif
#if DIMENSION==3
    triangulation.refine_global(4);
    double isovalues[] = {16,15,14,13};
#endif
    
    for( unsigned int step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
          if( Potential.value(p)  < isovalues[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }

    unsigned int tmp1[2], tmp2[2];
    tmp1[0] = triangulation.n_cells();
    tmp1[1] = triangulation.n_active_cells();

    MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(vec);
    m_workspace_1=vec;
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<Tensor<1, dim> > vec_grad(n_q_points);

    double JxW, vec_val_q;
    double T1=0.0, N1=0.0, W1=0.0;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vec_vals);
        fe_values.get_function_gradients(m_workspace_1, vec_grad);
        for (unsigned int q_point = 0; q_point<n_q_points; ++q_point)
        {
          JxW = fe_values.JxW(q_point);
          vec_val_q = vec_vals[q_point]*vec_vals[q_point];
          T1 += JxW*( vec_grad[q_point]*vec_grad[q_point] + Potential.value(fe_values.quadrature_point(q_point))*vec_val_q );
          N1 += JxW*vec_val_q;
          W1 += JxW*vec_val_q*vec_val_q;
        }
      }
    }
    MPI_Allreduce( &T1, &T, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &N1, &N, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &W1, &W, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    pcout << "Solving..." << endl;
    
    LA::MPI::Vector tmp_vec(locally_owned_dofs, mpi_communicator);
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, tmp_vec, system_rhs);

    constraints.distribute (tmp_vec);
    newton_update = tmp_vec;

    m_computing_timer.exit_section();
  }  
  
  template<int dim>
  void MySolver<dim>::do_superposition()
  {
    m_workspace_ng = 0;
    m_workspace_ng.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    constraints.distribute (m_workspace_ng);
    m_Psi_ref=m_workspace_ng;
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);
    m_workspace_1=m_Psi_1;
    m_workspace_2=m_Psi_2;
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );

    const QGauss<dim> quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim> > Psi_1_grad(n_q_points);
    vector<Tensor<1, dim> > Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, p12, p1q, p2q, Q;
    double I12 = 0.0;
    double T[2] = {};
    double W[5] = {};

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, Psi_1 );
        fe_values.get_function_values( m_workspace_2, Psi_2 );
        fe_values.get_function_gradients(m_workspace_1, Psi_1_grad);
        fe_values.get_function_gradients(m_workspace_2, Psi_2_grad);

        for (unsigned int qpt = 0; qpt < n_q_points; ++qpt)
        {
          JxW = fe_values.JxW(qpt);
          p12 = Psi_1[qpt]*Psi_2[qpt];
          p1q = Psi_1[qpt]*Psi_1[qpt];
          p2q = Psi_2[qpt]*Psi_2[qpt];
	        Q = Potential.value(fe_values.quadrature_point(qpt)) - m_mu[0]; 

          T[0] += JxW*(Psi_1_grad[qpt]*Psi_1_grad[qpt] + Q*p1q);
          T[1] += JxW*(Psi_2_grad[qpt]*Psi_2_grad[qpt] + Q*p2q);
          I12  += JxW*(Psi_1_grad[qpt]*Psi_2_grad[qpt] + Q*p12);
          W[0] += JxW*p1q*p1q;
          W[1] += JxW*p2q*p2q;
          W[2] += JxW*p1q*p2q;
          W[3] += JxW*p1q*p12;
          W[4] += JxW*p2q*p12;
        }  
      }
    }
    
    for( int i=0; i<5; i++ ) W[i] *= m_gs[0];

    MPI_Allreduce( T, m_T, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &I12, &m_I12, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( W, &m_W, 5, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
    /*
    if( m_root )
    {
      printf( "m_T[0] := %.15e;\n", m_T[0] );
      printf( "m_T[1] := %.15e;\n", m_T[1] );
      printf( "m_I12 := %.15e;\n", m_I12 );
      for( int i=0; i<5; i++ )
        printf( "m_W[%d] := %.15e;\n", i, m_W[i] );
    }
    */
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

      m_Psi_1.reinit (locally_owned_dofs, mpi_communicator);
      m_Psi_2.reinit (locally_owned_dofs, mpi_communicator);
    }

    m_Psi_i.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_f.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_ref.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    newton_update.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    system_rhs = 0;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    CompressedSimpleSparsityPattern csp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);    
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    m_computing_timer.enter_section(__func__);
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_1, "m_Psi_1");
    data_out.add_data_vector (m_Psi_2, "m_Psi_2");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu", mpi_communicator );

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    m_computing_timer.enter_section(__func__);
    
    string filename;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_i, "Psi_i");
    data_out.add_data_vector (m_Psi_f, "Psi_f");
    data_out.build_patches ();

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;

    do_superposition();
    assemble_rhs();
    m_res_old[0]=m_res[0];
    m_counter=0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      m_Psi_2.add( -m_t[1]/fabs(m_t[1]), newton_update); 

#ifdef __variant_1__
      find_ortho_min_2();
#endif
#ifdef __variant_2__
      find_ortho_min_2(false);
#endif

      do_superposition();
      assemble_rhs();
      
      m_resp[0] = m_res_old[0]-m_res[0];
      m_res_old[0] = m_res[0];
      m_res_over_resp[0] = fabs( m_res[0]/m_resp[0] );
	
      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res[0] );
      m_table.insert( cols, MyTable::RESP, m_resp[0] );
      m_table.insert( cols, MyTable::RES_OVER_RESP, m_res_over_resp[0] );
      m_table.insert( cols, MyTable::MU, m_mu[0] );
      m_table.insert( cols, MyTable::GS, m_gs[0] );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if( m_root ) cout << m_table;
      if( m_res[0] < m_epsilon ) { retval=Status::SUCCESS; break; }
      //if( m_res_over_resp[0] > 1e4 ) { retval=Status::SLOW_CONV; break; }
      if( l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

      m_counter++;
    }while( true );

    do_superposition();
    m_N[0] = Particle_Number(m_Psi_ref);
    
    if( m_N[0] < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);
    
    return retval;
  }

  template <int dim>
  int MySolver<dim>::DoIter_std_newton ( string path )
  {
    m_computing_timer.enter_section(__func__);    
    int retval = Status::SUCCESS;

    do_superposition();
    do
    {
      pcout << "--------------------------------------------------------------------------------\n";
      pcout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_ref.add( -1.0, newton_update);  
      constraints.distribute (m_Psi_ref);

      m_resp[0] = m_res_old[0]-m_res[0];
      m_res_over_resp[0] = fabs(m_res[0]/m_resp[0]);

      if( m_counter % m_NA == 0 ) output_results( path );

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res[0] );
      m_table.insert( cols, MyTable::RESP, m_resp[0] );
      m_table.insert( cols, MyTable::RES_OVER_RESP, m_res_over_resp[0] );
      m_table.insert( cols, MyTable::MU, m_mu[0] );
      m_table.insert( cols, MyTable::GS, m_gs[0] );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );      
      m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if( m_root ) cout << m_table;
      if( m_res[0] < m_epsilon ) { retval = Status::SUCCESS; break; }

      m_res_old[0] = m_res[0];
      m_counter++;
    }while( true );

    if( retval == Status::SUCCESS )
    {
      output_results( path, "final" );
    }

    if( m_root )
    {
      string filename = path + "//log.csv";
      m_table.dump_2_file(filename);
    }
    m_computing_timer.exit_section();
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run ()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    //make_grid();
    make_grid_custom();
    setup_system(true);

    m_Potential_str = "0.5*(x+3)^2+0.5*y^2";
    m_guess_str = "exp(-0.5*((x+3)^2+y^2))";
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> guess_fct;
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), m_guess_str , constants );

    VectorTools::interpolate (dof_handler, guess_fct, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
    
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu[0] = T/N+m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu[0] << endl;

    // inital state 
    status = DoIter();
    pcout << "L2_norm of m_Psi_ref: " << m_N[0] << endl;
    
    m_Psi_i = m_Psi_ref;

    // final state
    
    m_Potential_str = "0.5*(x-3)^2+0.5*y^2";
    m_guess_str = "exp(-0.5*((x-3)^2+y^2))";
    
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), m_guess_str , constants );

    VectorTools::interpolate (dof_handler, guess_fct, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;

    status = DoIter();
    pcout << "L2_norm of m_Psi_ref: " << m_N[0] << endl;
    
    m_Psi_f = m_Psi_ref;

    if( status == Status::SUCCESS )
    {
      output_results("","final");
      save("final.bin");
      dump_info_xml();
    }

    if( m_root )
    {
      ofstream ofs("log.txt");
      ofs << m_table;
    }
  }

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
//     LA::MPI::Vector tmp1, tmp2, tmp3;
//     tmp1.reinit(locally_owned_dofs, mpi_communicator);
//     tmp2.reinit(locally_owned_dofs, mpi_communicator);
//     tmp3.reinit(locally_owned_dofs, mpi_communicator);
//     tmp1=m_Psi_i;
//     tmp2=m_Psi_f;
//     tmp3=0;
    
    std::vector<const LA::MPI::Vector*> x_system (4);
    x_system[0] = &m_Psi_i;
    x_system[1] = &m_Psi_i;

    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(x_system);

    triangulation.save( filename.c_str() );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  using namespace BreedSolver;
  deallog.depth_console (0);

return EXIT_SUCCESS;
}
