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
#include <deal.II/lac/dynamic_sparsity_pattern.h>

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
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
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

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "my_table.h"
#include "ref_pt_list.h"
#include "MyParameterHandler.h"

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
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();
    void run2 ();
    void run2b ();
    void run2c ();

    double m_T[2];
    double m_W[5];
    double m_I12; 
  protected:
    int DoIter( string="" );

    void save( string );
    void save_one( string );
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    void do_superposition();
    void estimate_error( double& );
    void Interpolate_R_to_C( string );
    void compute_tangent();

    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_contributions();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );

    void output_results ( string, string = "step" );
    void output_guess ();

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    IndexSet locally_owned_dofs, locally_owned_dofs_2;
    IndexSet locally_relevant_dofs, locally_relevant_dofs_2;
    ConstraintMatrix constraints, constraints_2;

    LA::MPI::SparseMatrix m_system_matrix, m_system_matrix_2;
    LA::MPI::Vector m_system_rhs, m_system_rhs_2;
    LA::MPI::Vector m_newton_update;
    LA::MPI::Vector m_Psi_C;
    LA::MPI::Vector m_Psi_C_ghosted;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    MyTable m_table;
    MyTable m_results;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
    : 
    CBase<2>(xmlfilename),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe),2),
    dof_handler_2 (triangulation)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
    dof_handler_2.clear ();
  }
  
  #include "shared_1.h"
  #include "grid.h"  

  template<int dim>
  void MySolver<dim>::do_superposition()
  {
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    constraints.distribute (m_Psi_ref);
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);

    m_workspace_1 = m_Psi_1;
    m_workspace_2 = m_Psi_2;
    
    CPotential<dim> Potential_fct ( m_omega );
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
        fe_values.get_function_gradients( m_workspace_1, Psi_1_grad);
        fe_values.get_function_gradients( m_workspace_2, Psi_2_grad);

        for ( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          p12 = Psi_1[qp]*Psi_2[qp];
          p1q = Psi_1[qp]*Psi_1[qp];
          p2q = Psi_2[qp]*Psi_2[qp];
          Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu[0];

          T[0] += JxW*(Psi_1_grad[qp]*Psi_1_grad[qp] + Q*p1q);
          T[1] += JxW*(Psi_2_grad[qp]*Psi_2_grad[qp] + Q*p2q);
          I12  += JxW*(Psi_1_grad[qp]*Psi_2_grad[qp] + Q*p12);
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
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_computing_timer.enter_section(__func__);

    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    
    m_Psi_1.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_2.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_ref.reinit (locally_owned_dofs, mpi_communicator);
    m_newton_update.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi_1.local_size()  << endl;
     
    m_system_rhs=0;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    DynamicSparsityPattern csp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);
    
    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    m_Psi_C.reinit (locally_owned_dofs_2, mpi_communicator);
    m_Psi_C_ghosted.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);
    m_system_rhs_2.reinit(locally_owned_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2);
    constraints_2.close ();

    DynamicSparsityPattern csp_2 (locally_relevant_dofs_2);
    DoFTools::make_sparsity_pattern (dof_handler_2, csp_2, constraints_2, false);
    SparsityTools::distribute_sparsity_pattern (csp_2, dof_handler_2.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs_2);
    m_system_matrix_2.reinit (locally_owned_dofs_2, locally_owned_dofs_2, csp_2, mpi_communicator);    

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);
    
    m_workspace_1 = m_Psi_1;
    m_workspace_2 = m_Psi_2;
    
    CPotential<dim> Potential_fct ( m_omega );
    VectorTools::interpolate (dof_handler, Potential_fct, m_workspace_ng );    
    constraints.distribute(m_workspace_ng);
    m_Psi_ref = m_workspace_ng;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi_1");
    data_out.add_data_vector (m_workspace_2, "m_Psi_2");
    data_out.add_data_vector (m_Psi_ref, "Potential");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu", mpi_communicator );

    //DataOutBase::DataOutFilter data_filter(DataOutBase::DataOutFilterFlags(true,true));
    //data_out.write_hdf5_parallel(data_filter, "guess.h5", mpi_communicator);

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    string filename;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi_sol");
    data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.add_data_vector (subdomain, "subdomain");
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
    m_res_old[0] = m_res[0];
    m_counter=0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_2.add( -m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_2);

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

      if( m_root ) cout << m_table;
      if( m_res[0] < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

      m_counter++;
    }while( true );
    
    // Standard Newton 
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_ref.add( -0.1, m_newton_update); 
      constraints.distribute(m_Psi_ref);

      assemble_rhs();
      
      m_resp[0] = m_res_old[0]-m_res[0];
      m_res_old[0] = m_res[0];

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res[0] );
      m_table.insert( cols, MyTable::RESP, m_resp[0] );
      m_table.insert( cols, MyTable::MU, m_mu[0] );
      m_table.insert( cols, MyTable::GS, m_gs[0] );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N[0] );

      if( m_root ) cout << m_table;
      if( m_res[0] < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );
    
    m_N[0] = Particle_Number(m_Psi_ref);
    
    if( m_N[0] < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);
    
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    make_grid_custom();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    //CEigenfunctions<dim> Ef2( m_QN2, f );

    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    //VectorTools::interpolate (dof_handler, Ef2, m_workspace );
    //m_workspace *= 1.0/sqrt(Particle_Number(m_workspace));
    //m_Psi_1 += m_workspace;
    //m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
      
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu[0] = T/N+m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu[0] << endl;

    m_ti = sqrt((m_mu[0]*N-T)/(m_gs[0]*W)); // WARNING: NEHARI
    
    status = DoIter();
    
    pcout << "L2_norm of m_Psi_ref: " << m_N[0] << endl;


    if( status == Status::SUCCESS )
    {
      output_results("","final");
      save(  path + "final.bin" );
      Interpolate_R_to_C( path + "Cfinal.bin" );
      
      dump_info_xml();
    }

    if( m_root )
    {
      ofstream ofs("log.txt");
      ofs << m_table;
    }
  }

  template <int dim>
  void MySolver<dim>::run2 ()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    //make_grid();
    make_grid_custom();
    setup_system(true);

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

    status = DoIter();

    pcout << "L2_norm of m_Psi_ref: " << m_N[0] << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","final");
      dump_info_xml();
    }

    if( m_root )
    {
      ofstream ofs("log.txt");
      ofs << m_table;
    }
  }
  
  template <int dim>
  void MySolver<dim>::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu[0] = ceil(10.0*m_mu_0)/10.0 + m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu[0] << endl;

    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      // sqrt((m_mu[0]*N-T)/(m_gs[0]*W)); if m_Psi_2 == 0
      // sqrt((m_mu[0]*N-T)/(4.0*m_gs[0]*W)); if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu[0]*N-T)/(m_gs[0]*W));
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu[0] );
      m_results.insert( cols, MyTable::GS, m_gs[0] );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N[0] );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        save( path + "final.bin" );
        save_one( path + "final_one.bin" );
        vector<double> newgs = {m_gs[0]*m_N[0]}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        newgs[0] = m_gs[0];
        m_ph.Set_Physics( "gs_1", newgs );
        estimate_error(m_final_error);
        output_results(path,"final");
        dump_info_xml(path);
        Interpolate_R_to_C( path + "Cfinal.bin" );
        m_Psi_1 = m_Psi_ref;
        m_Psi_2 = 0;
      }
      else if( status == Status::SLOW_CONV )
      {
        m_Psi_2 = 0; 
      }
      else
      {
        break;
      }
      m_mu[0] += m_gs[0]/fabs(m_gs[0])*m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }

  template <int dim>
  void MySolver<dim>::run2c ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu[0] = ceil(10.0*m_mu_0)/10.0 + m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu[0] << endl;

    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      // sqrt((m_mu[0]*N-T)/(m_gs[0]*W)); if m_Psi_2 == 0
      // sqrt((m_mu[0]*N-T)/(4.0*m_gs[0]*W)); if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu[0]*N-T)/(m_gs[0]*W));
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu[0] );
      m_results.insert( cols, MyTable::GS, m_gs[0] );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N[0] );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        save( path + "final.bin" );
        save_one( path + "final_one.bin" );
        vector<double> newgs = {m_gs[0]*m_N[0]}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        newgs[0] = m_gs[0];
        m_ph.Set_Physics( "gs_1", newgs );
        estimate_error(m_final_error);
        output_results(path,"final");
        dump_info_xml(path);
        Interpolate_R_to_C( path + "Cfinal.bin" );

        compute_tangent();
        m_Psi_1 = m_Psi_ref;
        m_Psi_1.sadd( 0.1, m_workspace_ng );
        m_Psi_2 = 0;
      }
      else if( status == Status::SLOW_CONV )
      {
        m_Psi_2 = 0; 
      }
      else
      {
        break;
      }
      m_mu[0] += m_gs[0]/fabs(m_gs[0])*m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }  
  
  template<int dim>
  void MySolver<dim>::Interpolate_R_to_C( string filename )
  {
    m_computing_timer.enter_section(__func__);

    const QGauss<dim> quadrature_formula(fe.degree+1);
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    m_system_rhs_2=0;
    m_system_matrix_2=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);
    FEValues<dim> fe_values_2 (fe_2, quadrature_formula, update_values|update_JxW_values);

    const unsigned dofs_per_cell = fe_2.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    vector<double> vals(n_q_points);
    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    
    double JxW, Q1, tmp1;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator cell_2 = dof_handler_2.begin_active();
    for ( ; cell!=endc; ++cell, ++cell_2 )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values_2.reinit (cell_2);
        fe_values.get_function_values(m_workspace_1, vals);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values_2.JxW(qp);
          tmp1 = vals[qp]; 
          
          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) += JxW*tmp1*fe_values_2[rt].value(i,qp);
            for( unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j)+=JxW*(fe_values_2[rt].value(i,qp)*fe_values_2[rt].value(j,qp)+fe_values_2[it].value(i,qp)*fe_values_2[it].value(j,qp));
            }
          }
        }
        cell_2->get_dof_indices (local_dof_indices);
        constraints_2.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix_2, m_system_rhs_2);
      }
    }
    m_system_rhs_2.compress(VectorOperation::add);
    m_system_matrix_2.compress(VectorOperation::add);
    
    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix_2, m_Psi_C, m_system_rhs_2);
    constraints_2.distribute (m_Psi_C);
    
    m_Psi_C_ghosted=m_Psi_C;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler_2);
    solution_transfer.prepare_serialization(m_Psi_C_ghosted);
    triangulation.save( filename.c_str() );

/*    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler_2);
    data_out.add_data_vector (m_Psi_C_ghosted, "Psi");
    data_out.build_patches ();
    string filename = "C_final.vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);    
*/   
    m_computing_timer.exit_section();    
  }  

  template <int dim>
  void MySolver<dim>::compute_tangent ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    m_system_matrix=0;
    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_grad(n_q_points);
    vector<double> Psi(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi);
        fe_values.get_function_gradients(m_workspace_1, Psi_grad);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double pq = m_gs[0]*Psi[qp]*Psi[qp];
          double Q2 = Potential.value(fe_values.quadrature_point(qp)) - m_mu[0] + 3.0*pq;

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) -= JxW*Psi[qp]*fe_values.shape_value(i,qp);
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);
    m_system_matrix.compress(VectorOperation::add);

    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_workspace_ng, m_system_rhs);
    constraints.distribute (m_workspace_ng);    

    double N = Particle_Number(m_workspace_ng);
    m_workspace_ng = 1/sqrt(1+N); // +/- 1 / sqrt(1+N) 

    m_computing_timer.exit_section();

  }  

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    constraints.distribute(m_Psi_ref);
    m_workspace_1 = m_Psi_ref;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
  
  template<int dim>
  void MySolver<dim>::save_one( string filename )
  {
    double tmp = Particle_Number(m_Psi_ref);
    m_workspace_ng=m_Psi_ref;
    m_workspace_ng*=sqrt(1/tmp);
    m_workspace_1 = m_workspace_ng; 
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver<DIMENSION> solver("params.xml");
    solver.run2c ();
  }
return EXIT_SUCCESS;
}
