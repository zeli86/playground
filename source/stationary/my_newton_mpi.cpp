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
#include <deal.II/base/exceptions.h>
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

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "my_table.h"
//#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "muParser.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  #include "CBase.h"

  template <int dim>
  class MySolver : public CBase<dim,2>
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();
    void run2b ();
    void run2c ();

    double m_I[8]; 
  protected:
    int DoIter( string="" );

    void save( string );
    void save_one( string );
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void do_superposition();
    void estimate_error( double& );
    void compute_tangent();

    bool solve();
    void compute_contributions();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );

    void output_results ( string, string = "step" );
    void output_vector ( LA::MPI::Vector&, string );
    void output_guess ();

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_newton_update;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    MyTable m_table;
    MyTable m_results;

    double m_mu_punkt;

    using CBase<dim,2>::mpi_communicator;
    using CBase<dim,2>::m_xmin;
    using CBase<dim,2>::m_xmax;
    using CBase<dim,2>::m_ymin;
    using CBase<dim,2>::m_ymax;
    using CBase<dim,2>::m_zmin;
    using CBase<dim,2>::m_zmax;
    using CBase<dim,2>::m_root;
    using CBase<dim,2>::m_rank;
    using CBase<dim,2>::m_t;
    using CBase<dim,2>::m_ti;
    using CBase<dim,2>::m_N;
    using CBase<dim,2>::m_omega;
    using CBase<dim,2>::m_mu;
    using CBase<dim,2>::m_gs;
    using CBase<dim,2>::m_counter;
    using CBase<dim,2>::pcout;
    using CBase<dim,2>::m_computing_timer;
    using CBase<dim,2>::m_ph;
    using CBase<dim,2>::m_final_error;
    using CBase<dim,2>::m_NA;
    using CBase<dim,2>::m_Ndmu;
    using CBase<dim,2>::m_dmu;
    using CBase<dim,2>::m_QN1;
    using CBase<dim,2>::m_res;
    using CBase<dim,2>::m_resp;
    using CBase<dim,2>::m_res_old;
    using CBase<dim,2>::m_epsilon;
    using CBase<dim,2>::m_maxiter;
    using CBase<dim,2>::m_total_no_cells;
    using CBase<dim,2>::m_total_no_active_cells;
    using CBase<dim,2>::m_global_refinement;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
    : 
    CBase<dim,2>(xmlfilename),
    triangulation (CBase<dim,2>::mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
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

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1,dim>> Psi_1_grad(n_q_points);
    vector<Tensor<1,dim>> Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double locint[8] = {};

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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double p12 = Psi_1[qp]*Psi_2[qp];
          double p1q = Psi_1[qp]*Psi_1[qp];
          double p2q = Psi_2[qp]*Psi_2[qp];
          double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu;

          locint[0] += JxW*p1q*p1q;
          locint[1] += JxW*p2q*p2q;
          locint[2] += JxW*p1q*p2q;
          locint[3] += JxW*p1q*p12;
          locint[4] += JxW*p2q*p12;
          locint[5] += JxW*(Psi_1_grad[qp]*Psi_1_grad[qp] + Q*p1q);
          locint[6] += JxW*(Psi_2_grad[qp]*Psi_2_grad[qp] + Q*p2q);
          locint[7] += JxW*(Psi_1_grad[qp]*Psi_2_grad[qp] + Q*p12);
        }  
      }
    }
    for( int i=0; i<5; i++ ) 
      locint[i] *= m_gs;

    MPI_Allreduce( locint, m_I, 8, MPI_DOUBLE, MPI_SUM, mpi_communicator);
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
    m_workspace_3.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
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
    data_out.build_patches (gl_subdivisions);
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
    data_out.build_patches (gl_subdivisions);

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }

  template <int dim>
  void MySolver<dim>::output_vector ( LA::MPI::Vector& vec, string filename )
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(vec);
    m_workspace_1=vec;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "vec");
    data_out.build_patches (gl_subdivisions);
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }  

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    static int bla=0;

    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    
    CPotential<dim> Potential( m_omega );

    do_superposition();
    m_workspace_1 = m_Psi_ref;
    MyRealTools::MPI::AssembleRHS_L2gradient<dim>( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_res, m_system_rhs );
    m_res_old = m_res;
    m_counter=0;
    bool bempty;

    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      m_workspace_1 = m_Psi_ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_system_matrix);
      bool bsucc = solve();
      if( !bsucc ) { return Status::SINGULAR; }

      pcout << "m_newton_update.l2_norm() = " << m_newton_update.l2_norm() << endl;

/*
      m_workspace_1 = m_Psi_1;
      m_workspace_2 = m_newton_update;
      MyRealTools::MPI::orthonormalize(mpi_communicator, dof_handler, fe, constraints, m_workspace_2, m_workspace_1, m_workspace_ng );
      m_newton_update = m_workspace_ng;
*/
      double tau;
      m_workspace_2 = m_newton_update;
      MyRealTools::MPI::compute_stepsize( mpi_communicator, dof_handler, fe, Potential, m_workspace_1, m_workspace_2, m_mu, m_gs, tau );

      m_Psi_2.add( tau*m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_2);

      bempty = this->find_ortho_min();

      if( bempty ) return Status::FAILED;
        
      do_superposition();
      m_workspace_1 = m_Psi_ref;      
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_res, m_system_rhs );
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::STEPSIZE, tau );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );

      m_counter++;
      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( this->l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }
      if( m_counter == m_maxiter ) { retval=Status::MAXITER; break; }
      if( isnan(m_res) ) { retval=retval=Status::FAILED; break; }
   
    }while( true );
    
    // Standard Newton 
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      m_workspace_1 = m_Psi_ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_system_matrix);
      bool bsucc = solve();
      if( !bsucc ) { return Status::SINGULAR; }

      double tau;
      m_workspace_2 = m_newton_update;
      MyRealTools::MPI::compute_stepsize( mpi_communicator, dof_handler, fe, Potential, m_workspace_1, m_workspace_2, m_mu, m_gs, tau );

      m_Psi_ref.add( tau*m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_ref);

      m_workspace_1 = m_Psi_ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_res, m_system_rhs );
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::STEPSIZE, tau );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      m_counter++;

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }
      if( m_resp < 0 ) { retval=Status::SUCCESS; break; }
      if( m_counter == m_maxiter ) { retval=Status::MAXITER; break; }
    }while( true );
    
    m_workspace_1 = m_Psi_ref;
    m_N = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
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
    m_workspace_1 = m_Psi_1;
    m_Psi_1 *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 ));
    //VectorTools::interpolate (dof_handler, Ef2, m_workspace );
    //m_workspace *= 1.0/sqrt(Particle_Number(m_workspace));
    //m_Psi_1 += m_workspace;
    //m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
      
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    m_ti = sqrt((m_mu*N-T)/(m_gs*W)); // WARNING: NEHARI
    
    status = DoIter();
    
    pcout << "L2_norm of m_Psi_ref: " << m_N << endl;


    if( status == Status::SUCCESS )
    {
      output_results("","final");
      save(  path + "final.bin" );
      
      this->dump_info_xml();
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
    //CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_workspace_1 = m_Psi_1;
    m_Psi_1 *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 ));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    //output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      // sqrt((m_mu*N-T)/(m_gs*W)); if m_Psi_2 == 0
      // sqrt((m_mu*N-T)/(4.0*m_gs*W)); if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu*N-T)/(m_gs*W));
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        estimate_error(m_final_error);

        save( path + "final.bin" );
        save_one( path + "final_one.bin" );
        vector<double> newgs = {m_gs*m_N}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        newgs[0] = m_gs;
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params.xml" );       
        
        output_results(path,"final");
        this->dump_info_xml(path);
        m_Psi_1 = m_Psi_ref;
        m_Psi_2 = 0;
      }
      else if( status == Status::MAXITER || status == Status::SINGULAR )
      {
        m_Psi_ref = m_Psi_1;
        m_Psi_2 = 0; 
      }
      else
      {
        //m_Psi_2 = 0; 
        break;
      }
      m_mu += m_gs/fabs(m_gs)*m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;      
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
    //CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_workspace_1 = m_Psi_1;
    m_Psi_1 *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 ));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    //output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      // sqrt((m_mu*N-T)/(m_gs*W)); if m_Psi_2 == 0
      // sqrt((m_mu*N-T)/(4.0*m_gs*W)); if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu*N-T)/(m_gs*W));
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        estimate_error(m_final_error);

        save( path + "final.bin" );
        save_one( path + "final_one.bin" );
        vector<double> newgs = {m_gs*m_N}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        newgs[0] = m_gs;
        m_ph.Set_Physics( "gs_1", newgs );

        output_results(path,"final");
        this->dump_info_xml(path);

        compute_tangent();
        m_Psi_1 = m_Psi_ref;
        m_Psi_1.sadd( m_t[1]/fabs(m_t[1])*m_dmu, m_workspace_ng );
        m_Psi_2 = 0;
        //output_vector( m_workspace_ng, "tangent_" + to_string(i) + ".vtu"  );
        //output_vector( m_Psi_1, "Psi_1_" + to_string(i) + ".vtu"  );
      }
      else if( status == Status::MAXITER || status == Status::SINGULAR )
      {
        m_Psi_ref = m_Psi_1;
        m_Psi_2 = 0; 
      }
      else
      {
        m_Psi_2 = 0; 
        //break;
      }
//      m_mu += m_gs/fabs(m_gs)*0.02*m_mu_punkt;
      m_mu += m_gs/fabs(m_gs)*m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }  
    
  template <int dim>
  void MySolver<dim>::compute_tangent ()
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    CPotential<dim> Potential( m_omega );
    MyRealTools::MPI::AssembleSystem_tangent( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_system_matrix, m_system_rhs );
/*
    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_workspace_ng, m_system_rhs);
    constraints.distribute (m_workspace_ng);    
*/
    m_workspace_ng = 0;
    SolverControl solver_control (m_newton_update.size(), 1e-15);
    //PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
    PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);
    
    PETScWrappers::PreconditionBlockJacobi::AdditionalData adata;
    PETScWrappers::PreconditionBlockJacobi preconditioner(m_system_matrix,adata);    

    try 
    {
      solver.solve(m_system_matrix, m_workspace_ng, m_system_rhs, preconditioner);
    }
    catch( ExceptionBase& e )
    {
      //pcout << e.what() << endl;
      pcout << "Possible singular matrix!" << endl;
    }
    constraints.distribute (m_workspace_ng);

    m_workspace_1 = m_workspace_ng;
    double N = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
    //m_mu_punkt = 1/sqrt(1+N);    
    //m_workspace_ng *= m_mu_punkt; // +/- 1 / sqrt(1+N) 

    //pcout << "m_mu_punkt = " << m_mu_punkt << endl;
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
    m_workspace_1 = m_Psi_ref;
    double tmp = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
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
  //deallog.depth_console (2);

  MyParameterHandler params("params.xml");
  int dim=0;

  try
  {
    dim = int(params.Get_Mesh("DIM",0));
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    switch(dim)
    {
      case 2: { BreedSolver::MySolver<2> solver("params.xml");
                solver.run2b();
                break; }
      case 3: { BreedSolver::MySolver<3> solver("params.xml");
                solver.run2b();
                break; }
      default: cout << "You have found a new dimension!" << endl;
    }
  }
return EXIT_SUCCESS;
}
