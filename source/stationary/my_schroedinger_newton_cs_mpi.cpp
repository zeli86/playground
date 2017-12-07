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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "my_table.h"
#include "MyParameterHandler.h"
#include "CBase2.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  class MySolver : public CBase2
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run2b ();

    void set_t( const double a, const double b ) { m_t[0]=a; m_t[1]=b; };
    
    MPI_Comm mpi_communicator;

  protected:
    double l2norm_t() { return sqrt(m_t[0]*m_t[0]+m_t[1]*m_t[1]); };  
  
    int DoIter( string="" );

    void save( const string& );
    void load( const string& );
    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    void do_superposition();
    void do_superposition_U();
    void estimate_error( double& );
    void compute_U( LA::MPI::Vector&, LA::MPI::Vector&, LA::MPI::Vector& );
    void compute_contributions();

    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double&, double&, double& );

    void output_results ( string, string = "step" );
    void output_guess ();

    parallel::distributed::Triangulation<2> triangulation;
    FE_Q<2> fe;
    DoFHandler<2> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_newton_update;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_U;
    LA::MPI::Vector m_U1;
    LA::MPI::Vector m_U12;
    LA::MPI::Vector m_U2;
    LA::MPI::Vector m_ext_grav_pot;    
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;
    LA::MPI::Vector m_workspace_4;
    LA::MPI::Vector m_workspace_5;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    FunctionParser<2> m_guess;

    MyTable m_table;
    MyTable m_results;

    double m_upsilon;
  };

/*************************************************************************************************/
/**
 * Constructor
 */
  MySolver::MySolver ( const std::string xmlfilename ) 
    : 
    CBase2(xmlfilename),
    mpi_communicator(MPI_COMM_WORLD), 
    triangulation (mpi_communicator, Triangulation<2>::MeshSmoothing(Triangulation<2>::limit_level_difference_at_vertices|Triangulation<2>::eliminate_refined_inner_islands|Triangulation<2>::smoothing_on_refinement|Triangulation<2>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    m_guess(1)
  {
    std::map <std::string, double> constants;
    std::string variables = "x,y";

    try 
    { 
      constants["z_0"] = m_ph.Get_Physics("z_0",0);
      constants["M_L"] = m_ph.Get_Physics("M_L",0);
      constants["C"] = m_ph.Get_Physics("C",0);
      constants["mu_0"] = m_ph.Get_Physics("mu_0",0);
      constants["m"] = m_ph.Get_Physics("m",0);
      m_upsilon = 0.001;
      m_guess_str =  m_ph.Get_Parameter("guess");
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }

    m_guess.initialize(variables, m_guess_str, constants);
  }

  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }  

  void MySolver::do_superposition()
  {
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    constraints.distribute (m_Psi_ref);
  }

  void MySolver::do_superposition_U()
  {
    m_U=0;
    m_U.add(m_t[0]*m_t[0],m_U1,m_t[1]*m_t[1],m_U2);
    m_U.add(2*m_t[0]*m_t[1],m_U12);
    constraints.distribute (m_U);
  }

  void MySolver::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
    Point<2,double> pt1( m_xmin, m_ymin );
    Point<2,double> pt2( m_xmax, m_ymax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);

    m_computing_timer.exit_section();
  }

  void MySolver::setup_system()
  {
    m_computing_timer.enter_section(__func__);

    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    
    m_Psi_1.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_2.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_ref.reinit (locally_owned_dofs, mpi_communicator);
    m_ext_grav_pot.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_U.reinit (locally_owned_dofs, mpi_communicator);
    m_U1.reinit (locally_owned_dofs, mpi_communicator);
    m_U12.reinit (locally_owned_dofs, mpi_communicator);
    m_U2.reinit (locally_owned_dofs, mpi_communicator);
    m_newton_update.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_3.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_4.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_5.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi_1.local_size()  << endl;
     
    m_system_rhs=0;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<2>(), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    
    m_computing_timer.exit_section();
  }

  void MySolver::output_guess ()
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);
    
    m_workspace_1 = m_Psi_1;
    m_workspace_2 = m_Psi_2;
   
    
    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi_1");
    data_out.add_data_vector (m_workspace_2, "m_Psi_2");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu", mpi_communicator );

    m_computing_timer.exit_section();
  }

  void MySolver::output_results ( string path, string prefix )
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    constraints.distribute(m_U);
    m_workspace_2=m_U;

    m_workspace_3=m_newton_update;

    string filename;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi");
    data_out.add_data_vector (m_workspace_2, "U");
    data_out.add_data_vector (m_workspace_3, "L2 grad");
    data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }

  void MySolver::compute_contributions()
  {
    m_computing_timer.enter_section(__func__);

    compute_U( m_Psi_1, m_Psi_1, m_U1 );
    compute_U( m_Psi_1, m_Psi_2, m_U12 );
    compute_U( m_Psi_2, m_Psi_2, m_U2 );

    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);

    m_workspace_1 = m_Psi_1;
    m_workspace_2 = m_Psi_2;
    m_workspace_3 = m_U1;
    m_workspace_4 = m_U2;
    m_workspace_5 = m_U12;
    
    const QGauss<2> quadrature_formula(fe.degree+1);
    FEValues<2> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<double> Psi_1(n_q_points), Psi_2(n_q_points);
    vector<double> U1(n_q_points), U2(n_q_points), U12(n_q_points);
    vector<double> Ext_Grav_Pot(n_q_points);
    vector<Tensor<1,2> > Psi_1_grad(n_q_points);
    vector<Tensor<1,2> > Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, p12, p1q, p2q, Q;
    double I12 = 0.0;
    double T[2] = {};
    double W[5] = {};
    double V2[9] = {};

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, Psi_1 );
        fe_values.get_function_values( m_workspace_2, Psi_2 );
        fe_values.get_function_gradients( m_workspace_1, Psi_1_grad);
        fe_values.get_function_gradients( m_workspace_2, Psi_2_grad);
        fe_values.get_function_values( m_workspace_3, U1);
        fe_values.get_function_values( m_workspace_4, U2);
        fe_values.get_function_values( m_workspace_5, U12);
        fe_values.get_function_values( m_ext_grav_pot, Ext_Grav_Pot);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          p12 = Psi_1[qp]*Psi_2[qp];
          p1q = Psi_1[qp]*Psi_1[qp];
          p2q = Psi_2[qp]*Psi_2[qp];
          Q = Ext_Grav_Pot[qp] - m_mu;

          T[0] += JxW*(Psi_1_grad[qp]*Psi_1_grad[qp] + Q*p1q);
          T[1] += JxW*(Psi_2_grad[qp]*Psi_2_grad[qp] + Q*p2q);
          I12  += JxW*(Psi_1_grad[qp]*Psi_2_grad[qp] + Q*p12);
          W[0] += JxW*p1q*p1q;
          W[1] += JxW*p2q*p2q;
          W[2] += JxW*p1q*p2q;
          W[3] += JxW*p1q*p12;
          W[4] += JxW*p2q*p12;
          V2[0] += JxW*U1[qp]*p1q;
          V2[1] += JxW*2*U12[qp]*p1q;
          V2[2] += JxW*U2[qp]*p1q;
          V2[3] += JxW*U1[qp]*p2q;
          V2[4] += JxW*2*U12[qp]*p2q;
          V2[5] += JxW*U2[qp]*p2q;
          V2[6] += JxW*U1[qp]*p12;
          V2[7] += JxW*U12[qp]*p12;
          V2[8] += JxW*U2[qp]*p12;
        }  
      }
    }
    
    for( int i=0; i<5; i++ ) W[i] *= m_gs[0];
    
    MPI_Allreduce( T, m_T, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &I12, &m_I12, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( W, m_W, 5, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( V2, m_V2, 9, MPI_DOUBLE, MPI_SUM, mpi_communicator);
/*    
    pcout << "T[1] := " << m_T[0] << ";"<< endl;
    pcout << "T[2] := " << m_T[1] << ";"<< endl;
    pcout << "I12 := " << m_I12 << ";" << endl;
    for( int i=0; i<5; i++ )
      pcout << "W[" << to_string(i+1) << "] := " << m_W[i] << ";" << endl;
    for( int i=0; i<9; i++ )
      pcout << "V2[" << to_string(i+1) << "] := " << m_V2[i] << ";" << endl;
 */
    m_computing_timer.exit_section();
  }

  int MySolver::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    
    assemble_rhs();
    m_res_old = m_res;
    m_counter=0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();
     
      m_Psi_2.add( -0.01*m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_2);

      double err;
      estimate_error(err);    
      output_results(path);

      double T, T2, N, W, V2;
      compute_E_lin( m_Psi_ref, T, T2, N, W, V2 );
      pcout << "T = " << T << endl;
      pcout << "T2 = " << T2 << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "V2 = " << V2 << endl;

      find_new_t();
      assemble_rhs();
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

//      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs[0] );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      //if( l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

      m_counter++;
    }while( true );

    m_N = Particle_Number(m_Psi_ref);
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);

    pcout << "m_N == " << m_N << endl;
    pcout << "retval == " << retval << endl;
    
    return retval;
  }

  void MySolver::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, T2, N, W, V2;
    int status;

    load( "grav_pot.bin" );

    VectorTools::interpolate (dof_handler, m_guess, m_Psi_1 );

    m_Psi_1 *= sqrt(1.0/Particle_Number(m_Psi_1));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, T2, N, W, V2 );
    
    double m_mu_0 = T + m_gs[0]*W + V2;

    double sign_dmu = (m_gs[0]*W+V2)/fabs(m_gs[0]*W+V2);
    m_dmu *= sign_dmu;
    //m_dmu *= -1;

    pcout << "m_dmu = " << m_dmu << endl;
    pcout << "T = " << T << endl;
    pcout << "T2 = " << T2 << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "V2 = " << V2 << endl;
    pcout << "m_mu_0 == " << m_mu_0 << endl;
 
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_dmu;

    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      m_ti = sqrt(m_dmu*N/(m_gs[0]*W+V2));
      //m_ti = 1;
      pcout << "m_ti = " << m_ti << endl;
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs[0] );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        save( path + "final.bin" );
        estimate_error(m_final_error);
        output_results(path,"final");
        dump_info_xml(path);
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
      m_mu += m_gs[0]/fabs(m_gs[0])*m_dmu;
      compute_E_lin( m_Psi_ref, T, T2, N, W, V2 ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }

  void MySolver::save( const string& filename )
  {
    constraints.distribute(m_Psi_ref);
    parallel::distributed::SolutionTransfer<2,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi_ref);

    triangulation.save( filename.c_str() );
  }

  void MySolver::load( const string& filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    parallel::distributed::Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for (unsigned f=0; f < GeometryInfo<2>::faces_per_cell; ++f)
      {
        const Point<2> face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary() && face_center[0]==0 )  
        {  
          cell->face(f)->set_all_boundary_ids(1);
        }
      }
    }   

    setup_system();

    parallel::distributed::SolutionTransfer<2,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_workspace_ng);

    constraints.distribute(m_workspace_ng);
    m_ext_grav_pot = m_workspace_ng;
  }

  void MySolver::compute_E_lin( LA::MPI::Vector& vec, double& T, double& T2, double& N, double& W, double& V2 )
  {
    m_computing_timer.enter_section(__func__);
    
    do_superposition();
    compute_U(m_Psi_1,m_Psi_1,m_U);
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    m_workspace_2 = m_U;
    
    const QGauss<2>  quadrature_formula(fe.degree+1);
    FEValues<2> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<double> vec_vals_U(n_q_points);
    vector<double> Ext_Grav_Pot(n_q_points);
    vector<Tensor<1,2>> vec_grad(n_q_points);

    double JxW, vec_val_q;
    
    double tmp1[5]={}, res[5]={};
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vec_vals);
        fe_values.get_function_values(m_workspace_2, vec_vals_U);
        fe_values.get_function_values(m_ext_grav_pot, Ext_Grav_Pot);
        fe_values.get_function_gradients(m_workspace_1, vec_grad);
        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          vec_val_q = vec_vals[qp]*vec_vals[qp];
          tmp1[0] += JxW*vec_grad[qp]*vec_grad[qp];
          tmp1[1] += JxW*(Ext_Grav_Pot[qp] + vec_vals_U[qp])*vec_val_q;
          tmp1[2] += JxW*vec_val_q;
          tmp1[3] += JxW*vec_val_q*vec_val_q;
          tmp1[4] += JxW*vec_vals_U[qp]*vec_val_q;
        }
      }
    }
    MPI_Allreduce( tmp1, res, 5, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    T=res[0]; T2 = res[1]; N=res[2]; W=res[3]; V2=res[4]; 
    m_computing_timer.exit_section();
  }

  double MySolver::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0, retval=0;
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    const QGauss<2>  quadrature_formula(fe.degree+1);
    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );      

        for ( unsigned qp=0; qp<n_q_points; qp++ )
          tmp1 += fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*vec_vals[qp]*vec_vals[qp];
      }
    }
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }
  
  void MySolver::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);

    const QGauss<2> quadrature_formula(fe.degree+1);
/*
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    constraints.distribute(m_U);
    m_workspace_2=m_U;
*/ 
   
    m_system_matrix=0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1,2>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
    vector<double> U(n_q_points);
    vector<double> Ext_Grav_Pot(n_q_points);

    double JxW, Q2, tmp;
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi_ref);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_values(m_ext_grav_pot, Ext_Grav_Pot);
        fe_values.get_function_gradients(m_workspace_1, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          Q2 =  U[qp] + Ext_Grav_Pot[qp] - m_mu + 3.0*m_gs[0]*Psi_ref[qp]*Psi_ref[qp];

          for ( unsigned i=0; i<dofs_per_cell; ++i )
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, m_system_matrix);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }

  void MySolver::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);

    const QGauss<2> quadrature_formula(fe.degree+1);
    
    do_superposition();
    do_superposition_U();
    
    m_workspace_1=m_Psi_ref;
    m_workspace_2=m_U;

    m_system_rhs=0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, 2> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
    vector<double> U(n_q_points);
    vector<double> Ext_Grav_Pot(n_q_points);

    double JxW, Q1, pq;
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi_ref);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_values(m_ext_grav_pot, Ext_Grav_Pot);
        fe_values.get_function_gradients(m_workspace_1, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          pq = Psi_ref[qp]*Psi_ref[qp];
          Q1 = U[qp] + Ext_Grav_Pot[qp] - m_mu +  m_gs[0]*pq;

          for ( unsigned i=0; i<dofs_per_cell; ++i )
            cell_rhs(i) += JxW*(Psi_ref_grad[qp]*fe_values.shape_grad(i,qp) + Q1*Psi_ref[qp]*fe_values.shape_value(i,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_res = m_system_rhs.l2_norm();
    m_computing_timer.exit_section();
  }

  void MySolver::estimate_error ( double& err )
  {
    m_computing_timer.enter_section(__func__);
    
    const QGauss<2> quadrature_formula(fe.degree+1);
   
    compute_U(m_Psi_ref,m_Psi_ref,m_U);
    m_workspace_2=m_U;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<double> U(n_q_points);
    vector<Tensor<1,2>> grads(n_q_points);
    vector<double> Ext_Grav_Pot(n_q_points);

    double JxW, Q1, pq, tmp;
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_gradients(m_workspace_1, grads);
        fe_values.get_function_values(m_ext_grav_pot, Ext_Grav_Pot);        

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          Q1 = U[qp] + Ext_Grav_Pot[qp] - m_mu + m_gs[0]*(vals[qp]*vals[qp]);

          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; j++ )
              cell_matrix(i,j) += JxW*(fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_system_matrix.compress(VectorOperation::add);   

    solve();
    
    m_workspace_1=m_newton_update;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<2>(2), m_error_per_cell, QGauss<2>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD));     
  
    m_computing_timer.exit_section();
  }

  void MySolver::compute_U ( LA::MPI::Vector& vec1, LA::MPI::Vector& vec2, LA::MPI::Vector& U )
  {
    m_computing_timer.enter_section(__func__);
    
    const QGauss<2> quadrature_formula(fe.degree+1);
   
    constraints.distribute(vec1);
    m_workspace_1=vec1;

    constraints.distribute(vec2);
    m_workspace_2=vec2;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);

    double JxW;
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals_1);
        fe_values.get_function_values(m_workspace_2, vals_2);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];

          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) -= JxW*vals_1[qp]*vals_2[qp]*fe_values.shape_value(i,qp);
            for ( unsigned j=0; j<dofs_per_cell; j++ )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_system_matrix.compress(VectorOperation::add);   

/*
    U=0;
    SolverControl solver_control ( m_system_rhs.size(), 1e-15 );
    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);
    PETScWrappers::PreconditionBlockJacobi preconditioner(m_system_matrix);
    
    cg.solve (m_system_matrix, U, m_system_rhs, preconditioner);
    constraints.distribute (U);
*/

    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, U, m_system_rhs);
    constraints.distribute (U);

    U *= m_upsilon;

    m_computing_timer.exit_section();
  }
 
  void MySolver::solve ()
  {
    m_computing_timer.enter_section(__func__);
    pcout << "Solving..." << endl;
    
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_newton_update, m_system_rhs);
    constraints.distribute (m_newton_update);

    m_computing_timer.exit_section();
  }

} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver solver("params.xml");
    solver.run2b ();
  }
return EXIT_SUCCESS;
}
