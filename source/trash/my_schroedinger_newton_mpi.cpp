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
#include "CBase.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  template <int dim>
  class MySolver; // class forward declaration

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  using namespace std;
  using namespace dealii;
  
  template<int dim>
  int fun_f( const gsl_vector* x, void* params, gsl_vector* f )
  {
    MySolver<dim>* s = reinterpret_cast<MySolver<dim>*>(params); 
    double f0, f1;
    const double t0 = gsl_vector_get (x, 0);
    const double t1 = gsl_vector_get (x, 1);
    const double t0q = t0*t0;
    const double t1q = t1*t1;
    const double t0k = t0q*t0;
    const double t1k = t1q*t1;

    f0 = t0*s->m_T[0] + t1*s->m_I12;
    f0 += t0k*s->m_W[0] + 3*t0*t1q*s->m_W[2] + 3*t0q*t1*s->m_W[3] + t1k*s->m_W[4];
    f0 += t0k*s->m_V2[0] + t0q*t1*s->m_V2[6] + t0*t1q*s->m_V2[1] + t1k*s->m_V2[7] + t0q*t1*s->m_V2[2] + t0*t1q*s->m_V2[8];
    f1 = t1*s->m_T[1] + t0*s->m_I12; 
    f1 += t1k*s->m_W[1] + 3*t0*t1q*s->m_W[4] + 3*t0q*t1*s->m_W[2] + t0k*s->m_W[3];
    f1 += t0k*s->m_V2[6] + t0q*t1*s->m_V2[3] + t0*t1q*s->m_V2[7] + t1k*s->m_V2[4] + t0q*t1*s->m_V2[8] + t0*t1q*s->m_V2[5];
    
    gsl_vector_set (f, 0, f0);
    gsl_vector_set (f, 1, f1);
  return GSL_SUCCESS;
  }

  template<int dim>
  int fun_df (const gsl_vector* x, void* params, gsl_matrix* J )
  {
    MySolver<dim>* s = reinterpret_cast<MySolver<dim>*>(params); 
    double df00, df01, df10, df11;
    const double t0 = gsl_vector_get (x, 0);
    const double t1 = gsl_vector_get (x, 1);
    const double t0q = t0*t0;
    const double t1q = t1*t1;
    const double t0k = t0q*t0;
    const double t1k = t1q*t1;

    // df0/dt0
    df00 = s->m_T[0] + 3*t0q*s->m_W[0] + 3*t1q*s->m_W[2] + 6*t0*t1*s->m_W[3] + 3*t0q*s->m_V2[0] + 2*t0*t1*s->m_V2[6] + t1q*s->m_V2[1] + 2*t0*t1*s->m_V2[2] + t1q*s->m_V2[8];
    // df0/dt1
    df01 = s->m_I12 + 6*t0*t1*s->m_W[2] + 3*t0q*s->m_W[3] + 3*t1q*s->m_W[4] + t0q*s->m_V2[6] + 2*t0*t1*s->m_V2[1] + 3*t1q*s->m_V2[7] + t0q*s->m_V2[2] + 2*t0*t1*s->m_V2[8];
    // df1/dt0
    df10 = s->m_I12 + 3*t1q*s->m_W[4] + 6*t0*t1*s->m_W[2] + 3*t0q*s->m_W[3] + 3*t0q*s->m_V2[6] + 2*t0*t1*s->m_V2[3] + t1q*s->m_V2[7] + 2*t0*t1*s->m_V2[8] + t1q*s->m_V2[5];
    // df1/dt1
    df11 = s->m_T[1] + 3*t1q*s->m_W[1] + 6*t0*t1*s->m_W[4] + 3*t0q*s->m_W[2] + t0q*s->m_V2[3] + 2*t0*t1*s->m_V2[7] + 3*t1q*s->m_V2[4] + t0q*s->m_V2[8] + 2*t0*t1*s->m_V2[5];

    gsl_matrix_set (J, 0, 0, df00);
    gsl_matrix_set (J, 0, 1, df01);
    gsl_matrix_set (J, 1, 0, df10);
    gsl_matrix_set (J, 1, 1, df11);
    
//    printf( "f0 = %.10e, f1 = %.10e\n", f0, f1 );
  return GSL_SUCCESS;
  }

  template<int dim>
  int fun_fdf ( const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J )
  {
    fun_f<dim> (x, params, f);
    fun_df<dim> (x, params, J);
  return GSL_SUCCESS;
  }  
  
  template <int dim>
  class MySolver : public CBase2<2>
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();
    void run2b ();

    //void compute_J( double& );
    void set_t( const double a, const double b ) { m_t[0]=a; m_t[1]=b; };
    
    MPI_Comm mpi_communicator;

    double m_T[2];
    double m_W[5];
    double m_V2[9];
    double m_I12; 
  protected:
    double l2norm_t() { return sqrt(m_t[0]*m_t[0]+m_t[1]*m_t[1]); };  
  
    int DoIter( string="" );

    void save( string );
    void save_one( string );
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    void do_superposition();
    void do_superposition_U();
    //void find_min_J();
    void estimate_error( double& );
    void compute_U( LA::MPI::Vector&, LA::MPI::Vector&, LA::MPI::Vector& );
    void Interpolate_R_to_C( string );
    void compute_contributions();

    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double&, double& );
    //void find_min_J();

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
    LA::MPI::Vector m_U;
    LA::MPI::Vector m_U1;
    LA::MPI::Vector m_U12;
    LA::MPI::Vector m_U2;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;
    LA::MPI::Vector m_workspace_4;
    LA::MPI::Vector m_workspace_5;
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
    CBase2<2>(xmlfilename),
    mpi_communicator(MPI_COMM_WORLD), 
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe),2),
    dof_handler_2 (triangulation)
  {
    m_fun.f=&fun_f<dim>;
    m_fun.df=&fun_df<dim>;
    m_fun.fdf=&fun_fdf<dim>;
    m_fun.params=reinterpret_cast<void*>(this);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
    dof_handler_2.clear ();
  }  

  template<int dim>
  void MySolver<dim>::do_superposition()
  {
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    constraints.distribute (m_Psi_ref);
  }

  template<int dim>
  void MySolver<dim>::do_superposition_U()
  {
    m_U=0;
    m_U.add(m_t[0]*m_t[0],m_U1,m_t[1]*m_t[1],m_U2);
    m_U.add(2*m_t[0]*m_t[1],m_U12);
    constraints.distribute (m_U);
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
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
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

    CompressedSimpleSparsityPattern csp_2 (locally_relevant_dofs_2);
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
  void MySolver<dim>::compute_contributions()
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
    
    CPotential<dim> Potential_fct ( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<double> Psi_1(n_q_points), Psi_2(n_q_points);
    vector<double> U1(n_q_points), U2(n_q_points), U12(n_q_points);
    vector<Tensor<1, dim> > Psi_1_grad(n_q_points);
    vector<Tensor<1, dim> > Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, p12, p1q, p2q, Q;
    double I12 = 0.0;
    double T[2] = {};
    double W[5] = {};
    double V2[9] = {};

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
        fe_values.get_function_values( m_workspace_3, U1);
        fe_values.get_function_values( m_workspace_4, U2);
        fe_values.get_function_values( m_workspace_5, U12);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          p12 = Psi_1[qp]*Psi_2[qp];
          p1q = Psi_1[qp]*Psi_1[qp];
          p2q = Psi_2[qp]*Psi_2[qp];
          Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu;

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

    pcout << "T[1] := " << m_T[0] << ";"<< endl;
    pcout << "T[2] := " << m_T[1] << ";"<< endl;
    pcout << "I12 := " << m_I12 << ";" << endl;
    for( int i=0; i<5; i++ )
      pcout << "W[" << to_string(i+1) << "] := " << m_W[i] << ";" << endl;
    for( int i=0; i<9; i++ )
      pcout << "V2[" << to_string(i+1) << "] := " << m_V2[i] << ";" << endl;

    m_computing_timer.exit_section();
  }
  /*
  template <int dim>
  void MySolver<dim>::compute_J( double& retval )
  {
    m_computing_timer.enter_section(__func__);

    do_superposition();
    do_superposition_U();

    m_workspace_1 = m_Psi_ref;
    m_workspace_2 = m_U;
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<double> U(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);

    double psum=0, uq;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          uq = vals[qp]*vals[qp];
          psum += fe_values.JxW(qp)*( grads[qp]*grads[qp] + (Potential.value(fe_values.quadrature_point(qp))-m_mu+U[qp]+0.5*m_gs*uq)*uq ); 
        }
      }
    }
    psum *= 0.5;
    MPI_Allreduce( &psum, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }
  */
  /*
  template<int dim>
  void MySolver<dim>::find_min_J()
  {
    compute_U( m_Psi_1, m_Psi_1, m_U1 );
    compute_U( m_Psi_1, m_Psi_2, m_U12 );
    compute_U( m_Psi_2, m_Psi_2, m_U2 );
    
    size_t iter=0;
    int status;
    double size;

    gsl_multimin_function my_func;
    my_func.n = 2;
    my_func.f = fun_neu<dim>;
    my_func.params = this;

    gsl_vector * ss = gsl_vector_alloc (2);
    gsl_vector * x = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);
    //gsl_vector_set (x, 0, m_t[0]);
    //gsl_vector_set (x, 1, m_t[1]);
    gsl_vector_set (x, 0, -5.929406695);
    gsl_vector_set (x, 1, 2.804554713);

    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 2);

    gsl_multimin_fminimizer_set (s, &my_func, x, ss);

    do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      if (status) break;


      printf ("%7d %.7f %.7f %.7f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->f);

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size(size,1e-6);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    m_t[0] = gsl_vector_get (s->x, 0);
    m_t[1] = gsl_vector_get (s->x, 1);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free(ss);
    gsl_vector_free (x);
  }    
  */
  template <int dim>
  int MySolver<dim>::DoIter ( string path )
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
     
      m_Psi_2.add( -m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_2);

      compute_contributions();
     
      gsl_vector* x = gsl_vector_alloc (2);
      gsl_vector* f = gsl_vector_alloc (2);
      
      gsl_vector_set(x, 0, .2064756745 );
      gsl_vector_set(x, 1, .8249265324 );
      
      fun_f<2>( x, this, f );
      
      printf( "f = %.10e %10e\n", gsl_vector_get( f, 0 ), gsl_vector_get( f, 1 ) ); 
      
      gsl_vector_free(x);
      gsl_vector_free(f);

      //find_new_t();      
      return 0;
        
      pcout << "xxx--------------------------------------------------------------------------------" << endl;
      assemble_rhs();
      pcout << "xxxx--------------------------------------------------------------------------------" << endl;
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

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
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs[0] );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );
    
    m_N = Particle_Number(m_Psi_ref);
    
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
    double T, N, W, V2;

    make_grid_custom();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    //m_Psi_2 = 0;

    unsigned bla[] ={0,1,0};
    CEigenfunctions<dim> Ef2( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef2, m_Psi_2 );
    m_Psi_2 *= 1.0/sqrt(Particle_Number(m_Psi_2));
      
    do_superposition();
    compute_U(m_Psi_ref,m_Psi_ref, m_U);
    m_workspace_1=m_U;
      
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "U");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu", mpi_communicator );      
      
    return;
      
    compute_E_lin( m_Psi_1, T, N, W, V2 );
    output_guess();

    m_mu = T/N+m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    m_mu = 0;
    
    m_ti = 1; //sqrt((m_mu*N-T+V2)/(m_gs[0]*W)); // WARNING: NEHARI
    
    
    status = DoIter();
    
    pcout << "L2_norm of m_Psi_ref: " << m_N << endl;


    if( status == Status::SUCCESS )
    {
      output_results("","final");
      save(  path + "final.bin" );
      Interpolate_R_to_C( path + "C-final.bin" );
      
      //dump_info_xml();
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
    double T, N, W, V2;
    int status;

    make_grid_custom();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W, V2 );
    
    double m_mu_0 = (T+V2)/N;
    m_mu = 0; //ceil(10.0*m_mu_0)/10.0 + m_gs[0]/fabs(m_gs[0])*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "V2 = " << V2 << endl;
    pcout << "m_mu = " << m_mu << endl;

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
      m_ti = sqrt(-T/V2);
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
        save_one( path + "final_one.bin" );
        vector<double> newgs = {m_gs[0]*m_N}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        m_ph.Set_Physics( "gs_1", m_gs );
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
      m_mu += m_gs[0]/fabs(m_gs[0])*m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W, V2 ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
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

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    constraints.distribute(m_Psi_ref);
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi_ref);

    triangulation.save( filename.c_str() );
  }
  
  template<int dim>
  void MySolver<dim>::save_one( string filename )
  {
    double tmp = Particle_Number(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    m_workspace_1*=sqrt(1/tmp);
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }

  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W, double& V2 )
  {
    m_computing_timer.enter_section(__func__);
    
    do_superposition();
    compute_U(m_Psi_1,m_Psi_1,m_U);
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    m_workspace_2 = m_U;
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<double> vec_vals_U(n_q_points);
    vector<Tensor<1, dim> > vec_grad(n_q_points);

    double JxW, vec_val_q;
    
    double tmp1[4]={}, res[4]={};
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vec_vals);
        fe_values.get_function_values(m_workspace_2, vec_vals_U);
        fe_values.get_function_gradients(m_workspace_1, vec_grad);
        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          vec_val_q = vec_vals[qp]*vec_vals[qp];
          tmp1[0] += JxW*( vec_grad[qp]*vec_grad[qp] + (Potential.value(fe_values.quadrature_point(qp)) + vec_vals_U[qp])*vec_val_q );
          tmp1[1] += JxW*vec_val_q;
          tmp1[2] += JxW*vec_val_q*vec_val_q;
          tmp1[3] += JxW*vec_vals_U[qp]*vec_val_q;
        }
      }
    }
    MPI_Allreduce( tmp1, res, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    T=res[0]; N=res[1]; W=res[2]; V2=res[3];
    m_computing_timer.exit_section();
  }

  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0, retval=0;
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );      

        for ( unsigned qp=0; qp<n_q_points; qp++ )
          tmp1 += fe_values.JxW(qp)*vec_vals[qp]*vec_vals[qp];
      }
    }
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }
  
  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);

    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    constraints.distribute(m_U);
    m_workspace_2=m_U;
   
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
    vector<double> U(n_q_points);

    double JxW, Q2, tmp;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi_ref);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_gradients(m_workspace_1, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          Q2 = Potential.value(fe_values.quadrature_point(qp)) + U[qp] - m_mu + 3.0*m_gs[0]*Psi_ref[qp]*Psi_ref[qp];

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

  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    do_superposition();
    do_superposition_U();
    
    m_workspace_1=m_Psi_ref;
    m_workspace_2=m_U;

    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
    vector<double> U(n_q_points);

    double JxW, Q1, pq;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi_ref);
        fe_values.get_function_values(m_workspace_2, U);
        fe_values.get_function_gradients(m_workspace_1, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pq = Psi_ref[qp]*Psi_ref[qp];
          Q1 = Potential.value(fe_values.quadrature_point(qp)) + U[qp] - m_mu +  m_gs[0]*pq;

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

  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    m_computing_timer.enter_section(__func__);
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    compute_U(m_Psi_ref,m_Psi_ref,m_U);
    m_workspace_2=m_U;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<double> U(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);

    double JxW, Q1, pq, tmp;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) + U[qp] - m_mu + m_gs[0]*(vals[qp]*vals[qp]);

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
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD));     
  
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::compute_U ( LA::MPI::Vector& vec1, LA::MPI::Vector& vec2, LA::MPI::Vector& U )
  {
    m_computing_timer.enter_section(__func__);
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    constraints.distribute(vec1);
    m_workspace_1=vec1;

    constraints.distribute(vec2);
    m_workspace_2=vec2;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);

    double JxW;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
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
          JxW = fe_values.JxW(qp);

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

    m_computing_timer.exit_section();
  }
 
  template <int dim>
  void MySolver<dim>::solve ()
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

  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential<dim> Potential_fct ( m_omega );
    
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
    //triangulation.refine_global(6);    
    double isovalues[] = {32,30,28};
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
          if( Potential_fct.value(p)  < isovalues[step] )
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
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver<DIMENSION> solver("params.xml");
    solver.run2b ();
  }
return EXIT_SUCCESS;
}
