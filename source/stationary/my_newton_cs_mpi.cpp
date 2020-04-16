//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
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
#include <unistd.h>

#include "global.h"
#include "mpi.h"
#include "functions_cs.h"
#include "my_table.h"
#include "ref_pt_list.h"
#include "MyParameterHandler.h"

#define _FAKTOR_ 6.28318530718
//#define _FAKTOR_ 1

namespace BreedSolver
{
  template <int dim>
  class MySolver;

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

    double m_I[8];
  protected:
    int DoIter( string="" );
    
    void make_grid_custom();
    void setup_system( const bool );
    void assemble_rhs();
    void assemble_system();
    void do_superposition();
    void save( string );
    
    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_contributions();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );
    void estimate_error( double& );
    void Interpolate_R_to_C( string );

    void output_results ( string, string = "step" );
    void output_guess ();

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    IndexSet locally_owned_dofs, locally_owned_dofs_2;
    IndexSet locally_relevant_dofs, locally_relevant_dofs_2;
    AffineConstraints<double> constraints, constraints_2;

    LA::MPI::SparseMatrix m_system_matrix, m_system_matrix_2;
    LA::MPI::Vector m_system_rhs, m_system_rhs_2;
    LA::MPI::Vector m_Psi_C;
    LA::MPI::Vector m_Psi_C_ghosted;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_newton_update;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    MyTable m_table;
    MyTable m_results;

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
  
  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    CPotential Potential( m_omega, m_QN1[2] );
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
          vec_val_q = vec_vals[q_point]*vec_vals[q_point];
          JxW = fabs(fe_values.quadrature_point(q_point)[1])*fe_values.JxW(q_point);
          T1 += JxW*( vec_grad[q_point]*vec_grad[q_point] + Potential.value(fe_values.quadrature_point(q_point))*vec_val_q );
          N1 += JxW*vec_val_q;
          W1 += JxW*vec_val_q*vec_val_q;
        }
      }
    }
    
    T1 *= _FAKTOR_;
    N1 *= _FAKTOR_;
    W1 *= _FAKTOR_;
    
    MPI_Allreduce( &T1, &T, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &N1, &N, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce( &W1, &W, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
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
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals );      

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          tmp1 += fabs(fe_values.quadrature_point(q_point)[1])*fe_values.JxW(q_point)*vec_vals[q_point]*vec_vals[q_point];
      }
    }
    tmp1 *= _FAKTOR_;
    
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }
  
  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    m_computing_timer.enter_section(__func__);

    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    //const QGauss<dim-1> face_quadrature_formula(fe.degree+1);
   
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);
    //FEFaceValues<dim> fe_face_values ( fe, face_quadrature_formula, update_gradients|update_values|update_quadrature_points|update_normal_vectors|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    //const unsigned int n_face_q_points = face_quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi_ref, Psi_ref);
        fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          const double pq = m_gs*Psi_ref[qp]*Psi_ref[qp];
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + pq;

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxW*(Psi_ref_grad[qp]*fe_values.shape_grad(i,qp) + Q1*Psi_ref[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp);
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
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD)); 
    m_computing_timer.exit_section();
  }
  
  template<int dim>
  void MySolver<dim>::Interpolate_R_to_C ( string filename )
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

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double JxW = fe_values_2.JxW(qp);
          const double JxWxV = JxW*vals[qp];
          
          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxWxV*fe_values_2[rt].value(i,qp);
            for( unsigned j=0; j<dofs_per_cell; ++j )
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
    solution_transfer.prepare_for_serialization(m_Psi_C_ghosted);
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
  void MySolver<dim>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    double JxW, Q1;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi_ref, Psi_ref);
        fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*Psi_ref[qp]*Psi_ref[qp];

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
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    double JxW, Q2;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi_ref, Psi_ref);
        fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          Q2 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + 3.0*m_gs*Psi_ref[qp]*Psi_ref[qp];

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
    Point<dim,double> pt1(0,0); 
    Point<dim,double> pt2(m_xmax,m_ymax);
    
    CPotential Potential_fct ( m_omega, m_QN1[2] );
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(6);
    
    //double isovalues[] = {70,60,50,40};
    //double isovalues[] = {60,50,40};
    double isovalues[] = {60,50};

    for( unsigned int step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
	        if( fabs(Potential_fct.value(p))  < isovalues[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }
    
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        const Point<dim> face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary() && !(face_center[1]==0) )    
          cell->face(f)->set_all_boundary_ids(1);
      }
    }

    unsigned int tmp1[2], tmp2[2];
    tmp1[0] = triangulation.n_cells();
    tmp1[1] = triangulation.n_active_cells();

    MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);
/*
    GridOutFlags::Msh opt(true, true);
    string filename = "grid-" + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4);
    ofstream output ((filename + ".msh").c_str());
    GridOut grid_out;
    grid_out.set_flags(opt);
    grid_out.write_msh (triangulation,output);
*/
    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    m_computing_timer.exit_section();
  }
  
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
    
    CPotential Potential_fct ( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+2);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim> > Psi_1_grad(n_q_points);
    vector<Tensor<1, dim> > Psi_2_grad(n_q_points);
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
        fe_values.get_function_gradients(m_workspace_1, Psi_1_grad);
        fe_values.get_function_gradients(m_workspace_2, Psi_2_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
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

    for( int i=0; i<5; ++i ) locint[i] *= m_gs;
    for( int i=0; i<8; ++i ) locint[i] *= _FAKTOR_;
    
    MPI_Allreduce( locint, m_I, 8, MPI_DOUBLE, MPI_SUM, mpi_communicator);
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

      m_Psi_1.reinit (locally_owned_dofs, mpi_communicator);
      m_Psi_2.reinit (locally_owned_dofs, mpi_communicator);
    }

    m_Psi_ref.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_newton_update.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    int myrank;
    MPI_Comm_rank( mpi_communicator, &myrank );    
    cout << "(" << myrank << ") locally_owned_dofs = " << m_Psi_1.local_size()  << endl;
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), constraints);
    if( m_QN1[2] > 0 )
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
    VectorTools::interpolate_boundary_values (dof_handler_2, 1, ZeroFunction<dim>(2), constraints_2);
    if( m_QN1[2] > 0 )
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
    
    string filename, filename2;
    
    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);
    m_workspace_1=m_Psi_1;
    m_workspace_2=m_Psi_2;
    
    CPotential Potential_fct ( m_omega, m_QN1[2] );
    VectorTools::interpolate (dof_handler, Potential_fct, m_workspace_ng );
    constraints.distribute(m_workspace_ng);
    m_Psi_ref=m_workspace_ng;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi_1");
    data_out.add_data_vector (m_workspace_2, "m_Psi_2");
    data_out.add_data_vector (m_Psi_ref, "m_Potential");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu",mpi_communicator);

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    m_computing_timer.enter_section(__func__);

    string filename, filename2;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_ref, "Psi_sol");
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
    
    m_res=0;
    m_res_old=m_res;
    m_counter=0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      m_Psi_2.add( -1e-4*m_t[1]/fabs(m_t[1]), m_newton_update); 
      constraints.distribute(m_Psi_2);

      this->find_ortho_min();

      do_superposition();
      assemble_rhs();
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( this->l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

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
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );    

    do_superposition();
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
    double T, N, W;

    make_grid_custom();
    setup_system(true);

    CEigenfunctions Ef1( m_QN1, m_omega );
    //CEigenfunctions Ef2( m_QN2, f );

    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
      
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    status = DoIter();
    
    pcout << "L2_norm of m_Psi_ref: " << m_N << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","final");
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
    setup_system(true);

    CEigenfunctions Ef1( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    //m_Psi_2 = m_Psi_1;
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;
    
    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; ++i )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;
      
      // nehari
      //m_ti = sqrt((m_mu*N-T)/(4.0*m_gs*W)); // if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu*N-T)/(m_gs*W));

      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "m_mu = " << m_mu << endl;      
      pcout << "m_ti = " << m_ti << endl;      
      
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
        output_results(path,"final");
        save( path + "final.bin" );
        this->dump_info_xml(path);
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
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
      m_mu += m_gs/fabs(m_gs)*m_dmu;
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }
  
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    constraints.distribute(m_Psi_ref);    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_Psi_ref);

    triangulation.save( filename.c_str() );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run2b ();
  }
return EXIT_SUCCESS;
}
