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

/** 
 * Purpose: Real time propagation for the Gross-Pitaevskii equation (cartesian coordinates)
 * Method: fully implicit Crank-Nicolson 
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
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "global.h"
#include "mpi.h"
#include "my_table.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "MyComplexTools.h"
#include "muParser.h"

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string );
    ~MySolver();

    void run ();
    
  protected:
    void make_grid();
    void setup_system( const bool );
    void assemble_system();
    void assemble_rhs();
    
    void DoIter();
    void solve();
    void output_results ( string );
    void load( string );
    void save( string );    

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector newton_update;
    LA::MPI::Vector m_Psi; // Psi(t)
    LA::MPI::Vector m_Psi_t; // Psi(t)
    LA::MPI::Vector m_workspace;
    Vector<double> m_error_per_cell;

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_gs;
    double m_t;
    double m_dt;
    double m_dth; // 0.5*m_dt
    vector<double> m_omega;
    double m_res;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    unsigned m_NA;
    unsigned m_NK;
    
    MyTable m_table;    
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times )
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      m_zmin = m_ph.Get_Mesh("zrange",0);
      m_zmax = m_ph.Get_Mesh("zrange",1);

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }
    m_dth = 0.5*m_dt;
    m_t = 0;
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "shared_rt_prop.h"

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim> Potential ( m_omega );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    const double fak2 = 0.5*m_dt;
    const double fak4 = 0.25*m_gs*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double pot = Potential.value(fe_values.quadrature_point(qp)); 
  
          double sum_re = Psi[qp][0]+Psi_t[qp][0];
          double sum_im = Psi[qp][1]+Psi_t[qp][1];
          double sum_req = sum_re*sum_re;
          double sum_imq = sum_im*sum_im;
          double tmp1a = fak8*(sum_req + 3*sum_imq);
          double tmp1b = fak8*(sum_imq + 3*sum_req);
          double tmp2 = fak4*sum_re*sum_im;
  
          for( unsigned i=0; i<dofs_per_cell; i++ )
          {
            for( unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*(1.0-tmp2)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*(1.0+tmp2)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp));
              cell_matrix(i,j) += JxW*tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp));
            }
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix );
      }
    }
    system_matrix.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim> Potential ( m_omega );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    double JxW, pot=0, tmp1, sum_re, sum_im;

    const double fak2 = 0.5*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = Potential.value(fe_values.quadrature_point(qp)); 

          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];	  
          tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);	  
  
          for (unsigned int i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) -= JxW*fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp));
            cell_rhs(i) += JxW*fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp));
            cell_rhs(i) += JxW*((Psi_t[qp][0]-Psi[qp][0])*fe_values[rt].value(i,qp) - tmp1*sum_im*fe_values[rt].value(i,qp) + 
                                (Psi_t[qp][1]-Psi[qp][1])*fe_values[it].value(i,qp) + tmp1*sum_re*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global( cell_rhs, local_dof_indices, system_rhs );
      }
    }
    system_rhs.compress(VectorOperation::add);
    m_res = system_rhs.l2_norm();
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);

    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);

    constraints.distribute (newton_update);
    m_computing_timer.exit_section();
  }
  
  template<int dim> 
  void MySolver<dim>::DoIter()
  {
    m_computing_timer.enter_section(__func__);

    m_Psi_t = m_Psi;
    m_res = 0;
    assemble_rhs();
    pcout << "m_res = " << m_res << endl;       
    do
    {
      assemble_system();
      solve();
      
      system_rhs = m_Psi_t;
      system_rhs.add( -1, newton_update );
      constraints.distribute(system_rhs);
      m_Psi_t = system_rhs;   

      assemble_rhs();
      pcout << "m_res = " << m_res << endl;       
    }
    while( m_res >  1e-16 ); 
    m_t += m_dt;

    m_Psi = m_Psi_t;
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    load( "Cfinal.bin" );
    
    std::vector<double> p(dim);
    std::vector<double> pos(dim);
    std::vector<double> var(dim);

    output_results("");

    double N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
    pcout << "N == " << N << endl;
    
    MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, pos );
    MyComplexTools::MPI::Expectation_value_width( mpi_communicator, dof_handler, fe, m_Psi, pos, var );
    MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, p );

    pcout << "t == " << m_t << endl;
    pcout << "N == " << N << endl;
    pcout << "p == " << p[0]/N << ", " << p[1]/N << ", " << p[2]/N << endl;
    pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
    pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

    for( unsigned i=1; i<=m_NA; i++ )
    {
      for( unsigned j=1; j<=m_NK; j++ )
      {
        pcout << "t == " << m_t << endl;
        DoIter();
      }
      
      N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
      MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, pos );
      MyComplexTools::MPI::Expectation_value_width( mpi_communicator, dof_handler, fe, m_Psi, pos, var );
      MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, p );

      pcout << "N == " << N << endl;
      pcout << "p == " << p[0]/N << ", " << p[1]/N << ", " << p[2]/N << endl;
      pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
      pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

      output_results("");
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

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
      case 2: { realtime_propagation::MySolver<2> solver("params.xml");
                solver.run();
                break; }
      case 3: { realtime_propagation::MySolver<3> solver("params.xml");
                solver.run();
                break; }
      default: cout << "You have found a new dimension!" << endl;
    }    
  }
return EXIT_SUCCESS;
}
