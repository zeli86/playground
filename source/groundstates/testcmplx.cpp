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
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
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

using namespace std;
using namespace dealii;

class CTestfunction : public Function<2>
{
  public:
    CTestfunction () : Function<2>(2) { 
    }

    virtual double value ( const Point<2> &p, const unsigned int  component = 0) const;
};

double CTestfunction::value( const Point<2> &p, const unsigned int component ) const
{
  double retval = 0;

  switch( component )
  {
    case 0: retval = p[0]+p[1]; break;
    case 1: retval = 0; break;
  }
return retval;
}

class CTestfunction2 : public Function<2>
{
  public:
    CTestfunction2 () : Function<2>(2) { 
    }

    virtual double value ( const Point<2> &p, const unsigned int  component = 0) const;
};

double CTestfunction2::value( const Point<2> &p, const unsigned int component ) const
{
  double retval = 0;

  switch( component )
  {
    case 0: retval = 0; break;
    case 1: retval = 1; break;
  }
return retval;
}

class MySolver
{
public:
  MySolver();
  ~MySolver();

  void run ();
    
protected:
  void make_grid();
  void setup_system();
  void assemble();

  MPI_Comm mpi_communicator;
  parallel::distributed::Triangulation<2> triangulation;
  FESystem<2> fe;
  DoFHandler<2> dof_handler;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  ConstraintMatrix constraints;

  LA::MPI::SparseMatrix system_matrix;
  LA::MPI::Vector m_wf_1; 
  LA::MPI::Vector m_wf_2; 
  LA::MPI::Vector m_workspace; 

  ConditionalOStream pcout;

  bool m_root;
};

  MySolver::MySolver () 
    : 
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<2>::MeshSmoothing(Triangulation<2>::smoothing_on_refinement|Triangulation<2>::smoothing_on_coarsening)),
    fe (FE_Q<2>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {    
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  }

  void MySolver::make_grid ()
  {
    Point<2,double> pt1( -3, -3 );
    Point<2,double> pt2( 3, 3 );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(7);
  }  

  void MySolver::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_wf_1.reinit (locally_owned_dofs, mpi_communicator);
    m_wf_2.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  }  

  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  void MySolver::assemble ()
  {
    const QGauss<2> quadrature_formula(fe.degree+1);

    CTestfunction func1;
    CTestfunction2 func2;

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix = 0;
    m_wf_1 = 0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double a = func1.value(fe_values.quadrature_point(qp),0);
          double b = func1.value(fe_values.quadrature_point(qp),1);
          double c = func2.value(fe_values.quadrature_point(qp),0);
          double d = func2.value(fe_values.quadrature_point(qp),1);

          for (unsigned i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));
            }
            cell_rhs(i) += JxW*((a*c-b*d)*fe_values[rt].value(i,qp) + (b*c+a*d)*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, m_wf_1);
    }
  }
  system_matrix.compress(VectorOperation::add);
  m_wf_1.compress(VectorOperation::add);
}
  
void MySolver::run()
{
  make_grid();
  setup_system();
  assemble();
  //system_matrix.vmult( m_wf_2, m_wf_1 );

  SolverControl solver_control (m_wf_1.size(), 1e-15);
    
  m_wf_2=0;
  PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
  PETScWrappers::PreconditionSOR preconditioner(system_matrix);
  solver.solve (system_matrix, m_wf_2, m_wf_1, preconditioner);

  m_workspace = m_wf_1;

  vector<std::string> solution_names;
  solution_names.push_back ("Re_Psi");
  solution_names.push_back ("Im_Psi");    

  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (m_workspace, solution_names );
  data_out.build_patches ();
  data_out.write_vtu_in_parallel ( "cmplxdemo.vtu", mpi_communicator );  
}

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MySolver solver;
    solver.run();
  }
return EXIT_SUCCESS;
}
