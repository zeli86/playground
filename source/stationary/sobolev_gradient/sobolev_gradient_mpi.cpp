/*
 atus-pro testing - atus-pro testing playgroung
 Copyright (C) 2022 Želimir Marojević <zelimir.marojevic@gmail.com>

 This file is part of atus-pro testing.

 atus-pro testing is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 atus-pro testing is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "mpi.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "MyLogging.h"

#include "sobolev_gradient_mpi.hpp"

namespace BreedSolver
{

  template <int dim>
  MySolver<dim>::MySolver(const std::string& xmlfilename)
    :
    m_oParameters(xmlfilename),
    mpi_communicator(MPI_COMM_WORLD),
    m_oTriangulation(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    m_oFe(2),
    m_oDofHandler(m_oTriangulation)
  {
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  }


  template <int dim>
  MySolver<dim>::~MySolver()
  {
    m_oDofHandler.clear();
  }


  template<int dim>
  void MySolver<dim>::project_gradient()
  {
    double tmp1[] = {0, 0}, sum[] = {0, 0};

    m_oConstraints.distribute(m_vPhiSobolev);
    m_vWorkspace1 = m_vPhiSobolev;

    m_oConstraints.distribute(m_vSobGradient);
    m_vWorkspace2 = m_vSobGradient;

    m_oConstraints.distribute(m_vPhi);
    m_vWorkspace3 = m_vPhi;

    const QGauss<dim>  quadrature_formula(m_oFe.degree + 1);
    FEValues<dim> fe_values(m_oFe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    std::vector<double> vals_1(n_q_points);
    std::vector<double> vals_2(n_q_points);
    std::vector<double> vals_3(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(m_vWorkspace1, vals_1);
        fe_values.get_function_values(m_vWorkspace2, vals_2);
        fe_values.get_function_values(m_vWorkspace3, vals_3);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          tmp1[0] += JxW * (vals_3[qp] * vals_2[qp]);
          tmp1[1] += JxW * (vals_3[qp] * vals_1[qp]);
        }
      }
    }
    MPI_Allreduce(tmp1, sum, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_vSobGradient.add(-sum[0] / sum[1], m_vPhiSobolev);
  }


  template <int dim>
  void MySolver<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(m_oFe.degree + 1);

    m_oMatrix = 0;

    FEValues<dim> fe_values(m_oFe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = m_oFe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);

          for (unsigned i = 0; i < dofs_per_cell; ++i)
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        m_oConstraints.distribute_local_to_global(cell_matrix, local_dof_indices, m_oMatrix);
      }
    }
    m_oMatrix.compress(VectorOperation::add);
  }


  template <int dim>
  void MySolver<dim>::compute_psi_sobolev()
  {
    const QGauss<dim> quadrature_formula(m_oFe.degree + 1);

    m_oMatrix = 0;
    m_vRhs = 0;

    m_oConstraints.distribute(m_vPhi);
    m_vWorkspace1 = m_vPhi;

    FEValues<dim> fe_values(m_oFe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = m_oFe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<double> vals(n_q_points);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(m_vWorkspace1, vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double val = vals[qp];
          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) = JxW * val * fe_values.shape_value(i, qp);
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
            }
          }
        }
        cell->get_dof_indices(local_dof_indices);
        m_oConstraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_oMatrix, m_vRhs);
      }
    }
    m_oMatrix.compress(VectorOperation::add);
    m_vRhs.compress(VectorOperation::add);

    BOOST_LOG_TRIVIAL(info) << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_oMatrix, m_vPhiSobolev, m_vRhs);
    m_oConstraints.distribute(m_vPhiSobolev);

    /*
        m_Psi_sob=0;
        SolverControl solver_control ( m_Psi_sob.size(), 1e-15 );
        PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

        PETScWrappers::PreconditionSOR preconditioner;
        PETScWrappers::PreconditionSOR::AdditionalData data;
        preconditioner.initialize(m_system_matrix, data);

        cg.solve (m_system_matrix, m_Psi_sob, m_system_rhs, preconditioner);
    */
  }


  template <int dim>
  void MySolver<dim>::solve()
  {
    BOOST_LOG_TRIVIAL(info) << "Solving..." << endl;

    SolverControl solver_control;

    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_oMatrix, m_vSobGradient, m_vRhs);
    m_oConstraints.distribute(m_vSobGradient);
  }

  /*
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_sob_grad=0;
    SolverControl solver_control ( m_sob_grad.size(), 1e-15 );
    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(m_system_matrix, data);
    cg.solve (m_system_matrix, m_sob_grad, m_system_rhs, preconditioner);
    constraints.distribute (m_sob_grad);
  }
  */

  template <int dim>
  void MySolver<dim>::make_grid()
  {
    Point<dim, double> pt1;
    Point<dim, double> pt2;

    std::vector<double> oGridCornerOne{0, 0, 0};
    std::vector<double> oGridCornerTwo{0, 0, 0};

    m_oParameters.get("grid.grid_corner_one", oGridCornerOne);
    m_oParameters.get("grid.grid_corber_two", oGridCornerTwo);

    for (int i = 0; i < dim; i++)
    {
      pt1(i) = oGridCornerOne[i];
      pt2(i) = oGridCornerTwo[i];
    }

    GridGenerator::hyper_rectangle(m_oTriangulation, pt2, pt1);
    m_oTriangulation.refine_global(m_iGlobalRefinement);
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_oDofHandler.distribute_dofs(m_oFe);

    m_oLocallyOwnedDofs = m_oDofHandler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(m_oDofHandler, m_oLocallyRelevantDofs);
    m_vPhi.reinit(m_oLocallyOwnedDofs, mpi_communicator);
    m_vPhiSobolev.reinit(m_oLocallyOwnedDofs, mpi_communicator);
    m_vRhs.reinit(m_oLocallyOwnedDofs, mpi_communicator);
    m_vSobGradient.reinit(m_oLocallyOwnedDofs, mpi_communicator);
    m_vWorkspace1.reinit(m_oLocallyOwnedDofs, m_oLocallyRelevantDofs, mpi_communicator);
    m_vWorkspace2.reinit(m_oLocallyOwnedDofs, m_oLocallyRelevantDofs, mpi_communicator);
    m_vWorkspace3.reinit(m_oLocallyOwnedDofs, m_oLocallyRelevantDofs, mpi_communicator);

    cout << "(" << m_rank << ") locally_owned_dofs = " << m_vPhi.local_size()  << endl;

    m_oConstraints.clear();
    m_oConstraints.reinit(m_oLocallyRelevantDofs);
    dealii::DoFTools::make_hanging_node_constraints(m_oDofHandler, m_oConstraints);
    dealii::VectorTools::interpolate_boundary_values(m_oDofHandler, 0, ZeroFunction<dim>(), m_oConstraints);
    m_oConstraints.close();

    dealii::DynamicSparsityPattern csp(m_oLocallyRelevantDofs);
    dealii::DoFTools::make_sparsity_pattern(m_oDofHandler, csp, m_oConstraints, false);
    dealii::SparsityTools::distribute_sparsity_pattern(csp, m_oDofHandler.locally_owned_dofs(), mpi_communicator, m_oLocallyRelevantDofs);
    m_oMatrix.reinit(m_oLocallyOwnedDofs, m_oLocallyOwnedDofs, csp, mpi_communicator);
  }


  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    using namespace utils::real_wavefunction;

    int retval = Status::SUCCESS;

    CPotential<dim> Potential(m_omega);

    m_oConstraints.distribute(m_vPhi);
    m_vWorkspace1 = m_vPhi;
    assemble_L2gradient(dynamic_cast<IRealWavefunction<dim>*>(this), m_vWorkspace1, Potential, m_rMu, m_rG, m_vRhs);

    for (int iCounter = 0; iCounter < m_iMaxIter; ++iCounter)
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "- " << iCounter;

      assemble_system();
      solve();

      compute_psi_sobolev();
      project_gradient();
      m_res = m_vSobGradient.l2_norm();

      m_vPhi.add(-1e-3, m_vSobGradient);

      if (fabs(m_N - 1) > 1e-5)
      {
        m_vPhi *= 1 / sqrt(m_N);
      }

      m_oConstraints.distribute(m_vPhi);
      m_vWorkspace1 = m_vPhi;
      assemble_L2gradient(dynamic_cast<IRealWavefunction<dim>*>(this), m_vWorkspace1, Potential, m_rMu, m_rG, m_vRhs);

      m_N = particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), m_vWorkspace1, mpi_communicator);
      if (iCounter % m_NA == 0)
      {
        output_results(path);
      }

      if (m_res < m_rEpsilon)
      {
        retval = Status::SUCCESS;
        break;
      }
    }
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    using namespace utils::real_wavefunction;

    make_grid();
    setup_system();

    unsigned QN[] = {0, 0, 0};
    CEigenfunctions<dim> Ef1(QN, m_omega);
    VectorTools::interpolate(m_oDofHandler, Ef1, m_vPhi);
    CPotential<dim> Potential(m_omega);

    m_vWorkspace1 = m_vPhi;
    m_vPhi *= sqrt(1.0 / particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), m_vWorkspace1, mpi_communicator));

    auto status = DoIter("");

    if (status == Status::SUCCESS)
    {
      save("final.bin");
    }
  }


  template <int dim>
  void MySolver<dim>::output_results(string path, string prefix)
  {
    string filename;

    m_oConstraints.distribute(m_vPhi);
    m_vWorkspace1 = m_vPhi;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(m_oDofHandler);
    data_out.add_data_vector(m_vWorkspace1, "phi");
    data_out.build_patches();

    //filename = path + prefix + "-" + Utilities::int_to_string(m_counter, 5) + ".vtu";
    data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
  }


  template<int dim>
  void MySolver<dim>::save(string filename)
  {
    m_vWorkspace1 = m_vPhi;
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(m_oDofHandler);
    solution_transfer.prepare_for_serialization(m_vWorkspace1);
    m_oTriangulation.save(filename.c_str());
  }

  template class MySolver<2>;
  template class MySolver<3>;
} // end of namespace
