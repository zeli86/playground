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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "global.h"
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
    mpi_communicator(MPI_COMM_WORLD),
    m_oParameters(xmlfilename),
    m_root(Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    m_oTriangulation(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    fe(gl_degree_fe),
    dof_handler(m_oTriangulation)
  {
    try
    {
      m_oParameters.get("physics.omega", m_omega);
      m_oParameters.get("physics.gs_1", m_gs);

      // m_NA = int(m_ph.Get_Algorithm("NA", 0));
      // m_epsilon = m_ph.Get_Algorithm("epsilon");
      // m_guess_str = m_ph.Get_Parameter("guess_fct");
    }
    catch (const std::string& info)
    {
      std::cerr << info << endl;
      MPI_Abort(mpi_communicator, 0);
    }
  }

  template <int dim>
  MySolver<dim>::~MySolver()
  {
    dof_handler.clear();
  }


  template<int dim>
  void MySolver<dim>::Project_gradient()
  {
    double tmp1[] = {0, 0}, sum[] = {0, 0};

    constraints.distribute(m_Psi_sob);
    m_workspace_1 = m_Psi_sob;

    constraints.distribute(m_sob_grad);
    m_workspace_2 = m_sob_grad;

    constraints.distribute(m_Psi);
    m_workspace_3 = m_Psi;

    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);
    vector<double> vals_3(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(m_workspace_1, vals_1);
        fe_values.get_function_values(m_workspace_2, vals_2);
        fe_values.get_function_values(m_workspace_3, vals_3);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          tmp1[0] += JxW * (vals_3[qp] * vals_2[qp]);
          tmp1[1] += JxW * (vals_3[qp] * vals_1[qp]);
        }
      }
    }
    MPI_Allreduce(tmp1, sum, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);

    m_sob_grad.add(-sum[0] / sum[1], m_Psi_sob);
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;

    m_system_rhs = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    Vector<double> cell_rhs(dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1, dim>> grads(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) + m_rG * (vals[qp] * vals[qp]);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW * (grads[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);
  }


  template <int dim>
  void MySolver<dim>::assemble_system()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, m_system_matrix);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
  }


  template <int dim>
  void MySolver<dim>::compute_Psi_sob()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;
    m_system_rhs = 0;

    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<double> vals(n_q_points);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(m_workspace_1, vals);

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
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
    m_system_rhs.compress(VectorOperation::add);

    BOOST_LOG_TRIVIAL(info) << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_Psi_sob, m_system_rhs);
    constraints.distribute(m_Psi_sob);

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
    solver.solve(m_system_matrix, m_sob_grad, m_system_rhs);
    constraints.distribute(m_sob_grad);
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
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    m_Psi.reinit(locally_owned_dofs, mpi_communicator);
    m_Psi_sob.reinit(locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_sob_grad.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_3.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(m_oTriangulation.n_active_cells());

    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi.local_size()  << endl;

    vector<bool> mask(dof_handler.get_fe().n_components(), true);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::output_guess()
  {
    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(m_workspace_1, "m_Psi");
    data_out.build_patches();
    data_out.write_vtu_in_parallel("guess.vtu", mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::output_results(string path, string prefix)
  {
    string filename;

    Vector<float> subdomain(m_oTriangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
    {
      subdomain(i) = m_oTriangulation.locally_owned_subdomain();
    }

    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(m_workspace_1, "Psi");
    data_out.add_data_vector(m_error_per_cell, "error per cell");
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();

    filename = path + prefix + "-" + Utilities::int_to_string(m_counter, 5) + ".vtu";
    data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
  }

  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    int retval = Status::SUCCESS;

    assemble_rhs();

    m_res = 0;
    m_counter = 0;
    do
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "- " << m_counter;

      assemble_system();
      solve();

      compute_Psi_sob();
      Project_gradient();
      m_res = m_sob_grad.l2_norm();

      m_Psi.add(-1e-3, m_sob_grad);

      //compute_mu(m_rMu);
      // m_N = MyRealTools::MPI::Particle_Number(mpi_communicator, dof_handler, fe, m_Psi);

      if (fabs(m_N - 1) > 1e-5)
      {
        m_Psi *= 1 / sqrt(m_N);
      }

      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        output_results(path);
      }

      // columns& cols = m_table.new_line();
      // m_table.insert(cols, MyTable::COUNTER, double(m_counter));
      // m_table.insert(cols, MyTable::RES, m_res);
      // m_table.insert(cols, MyTable::RESP, m_resp);
      // m_table.insert(cols, MyTable::MU, m_rMu);
      // m_table.insert(cols, MyTable::GS, m_gs);
      // m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      // if (m_root)
      // {
      //   cout << m_table;
      // }
      if (m_res < m_epsilon[1])
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while (true);

    //TODO: m_N = MyRealTools::MPI::Particle_Number(mpi_communicator, dof_handler, fe, m_Psi);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }

    return retval;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    int status;
    double T, N, W;

    make_grid();
    setup_system();

    unsigned QN[] = {0, 0, 0};
    CEigenfunctions<dim> Ef1(QN, m_omega);
    VectorTools::interpolate(dof_handler, Ef1, m_Psi);
    CPotential<dim> Potential(m_omega);


    this->m_workspace_1 = this->m_Psi[0];
    //MyRealTools::MPI::Compute_E_lin(mpi_communicator, dof_handler, fe, constraints, Potential, m_workspace_1, T, N, W);
    m_Psi *= sqrt(1 / N);

    output_guess();

    BOOST_LOG_TRIVIAL(info) << setprecision(9);
    BOOST_LOG_TRIVIAL(info) << "T = " << T << endl;
    BOOST_LOG_TRIVIAL(info) << "N = " << N << endl;
    BOOST_LOG_TRIVIAL(info) << "W = " << W << endl;

    status = DoIter("");

    if (status == Status::SUCCESS)
    {
      //output_results("","final");
      //output_results("","final");
      save("final.bin");
      //dump_info_xml("");
    }

    if (m_root)
    {
      ofstream ofs("log.txt");
      //ofs << m_table;
    }
  }

  template<int dim>
  void MySolver<dim>::save(string filename)
  {
    m_workspace_1 = m_Psi;
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_workspace_1);
    m_oTriangulation.save(filename.c_str());
  }

  template class MySolver<2>;
} // end of namespace

