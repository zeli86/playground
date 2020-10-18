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

#include "default_includes.h"

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
#include "functions.h"
#include "MyParameterHandler.h"
#include "my_table.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

#include "cBaseMPI.h"

  template <int dim>
  class MySolver : public cBaseMPI<dim, 2, dealii::FESystem<dim>>
  {
  public:
    explicit MySolver(const std::string);
    virtual ~MySolver() = default;

    void run2b();

  protected:

    void compute_contributions();

    int DoIter(string = "");

    void make_grid_custom();
    void assemble_rhs();
    void assemble_system();
    void save_one(string);
    void dump_info_xml(const string);

    double Particle_Number(LA::MPI::Vector&);

    void solve();
    void compute_E_lin(LA::MPI::Vector&, double&, double&, double&);
    void estimate_error(double&);

    FESystem<dim> m_FESys;

    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::mpi_communicator;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_root;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_rank;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_t;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_ti;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_N;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_omega;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_mu;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_gs;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_counter;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::pcout;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_computing_timer;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_ph;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_final_error;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_NA;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_Ndmu;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_dmu;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_QN1;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_res;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_resp;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_res_old;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_epsilon;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_maxiter;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_total_no_cells;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_total_no_active_cells;
    using cBaseMPI<dim, 2, dealii::FESystem<dim>>::m_global_refinement;

    string m_guess_str;
    MyTable m_table;
    MyTable m_results;
  };

  /*************************************************************************************************/
  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const string xmlfilename)
    :
    m_FESys(FE_Q<dim>(gl_degree_fe), 2),
    cBaseMPI<dim, 2, dealii::FESystem<dim>> (xmlfilename, m_FESys)
  {
  }

  template <int dim>
  void MySolver<dim>::compute_E_lin(LA::MPI::Vector& vec, double& T, double& N, double& W)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->m_constraints.distribute(vec);
    this->m_Workspace[0] = vec;

    CPotential<dim> Potential(m_omega);
    const QGauss<dim>  quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> grad(n_q_points, vector<Tensor<1, dim>>(2));

    double tmp[] = {0, 0, 0}, res[] = {0, 0, 0};

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vals);
        fe_values.get_function_gradients(this->m_Workspace[0], grad);
        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double vec_val_q = vals[qp] * vals[qp];
          const double JxW = fe_values.JxW(qp);
          tmp[0] += JxW * (grad[qp][0] * grad[qp][0] + grad[qp][1] * grad[qp][1] + Potential.value(fe_values.quadrature_point(qp)) * vec_val_q);
          tmp[1] += JxW * vec_val_q;
          tmp[2] += JxW * vec_val_q * vec_val_q;
        }
      }
    }

    MPI_Allreduce(tmp, res, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    T = res[0];
    N = res[1];
    W = res[2];

    
  }

  template<int dim>
  double MySolver<dim>::Particle_Number(LA::MPI::Vector& vec)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    double tmp1 = 0, retval = 0;

    this->m_constraints.distribute(vec);
    this->m_Workspace[0] = vec;

    const QGauss<dim>  quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points, Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          tmp1 += fe_values.JxW(qp) * (vals[qp] * vals[qp]);
        }
      }
    }

    MPI_Allreduce(&tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
    return retval;
  }

  template <int dim>
  void MySolver<dim>::estimate_error(double& err)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt(0);
    const FEValuesExtractors::Scalar it(1);

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_constraints.distribute(this->m_Psi_Ref);
    this->m_Workspace[0] = this->m_Psi_Ref;

    this->m_System_RHS = 0;
    this->m_System_Matrix = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Vector<double>> vals(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> grads(n_q_points, vector<Tensor<1, dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_rhs = 0;
        cell_matrix = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vals);
        fe_values.get_function_gradients(this->m_Workspace[0], grads);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs * (vals[qp] * vals[qp]);

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW * (grads[qp][0] * fe_values[rt].gradient(i, qp) + Q1 * vals[qp][0] * fe_values[rt].value(i, qp) + grads[qp][1] * fe_values[it].gradient(i, qp) + Q1 * vals[qp][1] * fe_values[it].value(i, qp));
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) + fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
            }
          }
        }
        cell->get_dof_indices(local_dof_indices);
        this->m_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, this->m_System_Matrix, this->m_System_RHS);
      }
    }
    this->m_System_RHS.compress(VectorOperation::add);
    this->m_System_Matrix.compress(VectorOperation::add);

    solve();

    this->m_Workspace[0] = this->m_Search_Direction;
    VectorTools::integrate_difference(this->m_DOF_Handler, this->m_Workspace[0], ZeroFunction<dim>(2), this->m_error_per_cell, QGauss<dim>(this->m_FE.degree + 2), VectorTools::L2_norm);
    const double total_local_error = this->m_error_per_cell.l2_norm();
    err = std::sqrt(Utilities::MPI::sum(total_local_error * total_local_error, MPI_COMM_WORLD));
    
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt(0);
    const FEValuesExtractors::Scalar it(1);

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_System_RHS = 0;

    this->m_constraints.distribute(this->m_Psi_Ref);
    this->m_Workspace[0] = this->m_Psi_Ref;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    Vector<double> cell_rhs(dofs_per_cell);
    vector<Vector<double>> vals(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> grads(n_q_points, vector<Tensor<1, dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vals);
        fe_values.get_function_gradients(this->m_Workspace[0], grads);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs * (vals[qp] * vals[qp]);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += JxW * (grads[qp][0] * fe_values[rt].gradient(i, qp) + Q1 * vals[qp][0] * fe_values[rt].value(i, qp)
                                  + grads[qp][1] * fe_values[it].gradient(i, qp) + Q1 * vals[qp][1] * fe_values[it].value(i, qp));
        }
        cell->get_dof_indices(local_dof_indices);
        this->m_constraints.distribute_local_to_global(cell_rhs, local_dof_indices, this->m_System_RHS);
      }
    }
    this->m_System_RHS.compress(VectorOperation::add);
    m_res = this->m_System_RHS.l2_norm();
    
  }

  template <int dim>
  void MySolver<dim>::assemble_system()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt(0);
    const FEValuesExtractors::Scalar it(1);

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);
    //const QGauss<dim-1> face_quadrature_formula(this->m_FE.degree+1);

    this->m_constraints.distribute(this->m_Psi_Ref);
    this->m_Workspace[0] = this->m_Psi_Ref;
    this->m_System_Matrix = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);
    //FEFaceValues<dim> fe_face_values ( fe, face_quadrature_formula, update_gradients|update_values|update_quadrature_points|update_normal_vectors|update_JxW_values);

    const unsigned int dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    //const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<vector<Tensor<1, dim>>> Psi_ref_grad(n_q_points, vector<Tensor<1, dim>>(2));
    vector<Vector<double>> Psi_ref(n_q_points, Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], Psi_ref);
        fe_values.get_function_gradients(this->m_Workspace[0], Psi_ref_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double fak = m_gs * Psi_ref[qp][0] * Psi_ref[qp][1];
          const double Pot = Potential.value(fe_values.quadrature_point(qp)) - m_mu;
          const double req = Psi_ref[qp][0] * Psi_ref[qp][0];
          const double imq = Psi_ref[qp][1] * Psi_ref[qp][1];

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) += JxW * (fe_values[rt].gradient(i, qp) * fe_values[rt].gradient(j, qp) + (Pot + m_gs * (3 * req + imq)) * fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp)
                                          + fak * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp)
                                          + fak * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp)
                                          + fe_values[it].gradient(i, qp) * fe_values[it].gradient(j, qp) + (Pot + m_gs * (3 * imq + req)) * fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
        }
        cell->get_dof_indices(local_dof_indices);
        this->m_constraints.distribute_local_to_global(cell_matrix, local_dof_indices, this->m_System_Matrix);
      }
    }
    this->m_System_Matrix.compress(VectorOperation::add);
    
  }


  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->update_workspace();

    CPotential<dim> Potential_fct(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();

    vector<Vector<double>> Psi_0(n_q_points, Vector<double>(2));
    vector<Vector<double>> Psi_1(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> Psi_0_grad(n_q_points, vector<Tensor<1, dim>>(2));
    vector<vector<Tensor<1, dim>>> Psi_1_grad(n_q_points, vector<Tensor<1, dim>>(2));

    std::array<double, 8> local_contributions;
    std::array<double, 8> total_contributions;
    local_contributions.fill(0);
    total_contributions.fill(0);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], Psi_0);
        fe_values.get_function_values(this->m_Workspace[1], Psi_1);
        fe_values.get_function_gradients(this->m_Workspace[0], Psi_0_grad);
        fe_values.get_function_gradients(this->m_Workspace[1], Psi_1_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double rtp0 = Psi_0[qp][0];
          const double imp0 = Psi_0[qp][1];
          const double rtp1 = Psi_1[qp][0];
          const double imp1 = Psi_1[qp][1];
          const double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu;

          local_contributions[0] += JxW * (pow(rtp0, 4) + 2 * pow(rtp0, 2) * pow(imp0, 2) + pow(imp0, 4));
          local_contributions[1] += JxW * (pow(rtp1, 4) + 2 * pow(rtp1, 2) * pow(imp1, 2) + pow(imp1, 4));
          // local_contributions[2] += JxW * p0q * p1q;
          // local_contributions[3] += JxW * p0q * p01;
          // local_contributions[4] += JxW * p1q * p01;
          local_contributions[5] += JxW * (Psi_0_grad[qp][0] * Psi_0_grad[qp][0] + Psi_0_grad[qp][1] * Psi_0_grad[qp][1] + Q * (Psi_0[qp][0] * Psi_0[qp][0] + Psi_0[qp][1] * Psi_0[qp][1]));
          local_contributions[6] += JxW * (Psi_1_grad[qp][0] * Psi_1_grad[qp][0] + Psi_1_grad[qp][1] * Psi_1_grad[qp][1] + Q * (Psi_1[qp][0] * Psi_1[qp][0] + Psi_1[qp][1] * Psi_1[qp][1]));
          // local_contributions[7] += JxW * (Psi_0_grad[qp] * Psi_1_grad[qp] + Q * p01);
        }
      }
    }

    MPI_Allreduce(local_contributions.data(), total_contributions.data(), total_contributions.size(), MPI_DOUBLE, MPI_SUM, mpi_communicator);

    this->m_coeffs["t0^2"]      = 0.5 * total_contributions[5];
    this->m_coeffs["t0_t1"]     = total_contributions[7];
    this->m_coeffs["t1^2"]      = 0.5 * total_contributions[6];
    this->m_coeffs["t0^4"]      = 0.25 * m_gs * total_contributions[0];
    this->m_coeffs["t0^1_t1^3"] = 0.25 * m_gs * 4 * total_contributions[4];
    this->m_coeffs["t0^2_t1^2"] = 0.25 * m_gs * 6 * total_contributions[2];
    this->m_coeffs["t0^3_t1^1"] = 0.25 * m_gs * 4 * total_contributions[3];
    this->m_coeffs["t1^4"]      = 0.25 * m_gs * total_contributions[1];

    
  }


  /*
  template <int dim>
  void MySolver<dim>::solve ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    pcout << "Solving..." << endl;

    SolverControl solver_control;

    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(this->m_System_Matrix, this->m_Search_Direction, this->m_System_RHS);
    this->m_constraints.distribute (this->m_Search_Direction);
    
  }
  */

  template <int dim>
  void MySolver<dim>::solve()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->m_Search_Direction = 0;

    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(this->m_System_Matrix, data);

    SolverControl solver_control(this->m_DOF_Handler.n_dofs(), 1e-15);
    PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);

    solver.solve(this->m_System_Matrix, this->m_Search_Direction, this->m_System_RHS, preconditioner);
    this->m_constraints.distribute(this->m_Search_Direction);

    
  }

  template <int dim>
  void MySolver<dim>::make_grid_custom()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    Point<dim, double> pt1;
    Point<dim, double> pt2;

    double min[] = {m_ph.Get_Mesh("xrange", 0), m_ph.Get_Mesh("yrange", 0), m_ph.Get_Mesh("zrange", 0)};
    double max[] = {m_ph.Get_Mesh("xrange", 1), m_ph.Get_Mesh("yrange", 1), m_ph.Get_Mesh("zrange", 1)};

    for (int i = 0; i < dim; i++)
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    CPotential<dim> Potential_fct(m_omega);

    GridGenerator::hyper_rectangle(this->m_Triangulation, pt2, pt1);
    this->m_Triangulation.refine_global(5);

    double isovalues[] = {15, 13, 11};

    for (unsigned int step = 0; step < sizeof(isovalues) / sizeof(double); step++)
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->m_Triangulation.begin_active(),
                                                                               endc = this->m_Triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          Point<dim> p = cell->vertex(v);
          if (fabs(Potential_fct.value(p))  < isovalues[step])
          {
            cell->set_refine_flag();
            break;
          }
        }
      this->m_Triangulation.execute_coarsening_and_refinement();
    }

    unsigned int tmp1[2], tmp2[2];
    tmp1[0] = this->m_Triangulation.n_cells();
    tmp1[1] = this->m_Triangulation.n_active_cells();

    MPI_Allreduce(tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    
  }


  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    int retval = Status::SUCCESS;

    m_table.clear();

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    this->do_linear_superposition();
    assemble_rhs();

    m_res = 0;
    m_res_old = m_res;
    m_counter = 0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      this->m_Psi[1].add(-1e-3 * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi[1]);

      this->find_ortho_min();

      this->do_linear_superposition();
      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        this->output_results(path);
      }

      columns& cols = m_table.new_line();
      m_table.insert(cols, MyTable::COUNTER, double(m_counter));
      m_table.insert(cols, MyTable::RES, m_res);
      m_table.insert(cols, MyTable::RESP, m_resp);
      m_table.insert(cols, MyTable::MU, m_mu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if (m_root)
      {
        cout << m_table;
      }
      if (m_res < m_epsilon[0])
      {
        retval = Status::SUCCESS;
        break;
      }
      if (this->l2norm_t() < 1e-4)
      {
        retval = Status::ZERO_SOL;
        break;
      }

      m_counter++;
    }
    while (true);

    // Standard Newton
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      this->m_Psi_Ref.add(-0.1, this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi_Ref);

      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        this->output_results(path);
      }

      columns& cols = m_table.new_line();
      m_table.insert(cols, MyTable::COUNTER, double(m_counter));
      m_table.insert(cols, MyTable::RES, m_res);
      m_table.insert(cols, MyTable::RESP, m_resp);
      m_table.insert(cols, MyTable::MU, m_mu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if (m_root)
      {
        cout << m_table;
      }
      if (m_res < m_epsilon[1])
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while (true);

    this->do_linear_superposition();
    m_N = Particle_Number(this->m_Psi_Ref);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }

    string filename = path + "log.csv";
    if (m_root)
    {
      m_table.dump_2_file(filename);
    }

    return retval;
  }


  template <int dim>
  void MySolver<dim>::run2b()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    this->setup_system();

    CEigenfunctions<dim> Ef1(m_QN1, m_omega);
    VectorTools::interpolate(this->m_DOF_Handler, Ef1, this->m_Psi[0]);

    this->m_Psi[0] *= 1.0 / sqrt(Particle_Number(this->m_Psi[0]));
    this->m_Psi[1] = 0;

    compute_E_lin(this->m_Psi[0], T, N, W);
    double m_mu_0 = T / N;
    m_mu = ceil(10.0 * m_mu_0) / 10.0 + m_gs / fabs(m_gs) * m_dmu;

    this->output_guess();
    m_results.clear();
    for (int i = 0; i < m_Ndmu; ++i)
    {
      sprintf(shellcmd, "mkdir %.4d/", i);
      if (m_root)
      {
        system(shellcmd);
      }
      sprintf(shellcmd, "%.4d/", i);
      path = shellcmd;

      // nehari
      m_ti = sqrt((m_mu * N - T) / (m_gs * W));

      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "m_mu = " << m_mu << endl;
      pcout << "m_ti = " << m_ti << endl;

      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert(cols, MyTable::MU, m_mu);
      m_results.insert(cols, MyTable::GS, m_gs);
      m_results.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      m_results.insert(cols, MyTable::COUNTER, double(m_counter));
      m_results.insert(cols, MyTable::STATUS, double(status));

      if (status == Status::SUCCESS)
      {
        estimate_error(m_final_error);
        this->output_results(path, "Cfinal");
        string filename = path + "Cfinal.bin";
        this->save(filename);
        filename = path + "Cfinal-1.bin";
        save_one(filename);
        dump_info_xml(path);
        this->m_Psi[0] = this->m_Psi_Ref;
        this->m_Psi[1] = 0;
      }
      else if (status == Status::SLOW_CONV)
      {
        this->m_Psi[1] = 0;
      }
      else
      {
        break;
      }
      compute_E_lin(this->m_Psi_Ref, T, N, W);   // TODO: kommentier mich aus, falls ich kein nehari reset habe
      m_mu += m_gs / fabs(m_gs) * m_dmu;
    }
    if (m_root)
    {
      m_results.dump_2_file("results.csv");
    }
  }

  template<int dim>
  void MySolver<dim>::save_one(string filename)
  {
    double tmp = Particle_Number(this->m_Psi_Ref);
    this->m_Workspace[0] = this->m_Psi_Ref;
    this->m_Workspace[0] *= sqrt(1 / tmp);
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(this->m_DOF_Handler);
    solution_transfer.prepare_for_serialization(this->m_Workspace[0]);
    this->m_Triangulation.save(filename.c_str());
  }

  template <int dim>
  void MySolver<dim>::dump_info_xml(const string path)
  {
    // string filename = path + "info.xml";

    // wofstream fs2;
    // fs2.open(filename);

    // locale utf8_locale("en_US.UTF8");
    // fs2.imbue(utf8_locale);
    // fs2 << L"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<INFO>\n";
    // fs2 << L"<MU>" << m_mu << L"</MU>\n";
    // fs2 << L"<GS>" << m_gs << L"</GS>\n";
    // fs2 << L"<N>" << m_N << "L</N>\n";
    // fs2 << L"<XMIN>" << m_xmin << L"</XMIN>\n";
    // fs2 << L"<XMAX>" << m_xmax << L"</XMAX>\n";
    // fs2 << L"<YMIN>" << m_ymin << L"</YMIN>\n";
    // fs2 << L"<YMAX>" << m_ymax << L"</YMAX>\n";
    // fs2 << L"<ZMIN>" << m_zmin << L"</ZMIN>\n";
    // fs2 << L"<ZMAX>" << m_zmax << L"</ZMAX>\n";
    // fs2 << L"<FINAL_ERROR>" << m_final_error << L"</FINAL_ERROR>\n";
    // fs2 << L"<REVISION>" << STR2(GIT_SHA1) << L"</REVISION>\n";
    // fs2 << L"</INFO>\n";
    // fs2.close();
  }
} // end of namespace

int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run2b();
  }
  return EXIT_SUCCESS;
}