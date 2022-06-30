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

#include "cBaseMPI.h"

  template <int dim>
  class MySolver : public cBaseMPI<dim, 2>
  {
  public:
    explicit MySolver(const std::string);

    void run();
    void run2b();

    double m_I[8];
  protected:
    int DoIter(string = "");

    void make_grid_custom();

    void assemble_rhs();
    void assemble_system();

    double Particle_Number(LA::MPI::Vector&);

    void solve();
    void compute_contributions();
    void compute_E_lin(LA::MPI::Vector&, double&, double&, double&);
    void estimate_error(double&);

    MyTable m_table;
    MyTable m_results;

    using cBaseMPI<dim, 2>::mpi_communicator;
    using cBaseMPI<dim, 2>::m_root;
    using cBaseMPI<dim, 2>::m_rank;
    using cBaseMPI<dim, 2>::m_t;
    using cBaseMPI<dim, 2>::m_ti;
    using cBaseMPI<dim, 2>::m_N;
    using cBaseMPI<dim, 2>::m_omega;
    using cBaseMPI<dim, 2>::m_rMu;
    using cBaseMPI<dim, 2>::m_gs;
    using cBaseMPI<dim, 2>::m_counter;
    using cBaseMPI<dim, 2>::pcout;
    using cBaseMPI<dim, 2>::m_computing_timer;
    using cBaseMPI<dim, 2>::m_ph;
    using cBaseMPI<dim, 2>::m_final_error;
    using cBaseMPI<dim, 2>::m_NA;
    using cBaseMPI<dim, 2>::m_Ndmu;
    using cBaseMPI<dim, 2>::m_dmu;
    using cBaseMPI<dim, 2>::m_QN1;
    using cBaseMPI<dim, 2>::m_res;
    using cBaseMPI<dim, 2>::m_resp;
    using cBaseMPI<dim, 2>::m_res_old;
    using cBaseMPI<dim, 2>::m_epsilon;
    using cBaseMPI<dim, 2>::m_maxiter;
    using cBaseMPI<dim, 2>::m_total_no_cells;
    using cBaseMPI<dim, 2>::m_total_no_active_cells;
    using cBaseMPI<dim, 2>::m_global_refinement;
  };

  /*************************************************************************************************/
  /**
   * Constructor
   */

  template <int dim>
  MySolver<dim>::MySolver(const std::string xmlfilename) : cBaseMPI<dim, 2>(xmlfilename)
  {
  }

  template <int dim>
  void MySolver<dim>::compute_E_lin(LA::MPI::Vector& vec, double& T, double& N, double& W)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->m_constraints.distribute(vec);
    this->m_Workspace[0] = vec;

    CPotential Potential(m_omega, m_QN1[2]);
    const QGauss<dim>  quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<Tensor<1, dim>> vec_grad(n_q_points);

    double JxW, vec_val_q;
    double T1 = 0.0, N1 = 0.0, W1 = 0.0;

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vec_vals);
        fe_values.get_function_gradients(this->m_Workspace[0], vec_grad);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          vec_val_q = vec_vals[q_point] * vec_vals[q_point];
          JxW = fabs(fe_values.quadrature_point(q_point)[1]) * fe_values.JxW(q_point);
          T1 += JxW * (vec_grad[q_point] * vec_grad[q_point] + Potential.value(fe_values.quadrature_point(q_point)) * vec_val_q);
          N1 += JxW * vec_val_q;
          W1 += JxW * vec_val_q * vec_val_q;
        }
      }
    }

    T1 *= _FAKTOR_;
    N1 *= _FAKTOR_;
    W1 *= _FAKTOR_;

    MPI_Allreduce(&T1, &T, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce(&N1, &N, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    MPI_Allreduce(&W1, &W, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);

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
    vector<double> vec_vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], vec_vals);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          tmp1 += fabs(fe_values.quadrature_point(q_point)[1]) * fe_values.JxW(q_point) * vec_vals[q_point] * vec_vals[q_point];
        }
      }
    }
    tmp1 *= _FAKTOR_;

    MPI_Allreduce(&tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);

    return retval;
  }

  template <int dim>
  void MySolver<dim>::estimate_error(double& err)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential Potential(m_omega, m_QN1[2]);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);
    //const QGauss<dim-1> face_quadrature_formula(this->m_FE.degree+1);

    this->m_System_RHS = 0;
    this->m_System_Matrix = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);
    //FEFaceValues<dim> fe_face_values ( fe, face_quadrature_formula, update_gradients|update_values|update_quadrature_points|update_normal_vectors|update_JxW_values);

    const unsigned int dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    //const unsigned int n_face_q_points = face_quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_rhs = 0;
        cell_matrix = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Psi_Ref, Psi_ref);
        fe_values.get_function_gradients(this->m_Psi_Ref, Psi_ref_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp) * fabs(fe_values.quadrature_point(qp)[1]);
          const double pq = m_gs * Psi_ref[qp] * Psi_ref[qp];
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + pq;

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW * (Psi_ref_grad[qp] * fe_values.shape_grad(i, qp) + Q1 * Psi_ref[qp] * fe_values.shape_value(i, qp));
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp);
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
    VectorTools::integrate_difference(this->m_DOF_Handler, this->m_Workspace[0], ZeroFunction<dim>(), this->m_error_per_cell, QGauss<dim>(this->m_FE.degree + 2), VectorTools::L2_norm);
    const double total_local_error = this->m_error_per_cell.l2_norm();
    err = std::sqrt(Utilities::MPI::sum(total_local_error * total_local_error, MPI_COMM_WORLD));

  }


  template <int dim>
  void MySolver<dim>::assemble_rhs()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential Potential(m_omega, m_QN1[2]);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_System_RHS = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Psi_Ref, Psi_ref);
        fe_values.get_function_gradients(this->m_Psi_Ref, Psi_ref_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp) * fabs(fe_values.quadrature_point(qp)[1]);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + m_gs * Psi_ref[qp] * Psi_ref[qp];

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW * (Psi_ref_grad[qp] * fe_values.shape_grad(i, qp) + Q1 * Psi_ref[qp] * fe_values.shape_value(i, qp));
          }
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

    CPotential Potential(m_omega, m_QN1[2]);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_System_Matrix = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        cell_matrix = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Psi_Ref, Psi_ref);
        fe_values.get_function_gradients(this->m_Psi_Ref, Psi_ref_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp) * fabs(fe_values.quadrature_point(qp)[1]);
          const double Q2 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + 3.0 * m_gs * Psi_ref[qp] * Psi_ref[qp];

          for (unsigned i = 0; i < dofs_per_cell; ++i)
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        this->m_constraints.distribute_local_to_global(cell_matrix, local_dof_indices, this->m_System_Matrix);
      }
    }
    this->m_System_Matrix.compress(VectorOperation::add);

  }

  template <int dim>
  void MySolver<dim>::solve()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    pcout << "Solving..." << endl;

    SolverControl solver_control;

    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(this->m_System_Matrix, this->m_Search_Direction, this->m_System_RHS);
    this->m_constraints.distribute(this->m_Search_Direction);

  }

  template <int dim>
  void MySolver<dim>::make_grid_custom()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    Point<dim, double> pt1(0, 0);
    Point<dim, double> pt2(m_ph.Get_Mesh("xrange", 1), m_ph.Get_Mesh("yrange", 1));

    CPotential Potential_fct(m_omega, m_QN1[2]);

    GridGenerator::hyper_rectangle(this->m_Triangulation, pt2, pt1);
    this->m_Triangulation.refine_global(6);

    //double isovalues[] = {70,60,50,40};
    //double isovalues[] = {60,50,40};
    double isovalues[] = {60, 50};

    for (unsigned int step = 0; step < sizeof(isovalues) / sizeof(double); step++)
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->m_Triangulation.begin_active(), endc = this->m_Triangulation.end();
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

    typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->m_Triangulation.begin_active(), endc = this->m_Triangulation.end();
    for (; cell != endc; ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        const Point<dim> face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary() && !(face_center[1] == 0))
        {
          cell->face(f)->set_all_boundary_ids(1);
        }
      }
    }

    unsigned int tmp1[2], tmp2[2];
    tmp1[0] = this->m_Triangulation.n_cells();
    tmp1[1] = this->m_Triangulation.n_active_cells();

    MPI_Allreduce(tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);
    /*
        GridOutFlags::Msh opt(true, true);
        string filename = "grid-" + Utilities::int_to_string(this->m_Triangulation.locally_owned_subdomain(), 4);
        ofstream output ((filename + ".msh").c_str());
        GridOut grid_out;
        grid_out.set_flags(opt);
        grid_out.write_msh (this->m_Triangulation,output);
    */
    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];

  }


  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->update_workspace();

    CPotential Potential_fct(m_omega, m_QN1[2]);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 2);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim>> Psi_1_grad(n_q_points);
    vector<Tensor<1, dim>> Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    double locint[8] = {};

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], Psi_1);
        fe_values.get_function_values(this->m_Workspace[1], Psi_2);
        fe_values.get_function_gradients(this->m_Workspace[0], Psi_1_grad);
        fe_values.get_function_gradients(this->m_Workspace[1], Psi_2_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp) * fabs(fe_values.quadrature_point(qp)[1]);
          const double p12 = Psi_1[qp] * Psi_2[qp];
          const double p1q = Psi_1[qp] * Psi_1[qp];
          const double p2q = Psi_2[qp] * Psi_2[qp];
          const double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_rMu;

          locint[0] += JxW * p1q * p1q;
          locint[1] += JxW * p2q * p2q;
          locint[2] += JxW * p1q * p2q;
          locint[3] += JxW * p1q * p12;
          locint[4] += JxW * p2q * p12;
          locint[5] += JxW * (Psi_1_grad[qp] * Psi_1_grad[qp] + Q * p1q);
          locint[6] += JxW * (Psi_2_grad[qp] * Psi_2_grad[qp] + Q * p2q);
          locint[7] += JxW * (Psi_1_grad[qp] * Psi_2_grad[qp] + Q * p12);
        }
      }
    }

    for (int i = 0; i < 5; ++i)
    {
      locint[i] *= m_gs;
    }
    for (int i = 0; i < 8; ++i)
    {
      locint[i] *= _FAKTOR_;
    }

    MPI_Allreduce(locint, m_I, 8, MPI_DOUBLE, MPI_SUM, mpi_communicator);

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
      pcout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      this->m_Psi[1].add(-1e-4 * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
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
      m_table.insert(cols, MyTable::MU, m_rMu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      m_table.insert(cols, MyTable::total_no_cells, double(m_total_no_cells));
      m_table.insert(cols, MyTable::total_no_active_cells, double(m_total_no_active_cells));

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
      m_table.insert(cols, MyTable::MU, m_rMu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);

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
  void MySolver<dim>::run()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    make_grid_custom();
    this->setup_system();

    CEigenfunctions Ef1(m_QN1, m_omega);
    //CEigenfunctions Ef2( m_QN2, f );

    VectorTools::interpolate(this->m_DOF_Handler, Ef1, this->m_Psi[0]);
    this->m_Psi[0] *= 1.0 / sqrt(Particle_Number(this->m_Psi[0]));
    this->m_Psi[1] = 0;

    compute_E_lin(this->m_Psi[0], T, N, W);
    this->output_guess();

    m_rMu = T / N + m_gs / fabs(m_gs) * m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_rMu = " << m_rMu << endl;

    status = DoIter();

    pcout << "L2_norm of this->m_Psi_Ref: " << m_N << endl;

    if (status == Status::SUCCESS)
    {
      this->output_results("", "final");
      this->dump_info_xml();
    }

    if (m_root)
    {
      ofstream ofs("log.txt");
      //ofs << m_table;
    }
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

    CEigenfunctions Ef1(m_QN1, m_omega);
    VectorTools::interpolate(this->m_DOF_Handler, Ef1, this->m_Psi[0]);

    this->m_Psi[0] *= 1.0 / sqrt(Particle_Number(this->m_Psi[0]));
    //this->m_Psi[1] = m_Psi_1;
    this->m_Psi[1] = 0;

    compute_E_lin(this->m_Psi[0], T, N, W);
    double m_rMu_0 = T / N;
    m_rMu = ceil(10.0 * m_rMu_0) / 10.0 + m_gs / fabs(m_gs) * m_dmu;

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
      //m_ti = sqrt((m_rMu*N-T)/(4.0*m_gs*W)); // if this->m_Psi[1] == m_Psi_1
      m_ti = sqrt((m_rMu * N - T) / (m_gs * W));

      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "m_rMu = " << m_rMu << endl;
      pcout << "m_ti = " << m_ti << endl;

      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert(cols, MyTable::MU, m_rMu);
      m_results.insert(cols, MyTable::GS, m_gs);
      m_results.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      m_results.insert(cols, MyTable::COUNTER, double(m_counter));
      m_results.insert(cols, MyTable::STATUS, double(status));

      if (status == Status::SUCCESS)
      {
        estimate_error(m_final_error);
        this->output_results(path, "final");
        this->save(path + "final.bin");
        this->dump_info_xml(path);
        //Interpolate_R_to_C( path + "Cfinal.bin" );
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
      m_rMu += m_gs / fabs(m_gs) * m_dmu;
    }
    if (m_root)
    {
      m_results.dump_2_file("results.csv");
    }
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
