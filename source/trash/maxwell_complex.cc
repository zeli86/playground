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

/*
Solving time-harmonic Maxwell equations with complex coefficients.

\nabla \times \nabla \times E - i E = J in Omega,
n \nabla E = n \nabla g at boundary.

J = (1,1,1) + i (1,1,1),
g = (1,1,1) + i (1,1,1)

*/
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>
#include <complex>

using namespace dealii;

template <int dim>
class ComputeEr : public DataPostprocessorVector<dim>
{
  public:
    ComputeEr ();
    virtual void compute_derived_quantities_vector (
        const std::vector< Vector< double > > &uh,
        const std::vector< std::vector< Tensor< 1, dim > > > &duh,
        const std::vector< std::vector< Tensor< 2, dim > > > &dduh,
        const std::vector< Point< dim > > &normals,
        const std::vector<Point<dim> > &evaluation_points,
        std::vector< Vector< double > > &computed_quantities) const;
};

template <int dim>
ComputeEr<dim>::ComputeEr ()
  : DataPostprocessorVector<dim> ("Er", update_values)
{}

template <int dim>
void ComputeEr<dim>::compute_derived_quantities_vector (
    const std::vector< Vector< double > >                  &uh,
    const std::vector< std::vector< Tensor< 1, dim > > >  & /*duh*/,
    const std::vector< std::vector< Tensor< 2, dim > > >  & /*dduh*/,
    const std::vector< Point< dim > >                     & /*normals*/,
    const std::vector<Point<dim> >                        & /*evaluation_points*/,
    std::vector< Vector< double > >                        &computed_quantities
    ) const
{
  Assert(computed_quantities.size() == uh.size(),
      ExcDimensionMismatch (computed_quantities.size(), uh.size()));

  for (unsigned int i = 0; i < computed_quantities.size(); i++) {
    Assert(computed_quantities[i].size() == dim,
        ExcDimensionMismatch (computed_quantities[i].size(), dim));
    Assert(uh[i].size() == dim * 2, ExcDimensionMismatch (uh[i].size(), dim * 2));
    for (unsigned int j = 0; j < dim; ++j) {
      computed_quantities[i][j] = uh[i][j];
    }
  }
}

template <int dim>
class ComputeEi : public DataPostprocessorVector<dim>
{
  public:
    ComputeEi ();
    virtual void compute_derived_quantities_vector (
        const std::vector< Vector< double > > &uh,
        const std::vector< std::vector< Tensor< 1, dim > > > &duh,
        const std::vector< std::vector< Tensor< 2, dim > > > &dduh,
        const std::vector< Point< dim > > &normals,
        const std::vector<Point<dim> > &evaluation_points,
        std::vector< Vector< double > > &computed_quantities) const;
};

template <int dim>
ComputeEi<dim>::ComputeEi ()
  : DataPostprocessorVector<dim> ("Ei", update_values)
{}

template <int dim>
void ComputeEi<dim>::compute_derived_quantities_vector (
    const std::vector< Vector< double > >                  &uh,
    const std::vector< std::vector< Tensor< 1, dim > > >  & /*duh*/,
    const std::vector< std::vector< Tensor< 2, dim > > >  & /*dduh*/,
    const std::vector< Point< dim > >                     & /*normals*/,
    const std::vector<Point<dim> >                        & /*evaluation_points*/,
    std::vector< Vector< double > >                        &computed_quantities
    ) const
{
  Assert(computed_quantities.size() == uh.size(),
      ExcDimensionMismatch (computed_quantities.size(), uh.size()));

  for (unsigned int i = 0; i < computed_quantities.size(); i++) {
    Assert(computed_quantities[i].size() == dim,
        ExcDimensionMismatch (computed_quantities[i].size(), dim));
    Assert(uh[i].size() == dim * 2, ExcDimensionMismatch (uh[i].size(), dim * 2));
    for (unsigned int j = 0; j < dim; ++j) {
      computed_quantities[i][j] = uh[i][j + dim];
    }
  }
}

template <int dim>
class MaxwellProblem
{
  public:
    MaxwellProblem ();
    ~MaxwellProblem ();

    void run ();

  private:
    void create_coarse_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void postprocess();

    Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FESystem<dim> fe;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};

template <int dim>
MaxwellProblem<dim>::MaxwellProblem ()
  : dof_handler (triangulation), fe(FE_Nedelec<dim>(0), 1, FE_Nedelec<dim>(0), 1)
{
}

template <int dim>
MaxwellProblem<dim>::~MaxwellProblem ()
{
  dof_handler.clear ();
}

template <int dim>
void MaxwellProblem<dim>::create_coarse_grid ()
{
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.refine_global(dim == 2 ? 7 : 4);
}

template <int dim>
void MaxwellProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  VectorTools::project_boundary_values_curl_conforming(dof_handler, 0,
      ConstantFunction<dim>(1.0, fe.n_components()), 0, constraints);
  VectorTools::project_boundary_values_curl_conforming(dof_handler, dim,
      ConstantFunction<dim>(1.0, fe.n_components()), 0, constraints);
  constraints.close();

  CompressedSparsityPattern csp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
  sparsity_pattern.copy_from (csp);

  system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void MaxwellProblem<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
      update_values | update_gradients |
      update_quadrature_points | update_JxW_values);

  const size_t dofs_per_cell = fe.dofs_per_cell;
  const size_t n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const FEValuesExtractors::Vector Er(0);
  const FEValuesExtractors::Vector Ei(dim);

  size_t i, j, q_point;

  Point<dim> j_re;
  Point<dim> j_im;

  for (i = 0; i < dim; ++i) {
    j_re[i] = j_im[i] = 1.0;
  }

  auto cell = dof_handler.begin_active();
  for (; cell != dof_handler.end(); ++cell) {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;

    for (i = 0; i < dofs_per_cell; ++i) {
      auto base_index_i = fe.system_to_base_index(i);
      for (j = 0; j < dofs_per_cell; ++j) {
        auto base_index_j = fe.system_to_base_index(j);
        if (base_index_i.first.first == base_index_j.first.first) {
          for (q_point = 0; q_point < n_q_points; ++q_point) {
            cell_matrix(i, j) += (
                fe_values[Er].value(i, q_point) * fe_values[Er].value(j, q_point)
              + fe_values[Ei].value(i, q_point) * fe_values[Ei].value(j, q_point)
                ) * fe_values.JxW(q_point);
          }
        } else {
          for (q_point = 0; q_point < n_q_points; ++q_point) {
            cell_matrix(i, j) += (
              - fe_values[Er].value(i, q_point) * fe_values[Ei].value(j, q_point)
              + fe_values[Ei].value(i, q_point) * fe_values[Er].value(j, q_point)
                ) * fe_values.JxW(q_point);
          }
        }
      }

      for (q_point = 0; q_point < n_q_points; ++q_point) {
        cell_rhs(i) += fe_values[Er].value(i, q_point) *
          j_re * fe_values.JxW(q_point);
        cell_rhs(i) += fe_values[Ei].value(i, q_point) *
          j_im * fe_values.JxW(q_point);
      }
    }

    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global (cell_matrix, cell_rhs,
        local_dof_indices, system_matrix, system_rhs, false);
  }
}

template <int dim>
void MaxwellProblem<dim>::solve ()
{
  SparseDirectUMFPACK A_direct;
  A_direct.initialize(system_matrix);

  A_direct.vmult (solution, system_rhs);
  constraints.distribute (solution);
}

template <int dim>
void MaxwellProblem<dim>::postprocess ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler, QGauss<dim - 1>(2),
      typename FunctionMap<dim>::type(), solution, estimated_error_per_cell);

  ComputeEr<dim> Er;
  ComputeEi<dim> Ei;
  DataOut<dim,DoFHandler<dim> > data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(estimated_error_per_cell, "error");
  data_out.add_data_vector(solution, Er);
  data_out.add_data_vector(solution, Ei);
  data_out.build_patches();

  const std::string filename = "solution_complex.vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}

template <int dim>
void MaxwellProblem<dim>::run ()
{
  create_coarse_grid ();

  setup_system ();
  assemble_system ();
  solve ();
  postprocess();
}

int main (int argc, char** argv)
{
  MaxwellProblem<2> Maxwell_problem;
  Maxwell_problem.run ();

  return 0;
}
