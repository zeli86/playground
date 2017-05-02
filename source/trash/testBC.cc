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
 Solves curl(curl(u)) + u = f in 2D
 
 exact solution is given by:
 u(0) = cos(pi*x)*sin(pi*y) + C
 u(1) = -sin(pi*x)*cos(pi(y) + C
 
 f(0) = (2*pi^2 + 1)*cos(pi*x)*sin(pi*x) + C
 f(1) = -(2*pi^2 + 1)*sin(pi*x)*cos(pi(y) + C
 
 where C is some constant. To change C, edit bc_constant in ExactSolution and RightHandSide classes.
 
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_refinement.h>

using namespace dealii;
template <int dim>
class MaxwellProblem
{
public:
    MaxwellProblem (const unsigned int order);
    ~MaxwellProblem ();
    void run ();
private:
    double dotprod(const Tensor<1,dim> &A, const Tensor<1,dim> &B) const;
    double dotprod(const Tensor<1,dim> &A, const Vector<double> &B) const;
    void setup_system ();
    void assemble_system ();
    void solve ();
    void process_solution(const unsigned int cycle);
    void output_results_vtk (const unsigned int cycle) const;
    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;
    FE_Nedelec<dim>            fe;
    ConstraintMatrix     constraints;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;
    
    unsigned int p_order;
    unsigned int quad_order;
    
    ConvergenceTable	   convergence_table;
};
// EXACT SOLUTION CLASS
template<int dim>
class ExactSolution : public Function<dim>
{
public:
    ExactSolution() : Function<dim>() {}
    virtual double value (const Point<dim> &p,
                          const unsigned int component) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double> &result) const;
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component) const;
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &values) const;
private:
    const double PI = dealii::numbers::PI;
	const double bc_constant = 0.1;
};
// RIGHT HAND SIDE CLASS
template <int dim>
class RightHandSide :  public Function<dim>
{
public:
    RightHandSide ();
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
private:
    const double PI = dealii::numbers::PI;
	const double bc_constant = 0.1;
};
// DEFINE EXACT SOLUTION MEMBERS
template<int dim>
double ExactSolution<dim>::value(const Point<dim> &p,
                                 const unsigned int component) const
{
    Assert (dim >= 2, ExcNotImplemented());
    AssertIndexRange(component, dim);
    
    double val = -1000;
    switch(component) {
        case 0:	val = cos(PI*p(0))*sin(PI*p(1)) + bc_constant;
        case 1:	val = -sin(PI*p(0))*cos(PI*p(1)) + bc_constant;
    }
    return val;
    
}
template<int dim>
void ExactSolution<dim>::vector_value(const Point<dim> &p,
                                      Vector<double> &result) const
{
    Assert(dim >= 2, ExcNotImplemented());
    result(0) = cos(PI*p(0))*sin(PI*p(1)) + bc_constant;
    result(1) = -sin(PI*p(0))*cos(PI*p(1)) + bc_constant;
    
}
template <int dim>
void ExactSolution<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double> &values,
                                     const unsigned int component) const
{
    Assert (values.size() == points.size(), ExcDimensionMismatch(values.size(), points.size()));
    AssertIndexRange(component, dim);
    for (unsigned int i=0; i<points.size(); ++i)
    {
        const Point<dim> &p = points[i];
        switch(component)
        {
            case 0:
                values[i] = cos(PI*p(0))*sin(PI*p(1)) + bc_constant;
            case 1:
                values[i] = -sin(PI*p(0))*cos(PI*p(1)) + bc_constant;
        }
    }
}
template <int dim>
void ExactSolution<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> >   &values) const
{
    Assert (dim >= 2, ExcNotImplemented());
    Assert (values.size() == points.size(), ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0; i<points.size(); ++i)
    {
        const Point<dim> &p = points[i];
        values[i](0) = cos(PI*p(0))*sin(PI*p(1)) + bc_constant;
        values[i](1) = -sin(PI*p(0))*cos(PI*p(1)) + bc_constant;
    }
}
// END EXACT SOLUTION MEMBERS

// DEFINE RIGHT HAND SIDE MEMBERS
template <int dim>
RightHandSide<dim>::RightHandSide () :
Function<dim> (dim)
{}
template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
    Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));
    Assert (dim >= 2, ExcNotImplemented());
    
    //2D solution
    values(0) = (2*PI*PI + 1)*cos(PI*p(0))*sin(PI*p(1)) + bc_constant;
    values(1) = -(2*PI*PI + 1)*sin(PI*p(0))*cos(PI*p(1)) + bc_constant;
}
template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> >   &value_list) const
{
    Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));
    const unsigned int n_points = points.size();
    for (unsigned int p=0; p<n_points; ++p)
    {
        RightHandSide<dim>::vector_value (points[p], value_list[p]);
    }
}
// END RIGHT HAND SIDE MEMBERS

template <int dim>
MaxwellProblem<dim>::MaxwellProblem (const unsigned int order)
:
dof_handler (triangulation),
fe (order)
{
    p_order = order;
    quad_order = p_order+2;
}
template <int dim>
MaxwellProblem<dim>::~MaxwellProblem ()
{
    dof_handler.clear ();
}

template<int dim>
double MaxwellProblem<dim>::dotprod(const Tensor<1,dim> &A, const Tensor<1,dim> &B) const
{
    double return_val = 0;
    for(unsigned int k = 0; k < dim; k++) {
        return_val += A[k]*B[k];
    }
    return return_val;
}

template<int dim>
double MaxwellProblem<dim>::dotprod(const Tensor<1,dim> &A, const Vector<double> &B) const
{
    double return_val = 0;
    for(unsigned int k = 0; k < dim; k++) {
        return_val += A[k]*B(k);
    }
    return return_val;
}

template <int dim>
void MaxwellProblem<dim>::setup_system ()
{
    dof_handler.distribute_dofs (fe);
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    // FE_Nedelec boundary condition.
    VectorTools::project_boundary_values_curl_conforming(dof_handler, 0, ExactSolution<dim>(), 0, constraints);
    
    constraints.close ();
    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    c_sparsity,
                                    constraints,false);
    
    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void MaxwellProblem<dim>::assemble_system ()
{
    const QGauss<dim>  quadrature_formula(quad_order);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
    
    
    FEValuesExtractors::Vector A(0);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();
    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    RightHandSide<dim>      right_hand_side;
    std::vector<Vector<double> > rhs_values (n_q_points,Vector<double>(dim));
    Tensor<1,dim,double> rhst; 
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_matrix = 0;
        cell_rhs = 0;
        fe_values.reinit (cell);
        right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
	
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  cell_matrix(i,j) += fe_values.JxW(q_point)*(fe_values[A].curl(i,q_point)*fe_values[A].curl(j,q_point)+fe_values[A].value(i,q_point)*fe_values[A].value(j,q_point));
                }
		rhst[0] = rhs_values[q_point][0];
         	rhst[1] = rhs_values[q_point][1];
                cell_rhs(i) += rhst*fe_values[A].value(i,q_point)*fe_values.JxW(q_point);
            }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}

template <int dim>
void MaxwellProblem<dim>::solve ()
{
  /* CG:
    SolverControl           solver_control (1000, 1e-8);
    SolverCG<>              solver (solver_control);
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    solver.solve (system_matrix, solution, system_rhs,
                  preconditioner);
    constraints.distribute (solution);
    */
    
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
  
    A_direct.vmult (solution, system_rhs);
    constraints.distribute (solution);
    
}
template<int dim>
void MaxwellProblem<dim>::process_solution(const unsigned int cycle)
{
    const ExactSolution<dim> exact_solution;
    Vector<double> diff_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler, solution, exact_solution,
                                      diff_per_cell, QGauss<dim>(quad_order), VectorTools::L2_norm);
    const double L2_error = diff_per_cell.l2_norm();
    
    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", triangulation.n_active_cells());
    convergence_table.add_value("dofs", dof_handler.n_dofs());
    convergence_table.add_value("L2 Error", L2_error);
}
template <int dim>
void MaxwellProblem<dim>::output_results_vtk (const unsigned int cycle) const
{
    std::ostringstream filename;
    filename << "solution-" << cycle << ".vtk";
    std::ofstream output (filename.str().c_str());
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    std::vector<std::string> solution_names;
    
    switch (dim)
    {
        case 1:
            solution_names.push_back ("E");
            break;
        case 2:
            solution_names.push_back ("E1");
            solution_names.push_back ("E2");
            break;
        case 3:
            solution_names.push_back ("E1");
            solution_names.push_back ("E2");
            solution_names.push_back ("E3");
            break;
        default:
            Assert (false, ExcNotImplemented());
    }
    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_out.add_data_vector (solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.build_patches ();
    data_out.write_vtk (output);
}
template <int dim>
void MaxwellProblem<dim>::run ()
{
    for (unsigned int cycle=0; cycle<4; ++cycle)
    {
        std::cout << "Cycle " << cycle << ':' << std::endl;
        if (cycle == 0)
        {
            GridGenerator::hyper_cube (triangulation, -1, 1);
            triangulation.refine_global (2);
        }
        else
            triangulation.refine_global (1);
        std::cout << "   Number of active cells:       "
        << triangulation.n_active_cells()
        << std::endl;
        setup_system ();
        std::cout << "   Number of degrees of freedom: "
        << dof_handler.n_dofs()
        << std::endl;
        assemble_system ();
        solve ();
        process_solution (cycle);
        output_results_vtk (cycle);
    }
    convergence_table.set_precision("L2 Error",8);
    convergence_table.set_scientific("L2 Error",true);
    convergence_table.write_text(std::cout);
}
int main ()
{
    try
    {
        deallog.depth_console (0);
        MaxwellProblem<2> maxwell(4);
        maxwell.run ();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
        std::cerr << "Exception on processing: " << std::endl
        << exc.what() << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
        std::cerr << "Unknown exception!" << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
        return 1;
    }
    return 0;
}