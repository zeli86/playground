
#include <cassert>
#include <complex>
#include <cmath>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

namespace MyRealTools
{
  using namespace std;
  using namespace dealii;

  template<int dim>
  double Particle_Number(const DoFHandler<dim>& dof_handler,
                         const FE_Q<dim>& fe,
                         const Vector<double>& vec)
  {
    double retval = 0;
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();

    std::vector<double> vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(vec, vals);

      for (unsigned qp = 0; qp < n_q_points; qp++)
      {
        retval += fe_values.JxW(qp) * vals[qp] * vals[qp];
      }
    }
    return retval;
  }

  template <int dim>
  void compute_stepsize(const DoFHandler<dim>& dof_handler,
                        const FE_Q<dim>& fe,
                        const Function<dim>& Potential,
                        const Vector<double>& psi,
                        const Vector<double>& direction,
                        const double mu,
                        const double gs,
                        double& retval)
  {
    retval = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<double> u(n_q_points);
    vector<double> d(n_q_points);
    vector<Tensor<1, dim>> u_grad(n_q_points);
    vector<Tensor<1, dim>> d_grad(n_q_points);

    double total_int[4] = {};

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(psi, u);
        fe_values.get_function_values(direction, d);
        fe_values.get_function_gradients(psi, u_grad);
        fe_values.get_function_gradients(direction, d_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          double JxW = fe_values.JxW(qp);
          double Q = Potential.value(fe_values.quadrature_point(qp)) - mu;

          total_int[0] += JxW * (u_grad[qp] * d_grad[qp] + Q * u[qp] * d[qp] + gs * u[qp] * u[qp] * u[qp] * d[qp]); // tau^0
          total_int[1] += JxW * (d_grad[qp] * d_grad[qp] + Q * d[qp] * d[qp] + 3 * gs * u[qp] * u[qp] * d[qp] * d[qp]); // tau^1
          total_int[2] += JxW * d[qp] * d[qp] * d[qp] * u[qp]; // tau^2
          //local_int[3] += JxW*d[qp]*d[qp]*d[qp]*d[qp]; // tau^3
        }
      }
    }

    //total_int[1] = local_int[0];
    total_int[2] *= (3 * gs);
    //total_int[3] *= (gs);

    //printf( "%e, %e, %e, %e\n", total_int[0], total_int[1], total_int[2], total_int[2] );

    double xm = -0.5 * (total_int[1] - sqrt(total_int[1] * total_int[1] - 4 * total_int[2] * total_int[0])) / total_int[2];
    double xp = -0.5 * (total_int[1] + sqrt(total_int[1] * total_int[1] - 4 * total_int[2] * total_int[0])) / total_int[2];

    //retval=std::min(fabs(std::min( fabs(xp),fabs(xm))), 1.0);

    if (fabs(xp) < fabs(xm))
    {
      retval = xp;
    }
    else
    {
      retval = xm;
    }
    if (fabs(retval) > 1)
    {
      retval = retval / fabs(retval);
    }
    if (isnan(retval))
    {
      retval = 1;
    }
  }
}

namespace MyRealTools
{
  namespace MPI
  {
    namespace LA
    {
      using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;

    template<int dim>
    double Particle_Number(MPI_Comm mpi_communicator,
                           const DoFHandler<dim>& dof_handler,
                           const FE_Q<dim>& fe,
                           const LA::MPI::Vector& vec)
    {
      assert(vec.has_ghost_elements() == true);

      double tmp = 0;

      const QGauss<dim>  quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      vector<double> vec_vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            tmp += fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
          }
        }
      }

      double retval;
      MPI_Allreduce(&tmp, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      return retval;
    }

    template<int dim>
    void Expectation_value_position(MPI_Comm mpi_communicator,
                                    const DoFHandler<dim>& dof_handler,
                                    const FE_Q<dim>& fe,
                                    const LA::MPI::Vector& vec,
                                    vector<double>& retval)
    {
      assert(retval.size() == dim);
      assert(vec.has_ghost_elements() == true);

      Point<dim> tmp;

      const QGauss<dim>  quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      vector<double> vec_vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxWxn = fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
            Point<dim> spacept = fe_values.quadrature_point(qp);
            tmp += JxWxn * spacept;
          }
        }
      }

      vector<double> tmpv(dim, 0);
      for (unsigned i = 0; i < dim; i++)
      {
        tmpv[i] = tmp[i];
      }
      MPI_Allreduce(tmpv.data(), retval.data(), dim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }

    template<int dim>
    void Expectation_value_width(MPI_Comm mpi_communicator,
                                 const DoFHandler<dim>& dof_handler,
                                 const FE_Q<dim>& fe,
                                 const LA::MPI::Vector& vec,
                                 const vector<double>& pos,
                                 vector<double>& retval)
    {
      assert(pos.size() == dim);
      assert(retval.size() == dim);
      assert(vec.has_ghost_elements() == true);

      Point<dim> tmp;

      const QGauss<dim>  quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      vector<double> vec_vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            double JxWxn = fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
            Point<dim> spacept = fe_values.quadrature_point(qp);
            for (unsigned i = 0; i < dim; i++)
            {
              tmp[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
            }
          }
        }
      }
      vector<double> tmpv(dim, 0);
      for (unsigned i = 0; i < dim; i++)
      {
        tmpv[i] = tmp[i];
      }
      MPI_Allreduce(tmpv.data(), retval.data(), dim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }

    template <int dim>
    void AssembleSystem_Jacobian(const DoFHandler<dim>& dof_handler,
                                 const FE_Q<dim>& fe,
                                 const  AffineConstraints<double>& constraints,
                                 const LA::MPI::Vector& vec,
                                 const Function<dim>& Potential,
                                 const double mu,
                                 const double gs,
                                 LA::MPI::SparseMatrix& matrix)
    {
      assert(vec.has_ghost_elements() == true);

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      matrix = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);
            double Q2 = Potential.value(fe_values.quadrature_point(qp)) - mu + 3.0 * gs * vals[qp] * vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
              }
          }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix, local_dof_indices, matrix);
        }
      }
      matrix.compress(VectorOperation::add);
    }

    template <int dim>
    void AssembleRHS_L2gradient(const DoFHandler<dim>& dof_handler,
                                const FE_Q<dim>& fe,
                                const AffineConstraints<double>& constraints,
                                const LA::MPI::Vector& vec,
                                const Function<dim>& Potential,
                                const double mu,
                                const double gs,
                                double& res,
                                LA::MPI::Vector& rhs)
    {
      assert(vec.has_ghost_elements() == true);
      assert(rhs.has_ghost_elements() == false);

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      rhs = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      Vector<double> cell_rhs(dofs_per_cell);

      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      vector<Tensor<1, dim>> grad_vals(n_q_points);
      vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vals);
          fe_values.get_function_gradients(vec, grad_vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - mu + gs * vals[qp] * vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += JxW * (grad_vals[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_rhs, local_dof_indices, rhs);
        }
      }
      rhs.compress(VectorOperation::add);
      res = rhs.l2_norm();
    }

    template <int dim>
    void AssembleSystem_tangent(const DoFHandler<dim>& dof_handler,
                                const FE_Q<dim>& fe,
                                const  AffineConstraints<double>& constraints,
                                const LA::MPI::Vector& vec,
                                const Function<dim>& Potential,
                                const double mu,
                                const double gs,
                                LA::MPI::SparseMatrix& matrix,
                                LA::MPI::Vector& rhs)
    {
      assert(vec.has_ghost_elements() == true);
      assert(rhs.has_ghost_elements() == false);

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      matrix = 0;
      rhs = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      Vector<double> cell_rhs(dofs_per_cell);
      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      vector<Tensor<1, dim>> grad_vals(n_q_points);
      vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vals);
          fe_values.get_function_gradients(vec, grad_vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);
            double Q2 = Potential.value(fe_values.quadrature_point(qp)) - mu + 3.0 * gs * vals[qp] * vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += JxW * vals[qp] * fe_values.shape_value(i, qp);
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
              }
            }
          }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
        }
      }
      rhs.compress(VectorOperation::add);
      matrix.compress(VectorOperation::add);
    }

    template<int dim>
    void Compute_mu(MPI_Comm mpi_communicator,
                    const DoFHandler<dim>& dof_handler,
                    const FE_Q<dim>& fe,
                    const Function<dim>& Potential,
                    const LA::MPI::Vector& vec,
                    double& retval)
    {
      assert(vec.has_ghost_elements() == true);

      retval = 0;

      double tmp[] = {0, 0};

      const QGauss<dim>  quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned n_q_points = quadrature_formula.size();
      vector<double> vec_vals(n_q_points);
      vector<Tensor<1, dim>> vec_grads(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          fe_values.get_function_gradients(vec, vec_grads);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);
            double valq = vec_vals[qp] * vec_vals[qp];
            tmp[0] += JxW * (vec_grads[qp] * vec_grads[qp] + Potential.value(fe_values.quadrature_point(qp)) * valq);
            tmp[1] += JxW * valq;
          }
        }
      }

      double fak = 0;
      MPI_Allreduce(&(tmp[1]), &fak, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);

      tmp[1] = tmp[1] / fak;

      MPI_Allreduce(&(tmp[0]), &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }

    template <int dim>
    void compute_stepsize(MPI_Comm mpi_communicator,
                          const DoFHandler<dim>& dof_handler,
                          const FE_Q<dim>& fe,
                          const Function<dim>& Potential,
                          const LA::MPI::Vector& psi,
                          const LA::MPI::Vector& direction,
                          const double mu,
                          const double gs,
                          double& retval)
    {
      assert(psi.has_ghost_elements() == true);
      assert(direction.has_ghost_elements() == true);

      retval = 0;

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      vector<double> u(n_q_points);
      vector<double> d(n_q_points);
      vector<Tensor<1, dim>> u_grad(n_q_points);
      vector<Tensor<1, dim>> d_grad(n_q_points);
      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      double total_int[4] = {};
      double local_int[4] = {};

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(psi, u);
          fe_values.get_function_values(direction, d);
          fe_values.get_function_gradients(psi, u_grad);
          fe_values.get_function_gradients(direction, d_grad);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);
            double Q = Potential.value(fe_values.quadrature_point(qp)) - mu;

            local_int[0] += JxW * (u_grad[qp] * d_grad[qp] + Q * u[qp] * d[qp] + gs * u[qp] * u[qp] * u[qp] * d[qp]); // tau^0
            local_int[1] += JxW * (d_grad[qp] * d_grad[qp] + Q * d[qp] * d[qp] + 3 * gs * u[qp] * u[qp] * d[qp] * d[qp]); // tau^1
            local_int[2] += JxW * d[qp] * d[qp] * d[qp] * u[qp]; // tau^2
            //local_int[3] += JxW*d[qp]*d[qp]*d[qp]*d[qp]; // tau^3
          }
        }
      }

      //local_int[1] = local_int[0];
      local_int[2] *= (3 * gs);
      //local_int[3] *= (gs);

      MPI_Allreduce(local_int, total_int, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      //printf( "%e, %e, %e, %e\n", total_int[0], total_int[1], total_int[2], total_int[2] );

      double xm = -0.5 * (total_int[1] - sqrt(total_int[1] * total_int[1] - 4 * total_int[2] * total_int[0])) / total_int[2];
      double xp = -0.5 * (total_int[1] + sqrt(total_int[1] * total_int[1] - 4 * total_int[2] * total_int[0])) / total_int[2];

      //retval=std::min(fabs(std::min( fabs(xp),fabs(xm))), 1.0);

      if (fabs(xp) < fabs(xm))
      {
        retval = xp;
      }
      else
      {
        retval = xm;
      }
      if (fabs(retval) > 1)
      {
        retval = retval / fabs(retval);
      }
      if (isnan(retval))
      {
        retval = 1;
      }

    }

    template <int dim>
    void orthonormalize(MPI_Comm mpi_communicator,
                        const DoFHandler<dim>& dof_handler,
                        const FE_Q<dim>& fe,
                        const  AffineConstraints<double>& constraints,
                        const LA::MPI::Vector& vec1,
                        const LA::MPI::Vector& vec2,
                        LA::MPI::Vector& retval)
    {
      assert(vec1.has_ghost_elements() == true);
      assert(vec2.has_ghost_elements() == true);
      assert(retval.has_ghost_elements() == false);

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      vector<double> u1(n_q_points);
      vector<double> u2(n_q_points);
      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      double local_int[2] = {};
      double total_int[2] = {};

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec1, u1);
          fe_values.get_function_values(vec2, u2);

          for (unsigned int qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);

            local_int[0] += JxW * u1[qp] * u2[qp];
            local_int[1] += JxW * u2[qp] * u2[qp];
          }
        }
      }
      MPI_Allreduce(local_int, total_int, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);

      double fak = total_int[0] / total_int[1];
      retval = 0;
      retval.add(1, vec1, -fak, vec2);
      constraints.distribute(retval);
    }
  }
}
