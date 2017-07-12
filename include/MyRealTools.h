
#include <cassert>
#include <complex>
#include <cmath>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>

namespace MyRealTools
{
    using namespace dealii;
}

namespace MyRealTools { namespace MPI
{
    namespace LA
    {
        using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;

    template<int dim>
    double Particle_Number( MPI_Comm mpi_communicator, 
                            const DoFHandler<dim>& dof_handler,
                            const FE_Q<dim>& fe,
                            const LA::MPI::Vector& vec )
    {
        assert( vec.has_ghost_elements() == true );
       
        double tmp = 0;
        
        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<double> vec_vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
           if( cell->is_locally_owned() )
           {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    tmp += fe_values.JxW(qp)*vec_vals[qp]*vec_vals[qp];
                }
            }
        }

        double retval;
        MPI_Allreduce( &tmp, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    return retval;
    }        

    template<int dim>
    void Expectation_value_position( MPI_Comm mpi_communicator, 
                                     const DoFHandler<dim>& dof_handler,
                                     const FE_Q<dim>& fe,
                                     const LA::MPI::Vector& vec,
                                     vector<double>& retval )
    {
        assert( retval.size() == dim );
        assert( vec.has_ghost_elements() == true );
       
        Point<dim> tmp;
        
        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<double> vec_vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
           if( cell->is_locally_owned() )
           {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxWxn = fe_values.JxW(qp)*vec_vals[qp]*vec_vals[qp];
                    Point<dim> spacept = fe_values.quadrature_point(qp);
                    tmp += JxWxn*spacept;
                }
            }
        }

        vector<double> tmpv(dim,0);
        for( unsigned i=0; i<dim; i++ ) tmpv[i] = tmp[i];
        MPI_Allreduce( tmpv.data(), retval.data(), dim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }        

    template<int dim>
    void Expectation_value_width( MPI_Comm mpi_communicator, 
                                  const DoFHandler<dim>& dof_handler,
                                  const FE_Q<dim>& fe,
                                  const LA::MPI::Vector& vec,
                                  const vector<double>& pos,
                                  vector<double>& retval )
    {
        assert( pos.size() == dim );
        assert( retval.size() == dim );
        assert( vec.has_ghost_elements() == true );
       
        Point<dim> tmp;

        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<double> vec_vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxWxn = fe_values.JxW(qp)*vec_vals[qp]*vec_vals[qp];
                    Point<dim> spacept = fe_values.quadrature_point(qp);
                    for( unsigned i=0; i<dim; i++ ) 
                        tmp[i] += JxWxn*(spacept[i]-pos[i])*(spacept[i]-pos[i]);
                }
            }
        }
        vector<double> tmpv(dim,0);
        for( unsigned i=0; i<dim; i++ ) tmpv[i] = tmp[i];
        MPI_Allreduce( tmpv.data(), retval.data(), dim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }

    template <int dim>
    void AssembleSystem_Jacobian( const DoFHandler<dim>& dof_handler,
                                  const FE_Q<dim>& fe,
                                  const ConstraintMatrix& constraints,
                                  const LA::MPI::Vector& vec, 
                                  const Function<dim>& Potential,
                                  const double mu,
                                  const double gs,
                                  LA::MPI::SparseMatrix& matrix )
    {
        assert( vec.has_ghost_elements() == true );

        const QGauss<dim> quadrature_formula(fe.degree+1);
        
        matrix = 0;

        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        vector<double> vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; ++cell )
        {
            if( cell->is_locally_owned() )
            {
                cell_matrix = 0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec,vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    double Q2 = Potential.value(fe_values.quadrature_point(qp)) - mu + 3.0*gs*vals[qp]*vals[qp];

                    for ( unsigned i=0; i<dofs_per_cell; ++i )
                        for ( unsigned j=0; j<dofs_per_cell; ++j )
                            cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, local_dof_indices, matrix);
            }
        }
        matrix.compress(VectorOperation::add);
    }

    template <int dim>
    void AssembleRHS_L2gradient( const DoFHandler<dim>& dof_handler,
                                 const FE_Q<dim>& fe,
                                 const ConstraintMatrix& constraints,
                                 const LA::MPI::Vector& vec, 
                                 const Function<dim>& Potential,
                                 const double mu,
                                 const double gs,
                                 double& res,
                                 LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );

        const QGauss<dim> quadrature_formula(fe.degree+1);
        
        rhs=0;

        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        Vector<double> cell_rhs (dofs_per_cell);

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        vector<Tensor<1, dim> > grad_vals(n_q_points);
        vector<double> vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; ++cell )
        {
            if( cell->is_locally_owned() )
            {
                cell_rhs = 0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);
                fe_values.get_function_gradients(vec, grad_vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                double JxW = fe_values.JxW(qp);
                double Q1 = Potential.value(fe_values.quadrature_point(qp)) - mu + gs*vals[qp]*vals[qp];

                for ( unsigned i=0; i<dofs_per_cell; ++i )
                    cell_rhs(i) += JxW*(grad_vals[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_rhs, local_dof_indices, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        res = rhs.l2_norm();
    }

    template <int dim>
    void AssembleSystem_tangent( const DoFHandler<dim>& dof_handler,
                                 const FE_Q<dim>& fe,
                                 const ConstraintMatrix& constraints,
                                 const LA::MPI::Vector& vec, 
                                 const Function<dim>& Potential,
                                 const double mu,
                                 const double gs,
                                 LA::MPI::SparseMatrix& matrix,
                                 LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        const QGauss<dim> quadrature_formula(fe.degree+1);
        
        matrix=0;
        rhs=0;

        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        Vector<double> cell_rhs (dofs_per_cell);
        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        vector<Tensor<1,dim>> grad_vals(n_q_points);
        vector<double> vals(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; ++cell )
        {
            if( cell->is_locally_owned() )
            {
                cell_rhs = 0;
                cell_matrix = 0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);
                fe_values.get_function_gradients(vec, grad_vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    double Q2 = Potential.value(fe_values.quadrature_point(qp)) - mu + 3.0*gs*vals[qp]*vals[qp];

                    for ( unsigned i=0; i<dofs_per_cell; ++i )
                    {
                        cell_rhs(i) += JxW*vals[qp]*fe_values.shape_value(i,qp);
                        for ( unsigned j=0; j<dofs_per_cell; ++j )
                            cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);   
    }               

}}
