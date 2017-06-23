
#include <cassert>
#include <complex>
#include <cmath>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>

namespace MPI { namespace MyComplexTools
{
    namespace LA
    {
        using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;

    template <int dim>
    void AssembleSystem_mulvz( const DoFHandler<dim>& dof_handler,
                               const FESystem<dim>& fe,
                               const ConstraintMatrix& constraints,
                               const LA::MPI::Vector& vec, 
                               const std::complex<double> z, 
                               LA::MPI::SparseMatrix& matrix, 
                               LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        matrix=0;
        rhs=0;

        const QGauss<dim> quadrature_formula(fe.degree+1);
        const FEValuesExtractors::Scalar rt (0);
        const FEValuesExtractors::Scalar it (1);
        
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs (dofs_per_cell);
        vector<Vector<double>> vals(n_q_points,Vector<double>(2));

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        const double a = std::real(z);
        const double b = std::imag(z);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                cell_rhs=0;
                cell_matrix=0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);

                    for ( unsigned i=0; i<dofs_per_cell; i++ )
                    {
                        for ( unsigned j=0; j<dofs_per_cell; j++ )
                        {
                            cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));                        
                        }
                        double c = vals[qp][0];
                        double d = vals[qp][1];
                        cell_rhs(i) += JxW*((a*c-b*d)*fe_values[rt].value(i,qp) + (b*c+a*d)*fe_values[it].value(i,qp));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);
    }

    template <int dim>
    void AssembleSystem_NL_Step( const DoFHandler<dim>& dof_handler,
                                 const FESystem<dim>& fe,
                                 const ConstraintMatrix& constraints,
                                 const LA::MPI::Vector& vec, 
                                 const double gamdt, 
                                 LA::MPI::SparseMatrix& matrix, 
                                 LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        matrix=0;
        rhs=0;

        const QGauss<dim> quadrature_formula(fe.degree+1);
        const FEValuesExtractors::Scalar rt (0);
        const FEValuesExtractors::Scalar it (1);
        
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs (dofs_per_cell);
        vector<Vector<double>> vals(n_q_points,Vector<double>(2));

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                cell_rhs=0;
                cell_matrix=0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    double c = vals[qp][0];
                    double d = vals[qp][1];
                    double phi = -gamdt * (vals[qp][0]*vals[qp][0] + vals[qp][1]*vals[qp][1]);
                    double a, b;
                    sincos( phi, &b, &a );

                    for ( unsigned i=0; i<dofs_per_cell; i++ )
                    {
                        for ( unsigned j=0; j<dofs_per_cell; j++ )
                        {
                            cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));                        
                        }
                        cell_rhs(i) += JxW*((a*c-b*d)*fe_values[rt].value(i,qp) + (b*c+a*d)*fe_values[it].value(i,qp));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);
    }   

    template <int dim>
    void AssembleSystem_NL_Step( const DoFHandler<dim>& dof_handler,
                                 const FESystem<dim>& fe,
                                 const ConstraintMatrix& constraints,
                                 const LA::MPI::Vector& vec, 
                                 const Function<dim>& Potential,                                  
                                 const double dt,
                                 const double gam, 
                                 LA::MPI::SparseMatrix& matrix, 
                                 LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        matrix=0;
        rhs=0;

        const QGauss<dim> quadrature_formula(fe.degree+1);
        const FEValuesExtractors::Scalar rt (0);
        const FEValuesExtractors::Scalar it (1);
        
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs (dofs_per_cell);
        vector<Vector<double>> vals(n_q_points,Vector<double>(2));

        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                cell_rhs=0;
                cell_matrix=0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    double c = vals[qp][0];
                    double d = vals[qp][1];
                    double phi = -dt * ( gam * (vals[qp][0]*vals[qp][0] + vals[qp][1]*vals[qp][1]) - Potential.value(fe_values.quadrature_point(qp)) );
                    double a, b;
                    sincos( phi, &b, &a );

                    for ( unsigned i=0; i<dofs_per_cell; i++ )
                    {
                        for ( unsigned j=0; j<dofs_per_cell; j++ )
                        {
                            cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));                        
                        }
                        cell_rhs(i) += JxW*((a*c-b*d)*fe_values[rt].value(i,qp) + (b*c+a*d)*fe_values[it].value(i,qp));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);
    }       

    template <int dim>
    void AssembleSystem_LIN_Step( const DoFHandler<dim>& dof_handler,
                                  const FESystem<dim>& fe,
                                  const ConstraintMatrix& constraints,
                                  const LA::MPI::Vector& vec, 
                                  const Function<dim>& Potential, 
                                  const double dt, 
                                  LA::MPI::SparseMatrix& matrix, 
                                  LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        matrix=0;
        rhs=0;

        const QGauss<dim> quadrature_formula(fe.degree+1);
        const FEValuesExtractors::Scalar rt (0);
        const FEValuesExtractors::Scalar it (1);
        const double dth = 0.5*dt;
        
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs (dofs_per_cell);
        vector<Vector<double>> vals(n_q_points,Vector<double>(2));
        vector<vector<Tensor<1,dim>>> vals_grad(n_q_points, vector<Tensor<1,dim>>(2));
        
        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                cell_matrix = 0;
                cell_rhs = 0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);
                fe_values.get_function_gradients(vec, vals_grad);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    double pot = Potential.value(fe_values.quadrature_point(qp));

                    for( unsigned i=0; i<dofs_per_cell; i++ )
                    {
                        for( unsigned j=0; j<dofs_per_cell; j++ )
                        {
                        cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) - dth*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp)) + 
                                                    fe_values[it].value(i,qp)*fe_values[it].value(j,qp) + dth*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)) );
                        }
                        cell_rhs(i) += JxW*(vals[qp][0]*fe_values[rt].value(i,qp) + dth*(vals_grad[qp][1]*fe_values[rt].gradient(i,qp) + pot*vals[qp][1]*fe_values[rt].value(i,qp)) + 
                                            vals[qp][1]*fe_values[it].value(i,qp) - dth*(vals_grad[qp][0]*fe_values[it].gradient(i,qp) + pot*vals[qp][0]*fe_values[it].value(i,qp)));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);
    }    

    template <int dim>
    void AssembleSystem_LIN_Step( const DoFHandler<dim>& dof_handler,
                                  const FESystem<dim>& fe,
                                  const ConstraintMatrix& constraints,
                                  const LA::MPI::Vector& vec, 
                                  const double dt, 
                                  LA::MPI::SparseMatrix& matrix, 
                                  LA::MPI::Vector& rhs )
    {
        assert( vec.has_ghost_elements() == true );
        assert( rhs.has_ghost_elements() == false );
        
        matrix=0;
        rhs=0;

        const QGauss<dim> quadrature_formula(fe.degree+1);
        const FEValuesExtractors::Scalar rt (0);
        const FEValuesExtractors::Scalar it (1);
        const double dth = 0.5*dt;
        
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

        const unsigned dofs_per_cell = fe.dofs_per_cell;
        const unsigned n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs (dofs_per_cell);
        vector<Vector<double>> vals(n_q_points,Vector<double>(2));
        vector<vector<Tensor<1,dim>>> vals_grad(n_q_points, vector<Tensor<1,dim>>(2));
        
        vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for ( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                cell_matrix = 0;
                cell_rhs = 0;

                fe_values.reinit (cell);
                fe_values.get_function_values(vec, vals);
                fe_values.get_function_gradients(vec, vals_grad);

                for ( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);

                    for( unsigned i=0; i<dofs_per_cell; i++ )
                    {
                        for( unsigned j=0; j<dofs_per_cell; j++ )
                        {
                        cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) - dth*fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + 
                                                    fe_values[it].value(i,qp)*fe_values[it].value(j,qp) + dth*fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp));
                        }
                        cell_rhs(i) += JxW*(vals[qp][0]*fe_values[rt].value(i,qp) + dth*vals_grad[qp][1]*fe_values[rt].gradient(i,qp) + 
                                            vals[qp][1]*fe_values[it].value(i,qp) - dth*vals_grad[qp][0]*fe_values[it].gradient(i,qp));
                    }
                }
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
            }
        }
        rhs.compress(VectorOperation::add);   
        matrix.compress(VectorOperation::add);
    }        

    template<int dim>
    double Particle_Number( MPI_Comm mpi_communicator, 
                            const DoFHandler<dim>& dof_handler,
                            const FESystem<dim>& fe,
                            const LA::MPI::Vector& vec )
    {
        assert( vec.has_ghost_elements() == true );
       
        double tmp = 0;
        
        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
           if( cell->is_locally_owned() )
           {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    tmp += fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
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
                                     const FESystem<dim>& fe,
                                     const LA::MPI::Vector& vec,
                                     vector<double>& retval )
    {
        assert( retval.size() == dim );
        assert( vec.has_ghost_elements() == true );
       
        Point<dim> tmp;
        
        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
           if( cell->is_locally_owned() )
           {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxWxn = fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
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
                                  const FESystem<dim>& fe,
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
        vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
            if( cell->is_locally_owned() )
            {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxWxn = fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
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

    template<int dim>
    void Expectation_value_momentum( MPI_Comm mpi_communicator, 
                                     const DoFHandler<dim>& dof_handler,
                                     const FESystem<dim>& fe,
                                     const LA::MPI::Vector& vec,
                                     vector<double>& retval )
    {
        assert( retval.size() == dim );
        assert( vec.has_ghost_elements() == true );
       
        Point<dim> tmp;
        
        const QGauss<dim>  quadrature_formula(fe.degree+1);
        FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values);

        const unsigned n_q_points = quadrature_formula.size();
        vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
        vector<vector<Tensor<1,dim>>> vec_grads(n_q_points, vector<Tensor<1,dim>>(2));
    
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for( ; cell!=endc; cell++ )
        {
           if( cell->is_locally_owned() )
           {
                fe_values.reinit (cell);
                fe_values.get_function_values( vec, vec_vals );
                for( unsigned qp=0; qp<n_q_points; qp++ )
                {
                    double JxW = fe_values.JxW(qp);
                    for( unsigned i=0; i<dim; i++ ) 
                        tmp[i] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][i] - vec_vals[qp][1]*vec_grads[qp][0][i]);
                }
            }
        }

        vector<double> tmpv(dim);
        for( unsigned i=0; i<dim; i++ ) tmpv[i] = tmp[i];
        MPI_Allreduce( tmpv.data(), retval.data(), dim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    }

}}