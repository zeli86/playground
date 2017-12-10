/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */

 template <int dim, int no_time_steps>
 void MySolver<dim,no_time_steps>::compute_beta()
 {
   for( int s=0; s<m_potential.get_no_lambdas(); s++ )
   {
     m_beta[s] = 0;
     for( int ti=1; ti<no_time_steps-1; ti++ )
     {
       //m_beta[s] += m_grad(ti,s)*m_grad(ti,s)/(m_old_direction(ti,s)*(m_grad(ti,s)-m_old_grad(ti,s)));
       m_beta[s] += m_grad(ti,s)*m_grad(ti,s)/(m_old_grad(ti,s)*m_old_grad(ti,s));
       //if( m_root ) printf( "%d\t%d\t%g\t%g\t%g\n", s, ti, m_grad(ti,s), m_old_direction(ti,s), m_old_grad(ti,s) );
     }      
   }    
 }

 template <int dim, int no_time_steps>
 double MySolver<dim,no_time_steps>::compute_dot_product( const LAPACKFullMatrix<double>& mat1, const LAPACKFullMatrix<double>& mat2 )
 {
   double retval=0;
   for( int s=0; s<m_potential.get_no_lambdas(); s++ )
   {
     for( int ti=1; ti<no_time_steps-1; ti++ )
     {
       retval += mat1(ti,s)*mat2(ti,s);
     }
   }
   return retval*m_dt;
 }

 template <int dim, int no_time_steps>
 void MySolver<dim,no_time_steps>::compute_correction( const int ex )
 {
   const int nolam = m_potential.get_no_lambdas();
   
   double retval;

   const QGauss<dim>  quadrature_formula(fe.degree+1);
   FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

   const unsigned n_q_points = quadrature_formula.size();
   vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
   vector<Vector<double>> p(n_q_points,Vector<double>(2));

   // loop over all lambdas
   for( int s=0; s<nolam; s++ )
   {
     m_grad(0,s) = 0;
     m_grad(no_time_steps-1,s) = 0;
     // loop over all time steps
     for( int ti=1; ti<no_time_steps-1; ti++ )
     {
       double tmp1=0;
 
       m_workspace = m_all_Psi[ti];
       m_workspace_2 = m_all_p[ti];
 
       m_potential.set_time( double(ti)*m_dt );
       typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
       for (; cell!=endc; cell++ )
       {
         if( cell->is_locally_owned() )
         {
           fe_values.reinit (cell);
           fe_values.get_function_values( m_workspace, Psi );
           fe_values.get_function_values( m_workspace_2, p );
           
           for ( unsigned qp=0; qp<n_q_points; qp++ )
           {
             tmp1 -= fe_values.JxW(qp) * m_potential.value(fe_values.quadrature_point(qp), s+1) * (Psi[qp][0]*p[qp][0]+Psi[qp][1]*p[qp][1]);
           }
         }
       }
       //printf( "(%d) %d %d %g\n",m_rank,  s, ti, tmp1 );
       //tmp1 -= m_potential.laplacian(s);
       
       m_grad(ti,s) = tmp1;
       //if( m_root ) printf( "%d %d %g\n", s, ti, m_grad(ti,s) );
     }
   }
   m_grad *= (m_dt*m_dt);
   
   // minus Laplace
   LAPACKFullMatrix<double> lap(no_time_steps,no_time_steps);
   for( int i=0; i<no_time_steps; i++ ) lap(i,i) = 2;
   for( int i=1; i<no_time_steps-1; i++ ) lap(i,i+1) = -1; // rechte Nebendiagonale 
   for( int i=1; i<no_time_steps-1; i++ ) lap(i,i-1) = -1; // linke Nebendiagonale
   
   lap.compute_lu_factorization();
   lap.apply_lu_factorization(m_grad,false);

   //vector<vector<double>> new_lambdas(nolam,vector<double>(no_time_steps));

   for( int s=0; s<nolam; s++ )
   {
     m_norm_grad[s] = 0;
     for( int ti=1; ti<no_time_steps; ti++ )
     {
       m_norm_grad[s] += m_grad(ti,s)*m_grad(ti,s);
     }
     m_norm_grad[s] = sqrt(m_norm_grad[s]*m_dt);   
   }    

   /*
   ofstream out( "direction_" + to_string(ex) + ".txt" );
   for( int ti=1; ti<no_time_steps; ti++ )
   {
     out << double(ti)*m_dt << "\t";
     for( int s=0; s<nolam; s++ )
     {
        out << m_direction(ti,s) << ( s+1 == nolam ? "\n" : "\t" );
     }
   }

   ofstream out2( "old_direction_" + to_string(ex) + ".txt" );
   for( int ti=1; ti<no_time_steps; ti++ )
   {
     out2 << double(ti)*m_dt << "\t";
     for( int s=0; s<nolam; s++ )
     {
        out2 << m_old_direction(ti,s) << ( s+1 == nolam ? "\n" : "\t" );
     }
   }
*/
   ofstream out3( "grad_" + to_string(ex) + ".txt" );
   for( int ti=1; ti<no_time_steps; ti++ )
   {
     out3 << double(ti)*m_dt << "\t";
     for( int s=0; s<nolam; s++ )
     {
        out3 << m_grad(ti,s) << ( s+1 == nolam ? "\n" : "\t" );
     }
   }    
/*
   ofstream out4( "oldgrad_" + to_string(ex) + ".txt" );
   for( int ti=1; ti<no_time_steps; ti++ )
   {
     out4 << double(ti)*m_dt << "\t";
     for( int s=0; s<nolam; s++ )
     {
        out4 << m_old_grad(ti,s) << ( s+1 == nolam ? "\n" : "\t" );
     }
   }    
*/
 }

