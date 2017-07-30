  
#pragma once

#include "exprtk.hpp"

using namespace std;

  template <int dim, int N>
  class CPotential : public dealii::Function<dim> 
  {
    public:
      // Constructor
      CPotential () 
        : 
        Function<dim>(),
        m_t(N),
        m_pos_val(3)
      { 
      }

      // Copy Constructor
      CPotential ( const CPotential& rhs ) 
        : 
        Function<dim>(),
        m_t(N),
        m_pos_val(3)
      { 
        init( rhs.m_all_lambdas, rhs.m_all_potential, rhs.m_constants, rhs.m_T );
      }
      
      // Destructor
      virtual ~CPotential()
      {
        for( auto& i : m_lambdas )
        {
          delete i;
        }
        for( auto& i : m_pot )
        {
          delete i;
        }
      }

      void init(  const vector<string>& all_lambdas, const vector<string>& all_potential, const map<string,double>& constants, const double T )
      {
        assert( all_lambdas.size() != 0 );
        assert( all_potential.size() != 0 );
        assert( all_potential.size() == all_lambdas.size()+1 );

        m_all_lambdas = all_lambdas; 
        m_all_potential = all_potential; 
        m_constants = constants;

        m_no_lam = all_lambdas.size();
        m_dt = T/double(N-1);
        m_T = T;
        m_fak = 1/(m_dt*m_dt);

        for( int i=0; i<N; i++ )
        {
          m_t[i] = double(i)*m_dt;
        }

        try
        {
          // Setup all initial lambda_i(t) with the guess from all_lambdas
          vector<double> tmpvec(N);
          for( auto str : all_lambdas )
          {
            FunctionParser<1> lam;
            lam.initialize( "t", str, constants );
            
            for( int i=0; i<N; i++ )
            {
              Point<1> t(double(i)*m_dt);
              tmpvec[i] = lam.value(t);
            }
            m_lambdas.push_back( new Functions::CSpline<1>(m_t, tmpvec) );
          }

          m_lam_val.resize(all_lambdas.size());
          // Setup the potential and all derivatives with respect to lambda_i(t)
          for( auto str : all_potential )
          {
            m_pot.push_back( new mu::Parser() );
            mu::Parser * p = m_pot.back();
            
            p->SetExpr(str);

            p->DefineVar("x", &(m_pos_val.data()[0]));
            p->DefineVar("y", &(m_pos_val.data()[1]));
            p->DefineVar("z", &(m_pos_val.data()[2]));

            for( int i=0; i<all_lambdas.size(); i++ )
            {
              string tmp = "lam_" + to_string(i);
              p->DefineVar(tmp, &(m_lam_val.data()[i]));
            }
            for( auto i : constants )
            {
              p->DefineConst(i.first, i.second);
            }
          }
        }
        catch (mu::Parser::exception_type &e)
        {
          std::cout << e.GetMsg() << std::endl;
        }
      }

      void reinit( const vector<vector<double>>& new_lambdas )
      {
        assert( m_no_lam == new_lambdas.size() );

        for( auto& i : m_lambdas )
        {
          delete i;
        }          

        for( int i=0; i<m_no_lam; i++ )
        {
          m_lambdas[i] = new Functions::CSpline<1>(m_t, new_lambdas[i]);
        }
      }

      virtual double value ( const Point<dim> &p, const unsigned component = 0) const 
      {
        Point<1> pt(this->get_time());
        double retval=0;

        for( int i=0; i<dim; i++ )
        {
          m_pos_val[i] = p[i]; // setting the spatial coordinate
        }
        int s=0;
        for( auto i : m_lambdas )
        {
          assert( this->get_time() >= 0 && this->get_time() <= m_T );
          m_lam_val[s] = i->value(pt);  // setting lam_i(t)
          s++;
        }        

        try
        {
          retval = m_pot[component]->Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
          std::cout << e.GetMsg() << std::endl;
        }        
      return retval;
      }    

      void output( const std::string& filename )
      {
        ofstream out( filename );

        for( int i=0; i<N; i++ )
        {
          Point<1> pt(double(i)*m_dt);
          out << pt[0] << "\t";
          for( int j=0; j<m_no_lam; j++ )
          {
            out << m_lambdas[j]->value(pt) << ( j+1 == m_no_lam ? "\n" : "\t" );
          }
        }
      }

      bool save( const string& filename )
      {
        std::ofstream out( filename, std::ifstream::binary );

        if( !out.is_open() ) return false;

        out.write(reinterpret_cast<char*>(&m_no_lam), sizeof(int));

        std::vector<double> tmp(N);

        for( int j=0; j<m_no_lam; j++ )
        {
          for( int i=0; i<N; i++ )
          {
            tmp[i] = m_lambdas[j]->value(Point<1>(double(i)*m_dt));
          }
          out.write(reinterpret_cast<char*>(tmp.data()),sizeof(double)*N);
        }
        return true;
      }

      bool load( const string& filename )
      {
        std::ifstream in(filename);

        if( !in.is_open() ) return false;

        int no_lam;
        vector<vector<double>> tmp( m_no_lam, vector<double>(N) );
        in.read( reinterpret_cast<char*>(&no_lam), sizeof(int));

        assert( no_lam == m_no_lam );

        for( int j=0; j<m_no_lam; j++ )
        {
          in.read( reinterpret_cast<char*>(tmp[j].data()), sizeof(double)*N);
        }

        reinit( tmp );

        return true;
      }

      double get_no_lambdas() { return m_no_lam; };

      vector<dealii::Functions::CSpline<1>*> m_lambdas;
    protected:
      int m_no_lam;
      double m_dt;
      double m_T;
      double m_fak;
      vector<mu::Parser*> m_pot;
      vector<mu::Parser*> m_pot;
      vector<double> m_t;
      mutable vector<double> m_pos_val;
      mutable vector<double> m_lam_val;
      vector<string> m_all_lambdas; 
      vector<string> m_all_potential; 
      map<string,double> m_constants;
  };