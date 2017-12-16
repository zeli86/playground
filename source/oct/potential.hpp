  
#pragma once

using namespace std;

  template <int dim>
  class CPotential : public dealii::Function<dim> 
  {
    public:
      // Constructor
      CPotential () 
        : 
        Function<dim>(),
        m_pos_val(3)
      { 
      }

      // Copy Constructor
      CPotential ( const CPotential& rhs ) 
        : 
        Function<dim>(),
        m_pos_val(3)
      { 
        init( rhs.m_all_lambdas, rhs.m_all_potential, rhs.m_constants, rhs.m_T, rhs.m_N );
      }
      
      CPotential& operator=( const CPotential& rhs )
      {
        init( rhs.m_all_lambdas, rhs.m_all_potential, rhs.m_constants, rhs.m_T, rhs.m_N );
        return *this;
      }

      void init(  const vector<string>& all_lambdas, const vector<string>& all_potential, const map<string,double>& constants, const double T, const int N )
      {
        assert( all_lambdas.size() != 0 );
        assert( all_potential.size() != 0 );
        assert( all_potential.size() == all_lambdas.size()+1 );

        m_all_lambdas = all_lambdas; 
        m_all_potential = all_potential; 
        m_constants = constants;

        m_N = N;
        m_no_lam = all_lambdas.size();
        m_dt = T/double(m_N-1);
        m_T = T;
        m_fak = 1/(m_dt*m_dt);

        m_lambdas.reinit(m_N,m_no_lam);

        try
        {
          // Setup all initial lambda_i(t) with the guess from all_lambdas
          for( int s=0; s<m_no_lam; s++ )
          {
            FunctionParser<1> lam;
            lam.initialize( "t", m_all_lambdas[s], constants );
            
            for( int i=0; i<N; i++ )
            {
              m_lambdas(i,s) = lam.value(Point<1>(double(i)*m_dt));
            }
          }

          m_lam_val.resize(all_lambdas.size());
          // Setup the potential and all derivatives with respect to lambda_i(t)
          for( auto str : all_potential )
          {
            mu::Parser newp;  
            
            newp.SetExpr(str);

            newp.DefineVar("x", &(m_pos_val.data()[0]));
            newp.DefineVar("y", &(m_pos_val.data()[1]));
            newp.DefineVar("z", &(m_pos_val.data()[2]));

            for( int i=0; i<all_lambdas.size(); i++ )
            {
              string tmp = "lam_" + to_string(i);
              newp.DefineVar(tmp, &(m_lam_val.data()[i]));
            }
            for( auto i : constants )
            {
              newp.DefineConst(i.first, i.second);
            }
            m_pot.push_back(newp);
          }
        }
        catch (mu::Parser::exception_type &e)
        {
          std::cout << e.GetMsg() << std::endl;
        }
      }

      virtual double value ( const Point<dim> &p, const unsigned component = 0) const 
      {
        double retval=0;
        int ti = this->get_time() / m_dt;

        for( int i=0; i<dim; i++ )
        {
          m_pos_val[i] = p[i]; // setting the spatial coordinate
        }
        for( int s=0 ; s<m_no_lam; s++ )
        {
          m_lam_val[s] = m_lambdas(ti,s);  // setting lam_i(t)
        }        

        try
        {
          retval = m_pot[component].Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
          std::cout << e.GetMsg() << std::endl;
        }        
      return retval;
      }

      bool save( const string& filename )
      {
        std::ofstream out( filename, std::ifstream::binary );

        if( !out.is_open() ) return false;

        out.write(reinterpret_cast<char*>(&m_no_lam), sizeof(int));
        out.write(reinterpret_cast<char*>(&m_N), sizeof(int));
        
        std::vector<double> tmp(m_N*m_no_lam);

        for( int s=0; s<m_no_lam; s++ )
        {
          for( int i=0; i<m_N; i++ )
          {
            tmp[i+m_N*s] = m_lambdas(i,s);
          }
        }
        out.write(reinterpret_cast<char*>(tmp.data()),sizeof(double)*m_N*m_no_lam);
        return true;
      }

      bool load( const string& filename )
      {
        std::ifstream in(filename);

        if( !in.is_open() ) return false;

        int no_lam, Nt;
        in.read( reinterpret_cast<char*>(&no_lam), sizeof(int));
        in.read( reinterpret_cast<char*>(&Nt), sizeof(int));

        assert( no_lam == m_no_lam );
        
        std::vector<double> tmp(m_N*m_no_lam);
       
        in.read( reinterpret_cast<char*>(tmp.data()), sizeof(double)*m_no_lam*m_N);

        for( int i=0; i<m_N*m_no_lam; i++ )
        {
          int s = i / m_N;
          int ti = i - s*m_N;
          m_lambdas(ti,s) = tmp[i];
        }
        return true;
      }          

      void output( const std::string& filename )
      {
        ofstream out( filename );

        for( int i=0; i<m_N; i++ )
        {
          Point<1> pt(double(i)*m_dt);
          out << pt[0] << "\t";
          for( int j=0; j<m_no_lam; j++ )
          {
            out << m_lambdas(i,j) << ( j+1 == m_no_lam ? "\n" : "\t" );
          }
        }
      }

      void add( const double tau, const LAPACKFullMatrix<double>& direction )
      {
        vector<vector<double>> new_lambdas(m_no_lam,vector<double>(m_N,1));

        for( int s=0; s<m_no_lam; s++ )
          for( int ti=0; ti<m_N; ti++ )
          {
            m_lambdas(ti,s) += tau*direction(ti,s);
          }
      }

      double get_no_lambdas() { return m_no_lam; };

    protected:
      int m_no_lam; // number of lambdas
      int m_N; // total number of time steps
      double m_dt;
      double m_T;
      double m_fak;
      vector<mu::Parser> m_pot; 
      mutable vector<double> m_pos_val;
      mutable vector<double> m_lam_val;
      vector<string> m_all_lambdas; 
      vector<string> m_all_potential; 
      map<string,double> m_constants;
      LAPACKFullMatrix<double> m_lambdas;
    };