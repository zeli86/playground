
template <int dim, int N, class T, bool bComplex>
class cBaseMPI
{
public:
  virtual ~cBaseMPI()
  {
    m_DOF_Handler.clear();
  };

protected:

  typedef typename dealii::LinearAlgebraPETSc::MPI::Vector aVector;
  typedef typename dealii::LinearAlgebraPETSc::MPI::SparseMatrix aSparseMatrix;
  typedef typename dealii::parallel::distributed::Triangulation<dim> aTriangulation;

  friend double myfunc<dim>(const std::vector<double>&, std::vector<double>&, void*);

  aVector m_Psi_Ref; // non ghosted
  aVector m_System_RHS; // non ghosted
  aVector m_Search_Direction; // non ghosted
  aVector m_Workspace_NG;
  dealii::Vector<double> m_error_per_cell;

  std::array<aVector, N> m_Psi; // non ghosted
  std::array<aVector, N> m_Workspace; // ghosted

  aSparseMatrix m_System_Matrix;
  std::vector<int> m_QN1;
};

template <int dim, int N, class T, bool bComplex>
void cBaseMPI<dim, N, T, bComplex>::update_workspace()
{
  for (int i = 0; i < N; i++)
  {
    m_constraints.distribute(m_Psi[i]);
    m_Workspace[i] = m_Psi[i];
  }
}

template <int dim, int N, class T, bool bComplex>
void cBaseMPI<dim, N, T, bComplex>::save(const std::string& filename)
{
  // m_constraints.distribute(m_Psi_Ref);
  // parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(m_DOF_Handler);
  // solution_transfer.prepare_for_serialization(m_Psi_Ref);
  // m_Triangulation.save(filename.c_str());
}

template <int dim, int N, class T, bool bComplex>
void cBaseMPI<dim, N, T, bComplex>::output_results(const std::string& path, std::string prefix)
{
  // dealii::TimerOutput::Scope timing_section(m_computing_timer, "");

  // std::string filename;

  // dealii::Vector<float> subdomain(m_Triangulation.n_active_cells());
  // for (unsigned int i = 0; i < subdomain.size(); ++i)
  // {
  //   subdomain(i) = m_Triangulation.locally_owned_subdomain();
  // }

  // dealii::DataOut<dim> data_out;
  // data_out.attach_dof_handler(m_DOF_Handler);
  // data_out.add_data_vector(m_Psi_Ref, "Psi_sol");
  // data_out.add_data_vector(m_error_per_cell, "error per cell");
  // data_out.add_data_vector(subdomain, "subdomain");
  // data_out.build_patches();

  // filename = path + prefix + "-" + dealii::Utilities::int_to_string(m_counter, 5) + ".vtu";
  // data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
}

template <int dim, int N, class T, bool bComplex>
void cBaseMPI<dim, N, T, bComplex>::output_guess()
{
  // dealii::TimerOutput::Scope timing_section(m_computing_timer, "");

  // this->update_workspace();

  // // CPotential Potential_fct ( m_omega, m_QN1[2] );
  // // VectorTools::interpolate (this->m_DOF_Handler, Potential_fct, this->m_Workspace_NG );
  // // this->m_constraints.distribute(this->m_Workspace_NG);
  // // this->m_Psi_Ref=this->m_Workspace_NG;

  // dealii::DataOut<dim> data_out;
  // data_out.attach_dof_handler(this->m_DOF_Handler);
  // data_out.add_data_vector(this->m_Workspace[0], "Psi_0");  // todo : loop
  // data_out.add_data_vector(this->m_Workspace[1], "Psi_1");
  // // data_out.add_data_vector (this->m_Psi_Ref, "m_Potential");
  // data_out.build_patches();
  // data_out.write_vtu_in_parallel("guess.vtu", mpi_communicator);
}

