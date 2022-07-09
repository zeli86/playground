
// #define STR1(x) #x
// #define STR2(x) STR1(x)

template <int dim, int N>
cBase<dim, N>::cBase(const std::string& xmlfilename)
  :
  m_ph(xmlfilename),
  m_Triangulation(typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
  m_FE(2),
  m_DOF_Handler(m_Triangulation)
{
  // m_omega = m_ph.Get_Physics("omega");
  // m_gs = m_ph.Get_Physics("gs_1", 0);
  // m_QN1[0] = int(m_ph.Get_Physics("QN1", 0));
  // m_QN1[1] = int(m_ph.Get_Physics("QN1", 1));
  // m_QN1[2] = int(m_ph.Get_Physics("QN1", 2));

  // m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));

  // m_ti = m_ph.Get_Algorithm("ti", 0);
  // m_epsilon = m_ph.Get_Algorithm("epsilon");
  // m_t[0] = m_ti;
  // m_t[1] = m_ti;
  // m_t_guess[0] = m_ti;
  // m_t_guess[1] = m_ti;

  // m_NA = int(m_ph.Get_Algorithm("NA", 0));
  // m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
  // m_dmu = m_ph.Get_Algorithm("dmu", 0);

}

// template <int dim, int N>
// void cBase<dim, N>::save(const std::string& filename)
// {
//   ofstream out(filename);
//   this->m_Psi_Ref.block_write(out);
// }

