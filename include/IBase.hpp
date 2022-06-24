
template<typename tDoFHandler, typename tFe, typename tConstraints>
class IBase
{
public:
  virtual tDoFHandler& get_dof_handler() = 0;
  virtual const tFe& get_fe() = 0;
  virtual const tConstraints& get_affine_constraints() = 0;
};