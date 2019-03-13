
#ifndef VOTCA_CSG_CGINTERACTIONSTENCIL
#define VOTCA_CSG_CGINTERACTIONSTENCIL
#include <string>
#include <vector>

namespace votca {
  namespace csg {

  struct CGInteractionStencil {
    std::string type_;
    std::string group_;
    std::vector<std::string> bead_names_;
  };

  }
}
#endif // VOTCA_CSG_CGINTERACTIONSTENCIL
