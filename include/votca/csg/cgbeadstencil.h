
#ifndef VOTCA_CSG_CGBEADSTENCIL
#define VOTCA_CSG_CGBEADSTENCIL
#include <string>
#include <vector>
#include <votca/tools/types.h>

namespace TOOLS = votca::tools;
namespace votca {

  namespace csg {

  struct CGBeadStencil {
    std::string cg_name_;
    std::string cg_bead_type_;
    TOOLS::byte_t cg_symmetry_;
    std::string mapping_;
    std::vector<std::string> atomic_subbeads_;
    std::vector<double> subbead_weights_;
    std::vector<double> subbead_d_;
  };

  }
}
#endif // VOTCA_CSG_CGBEADSTENCIL
