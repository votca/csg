
#ifndef VOTCA_CSG_CGBEADINFO
#define VOTCA_CSG_CGBEADINFO
#include <string>
#include <vector>
#include <votca/tools/types.h>

namespace TOOLS = votca::tools;
namespace votca {

  namespace csg {

  struct CGBeadInfo {
    std::string cg_name_;
    std::string cg_bead_type_;
    TOOLS::byte_t cg_symmetry_;
    std::string mapping_;
    std::vector<std::string> atomic_subbeads_;
    std::vector<double> subbead_weights_;
    std::vector<double> subbead_d_;
  };

  struct CGInteractionInfo {
    std::string type_;
    std::string group_;
    std::vector<std::string> bead_names_;
  };

  }
}
#endif // VOTCA_CSG_CGBEADINFO
