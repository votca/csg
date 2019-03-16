/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
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

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_CGBEADSTENCIL
