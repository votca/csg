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

#ifndef VOTCA_CSG_CGMOLECULEDEF_H
#define VOTCA_CSG_CGMOLECULEDEF_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "exclusionlist.h"
#include "map.h"
#include "molecule.h"
#include <votca/tools/property.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BoundaryCondition;

struct CGBeadInfo {
  std::string cg_name_;
  std::string type_;
  byte_t symmetry_;
  std::string mapping_;
};

struct CGInteractionInfo {
  std::string type_;
  std::string group_;
  std::vector<std::string> bead_names_;
};

/**
 * @brief Definition of coarse grained molecule
 *
 * This class is to define a coarse grained molecule, which includes the
 * topology, mapping, ...
 */
class CGStencil {
 public:
  CGStencil(string cg_molecule_type, string atomistic_molecule_type)
      : cg_molecule_type_(cg_molecule_type),
        atomistic_molecule_type_(atomistic_molecule_type){};

  ~CGStencil(){};

  AddBeadInfo(vector<CGBeadInfo> bead_info) { bead_info_ = bead_info; }

  AddInteractionInfo(vector<CGInteractionInfo> interaction_info) {
    interaction_info_ = interaction_info;
  }

  const vector<CGBeadInfo> &getBeadInfo() { return bead_info_; }

  const vector<CGInteractionInfo> &getInteractionInfo() {
    return interaction_info_;
  }

  const std::string &getCGType() const { return cg_molecule_type_; }

  const std::string &getAtomisticType() const {
    return atomistic_molecule_type_;
  }

 private:
  std::vector<CGBeadInfo> bead_info_;
  std::vector<CGInteractionInfo> interaction_info_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGMOLECULEDEF_H
