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

#include "../../include/votca/csg/topologytypecontainer.h"
namespace votca {
namespace csg {

bool TopologyTypeContainer::MoleculeTypeExist(std::string molecule_type) const {
  return molecule_types_.count(molecule_type);
}
void TopologyTypeContainer::AddMoleculeType(std::string molecule_type) {
  if (molecule_types_.count(molecule_type) == 0) {
    molecule_types_[molecule_type] = molecule_types_.size();
  }
}

bool TopologyTypeContainer::ResidueTypeExist(std::string residue_type) const {
  return residue_types_.count(residue_type);
}
void TopologyTypeContainer::AddResidueType(std::string residue_type) {
  if (residue_types_.count(residue_type) == 0) {
    residue_types_[residue_type] = residue_types_.size();
  }
}
bool TopologyTypeContainer::BeadTypeExist(std::string bead_type) const {
  return bead_types_.count(bead_type);
}
void TopologyTypeContainer::AddBeadType(std::string bead_type) {
  if (bead_types_.count(bead_type) == 0) {
    bead_types_[bead_type] = bead_types_.size();
  }
}

}  // namespace csg
}  // namespace votca
