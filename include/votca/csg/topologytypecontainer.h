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

#ifndef _VOTCA_CSG_TOPOLOGYTYPECONTAINER_H
#define _VOTCA_CSG_TOPOLOGYTYPECONTAINER_H

#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>

namespace votca {
namespace csg {

/**
 * \brief Keeps track of the different topology types
 *
 * The purpose of this class is to track information related the molecule types
 * and the residue types and their associated numbers or ids
 **/
class TopologyTypeContainer {
 public:
  TopologyTypeContainer(){};

  bool MoleculeTypeExist(const std::string molecule_type) const;

  void AddMoleculeType(std::string molecule_type);

  size_t MoleculeTypeCount() const { return molecule_types_.size(); }

  const std::unordered_map<std::string, int>& getMoleculeTypes() const {
    return molecule_types_;
  }

  std::vector<std::string> getResidueTypes() const;

  bool ResidueTypeExist(const std::string residue_type) const;

  void AddResidueType(std::string residue_type);

  size_t ResidueTypeCount() const { return residue_types_.size(); }

  bool BeadTypeExist(const std::string bead_type) const;

  std::vector<std::string> getBeadTypes() const;

  void AddBeadType(std::string bead_type);

  int getBeadTypeId(std::string bead_type) const {
    assert(bead_types_.count(bead_type) && "Bead type is not recognized");
    return bead_types_.at(bead_type);
  }
  void Clear() {
    bead_types_.clear();
    residue_types_.clear();
    molecule_types_.clear();
  }

 private:
  std::unordered_map<std::string, int> bead_types_;
  std::unordered_map<std::string, int> residue_types_;
  std::unordered_map<std::string, int> molecule_types_;
};
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_TOPOLOGYTYPECONTAINER_H
