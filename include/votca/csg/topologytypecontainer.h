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

#ifndef VOTCA_CSG_TOPOLOGYTYPECONTAINER_H
#define VOTCA_CSG_TOPOLOGYTYPECONTAINER_H

#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>

namespace votca {
namespace csg {

/** \brief Keeps track of the different topology types
 *
 * The purpose of this class is to track information related the molecule types
 * and the residue types and their associated numbers or ids. All type tracking
 * is done and handled by this class.
 *
 * This class does not track all the molecules residues or beads of a
 * particular type, that is handled by the topology class. This class simply
 * keeps track of the different types that are being used. E.g. If 10 propane
 * molecules are being used the class does not record each molecule, it simply
 * registers that someonewhere in the system a molecule of type "propane"
 * exists.
 *
 **/
class TopologyTypeContainer {
 public:
  TopologyTypeContainer(){};

  /**
   * @brief Determines if a molecule of type **molecule_type** is present in
   * this type container
   *
   * @param[in] molecule_type the name of the molecular type that is to be
   * checked
   *
   * @return bool value indicating if a molecule of the specified type has been
   * registered or not, true it has been registered false it has not.
   */
  bool MoleculeTypeExist(const std::string& molecule_type) const;

  /**
   * @brief Register a particular molecule type
   *
   * @param[in] molecule_type to be registered e.g. "propane", "P3HT", "rubrene"
   */
  void AddMoleculeType(std::string molecule_type);

  /**
   * @brief Determine how many unique types of molecules exist
   *
   * E.g. if there are 10 "hexane", 2 "PTB7" and 1000 water molecules in the
   * system this class would only be responsible for knowing that there are
   * three different molecule types present. As such a count of 3 would be
   * returned.
   *
   * @return - the total number of unique molecular types that have been
   * registered
   */
  size_t MoleculeTypeCount() const { return molecule_types_.size(); }

  const std::unordered_map<std::string, int>& getMoleculeTypes() const {
    return molecule_types_;
  }

  std::vector<std::string> getResidueTypes() const;

  bool ResidueTypeExist(const std::string& residue_type) const;

  void AddResidueType(std::string residue_type);

  /**
   * @brief Returns the number of unique residues that have been registered
   *
   * @return size-t
   */
  size_t ResidueTypeCount() const { return residue_types_.size(); }

  /**
   * @brief Determine if a bead of the type `bead_type` has been registered
   *
   * @param[in] bead_type - these types typically have either elemental symbols
   * such as "H", "O", "Si" or are a type of coarse grained bead such as "CH2",
   * "CH3","OH"
   *
   * @return bool true if the bead type has been registered false otherwise.
   */
  bool BeadTypeExist(const std::string& bead_type) const;

  std::vector<std::string> getBeadTypes() const;

  void AddBeadType(std::string bead_type);

  int getBeadTypeId(std::string bead_type) const {
    assert(bead_types_.count(bead_type) && "Bead type is not recognized");
    return bead_types_.at(bead_type);
  }

  /**
   * @brief Removes memory of all the currently stored bead types
   */
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

#endif  // VOTCA_CSG_TOPOLOGYTYPECONTAINER_H
