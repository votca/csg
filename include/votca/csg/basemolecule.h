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

#ifndef _VOTCA_CSG_BASEMOLECULE_H
#define _VOTCA_CSG_BASEMOLECULE_H

#include <assert.h>
#include <cassert>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "basebead.h"
#include "beadstructure.h"
#include "topologyitem.h"

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

// class Interaction;

/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
template <class T>
class BaseMolecule : public BeadStructure<T> {
 public:
  /// get the molecule ID
  int getId() const { return id_.getId(); }

  /// get the name of the molecule
  const std::string &getName() const { return name_.getName(); }

  /// set the name of the molecule
  void setName(const std::string &name) { name_.setName(name); }

  // The bead already has a name handled by beadstructure, but we need to
  // override it
  void AddBead(T *bead);

  /**
   * @brief Returns the ids of the beads with the name `name`
   *
   * @param[in] name
   *
   * @return unordered set of all the ids
   */
  std::unordered_set<int> getBeadIdsByName(const std::string &name) const;

  /**
   * @brief Returns the ids of the beads that have the type given by name
   *
   * @param[in] name
   *
   * @return unordered set with the ids of the beads
   */
  std::unordered_set<int> getBeadIdsByType(const std::string &name) const;

  /**
   * @brief Returns the bead type of the bead given by id
   *
   * @param[in] id
   *
   * @return bead pointer to the bead
   */
  const std::string &getBeadType(const int &id) const;

  /**
   * @brief Grabs the position of the bead with id `id`
   *
   * @param[in] id
   *
   * @return vector describting the position of the bead
   */
  const TOOLS::vec &getBeadPosition(const int &id) const;

  /**
   * @brief Determines the beads name provided the id
   *
   * @param[in] id of the bead
   *
   * @return string name of the bead
   */
  const std::string getBeadName(int id) const;

 protected:
  TOOLS::Identity<int> id_;
  TOOLS::Name name_;

  std::unordered_map<std::string, std::unordered_set<int>> bead_name_and_ids_;
};

template <class T>
void BaseMolecule<T>::AddBead(T *bead) {
  assert(!BeadStructure<T>::beads_.count(bead->getId()) &&
         "Cannot add a bead to the basemolecule"
         " when it has been previously added.");

  BeadStructure<T>::AddBead(bead);
  bead_name_and_ids_[bead->getName()].insert(bead->getId());
  bead->setMoleculeId(getId());
}

template <class T>
const std::string BaseMolecule<T>::getBeadName(int id) const {
  assert(BeadStructure<T>::beads_.count(id) &&
         "Cannot get bead name for bead id because "
         "is is not stored in the base molecule.");
  return BeadStructure<T>::beads_.at(id)->getName();
}

template <class T>
const std::string &BaseMolecule<T>::getBeadType(const int &id) const {
  assert(BeadStructure<T>::beads_.count(id) &&
         "Cannot get bead type with id beacuse "
         "bead is not stored in base molecule.");
  return BeadStructure<T>::beads_.at(id)->getType();
}

template <class T>
const TOOLS::vec &BaseMolecule<T>::getBeadPosition(const int &id) const {
  assert(BeadStructure<T>::beads_.count(id) &&
         "Cannot get bead position with id because "
         "bead is not stored in the base molecule.");
  return BeadStructure<T>::beads_.at(id)->getPos();
}

template <class T>
std::unordered_set<int> BaseMolecule<T>::getBeadIdsByName(
    const std::string &name) const {
  assert(bead_name_and_ids_.count(name) &&
         "BaseMolecule does not contain any "
         "beads with name ");
  return bead_name_and_ids_.at(name);
}

template <class T>
std::unordered_set<int> BaseMolecule<T>::getBeadIdsByType(
    const std::string &type) const {

  std::unordered_set<int> bead_ids;
  for (const std::pair<const int, T *> &id_and_bead :
       BeadStructure<T>::beads_) {
    if (type.compare(id_and_bead.second->getType()) == 0) {
      bead_ids.insert(id_and_bead.first);
    }
  }
  return bead_ids;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BASEMOLECULE_H