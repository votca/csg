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

#pragma once
#ifndef VOTCA_CSG_MOLECULE_H
#define VOTCA_CSG_MOLECULE_H

#include "basemolecule.h"
#include "bead.h"
#include <assert.h>
#include <map>
#include <string>
#include <vector>

namespace votca {
namespace csg {

class Interaction;

/**
    \brief Information about molecules.

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule : public BaseMolecule<Bead> {
 public:
  Molecule(){};
  /**
   * @brief Grabs all beads that have the label given by `label`
   *
   * @param label string of the form ResidueName:ResidueNumber:AtomName
   *
   * @return an unordered set with the ids of the beads that match the label
   */
  std::vector<int> getBeadIdsByLabel(const std::string &label) const;

  /// Add an interaction to the molecule
  void AddInteraction(Interaction *ic) { _interactions.push_back(ic); }

  const std::vector<Interaction *> Interactions() const {
    return _interactions;
  }

  tools::StructureParameters getParameters() const;

 private:
  std::vector<Interaction *> _interactions;

  /// constructor
  Molecule(int id, std::string molecule_type) {
    id_.setId(id);
    type_.setName(molecule_type);
  }

  friend class CSG_Topology;
};

inline std::vector<int> Molecule::getBeadIdsByLabel(
    const std::string &label) const {
  std::vector<int> bead_ids;
  for (const std::pair<const int, Bead *> &id_and_bead : beads_) {
    if (label.compare(id_and_bead.second->getLabel()) == 0) {
      bead_ids.push_back(id_and_bead.first);
    }
  }
  return bead_ids;
}

inline tools::StructureParameters Molecule::getParameters() const {
  tools::StructureParameters params;
  params.set(tools::StructureParameter::MoleculeId, getId());
  params.set(tools::StructureParameter::MoleculeType, getType());
  return params;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_MOLECULE_H
