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

#ifndef _VOTCA_CSG_CSGTOPOLOGY_H
#define _VOTCA_CSG_CSGTOPOLOGY_H

#include <votca/tools/types.h>

#include "basetopology.h"
#include "bead.h"
#include "molecule.h"

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class CSG_Topology : public Topology<Bead, Molecule> {
 public:
  ~CSG_Topology(){};
  Molecule* CreateMolecule(int id, std::string molecule_type) {
    if (!type_container_.MoleculeTypeExist(molecule_type)) {
      type_container_.AddMoleculeType(molecule_type);
    }

    assert(!molecules_.count(id) &&
           "molecule with the provided id already exists within the topology!");

    Molecule molecule = Molecule(id, molecule_type);
    molecules_[id] = molecule;
    return &molecules_[id];
  }

  Bead* CreateBead(byte_t symmetry, std::string bead_type, int bead_id,
                   int molecule_id, int residue_id, std::string residue_type,
                   std::string element_symbol, double mass, double charge) {

    if (!type_container_.ResidueTypeExist(residue_type)) {
      type_container_.AddResidueType(residue_type);
    }
    if (!type_container_.BeadTypeExist(bead_type)) {
      type_container_.AddBeadType(bead_type);
    }

    assert(!beads_.count(bead_id) &&
           "bead with provided id already exists in the topology!");

    Bead bead = Bead(symmetry, bead_id, bead_type, residue_id, residue_type,
                     molecule_id, element_symbol, mass, charge);

    std::cout << "Creating bead " << bead_id << std::endl;
    beads_[bead_id] = bead;
    return &beads_.at(bead_id);
  }

  Interaction* CreateInteraction(Interaction::interaction_type type,
                                 std::string group, int bond_id,
                                 int molecule_id, std::vector<int> bead_ids) {
    assert(beads_.size() > 0 &&
           "Cannot create interactions before beads have been initialized");

    std::vector<const BaseBead*> beads;
    for (const int& bead_id : bead_ids) {
      std::cout << "Bead id " << bead_id << std::endl;
      assert(beads_.count(bead_id) &&
             "Cannot add interaction as there are no beads, create the beads "
             "before the interactions.");
      beads.push_back(&beads_.at(bead_id));
    }

    if (type == Interaction::bond) {
      interactions_.push_back(std::unique_ptr<IBond>(new IBond(beads)));
    } else if (type == Interaction::angle) {
      interactions_.push_back(std::unique_ptr<IAngle>(new IAngle(beads)));
    } else if (type == Interaction::dihedral) {
      interactions_.push_back(std::unique_ptr<IDihedral>(new IDihedral(beads)));
    } else {
      assert(!"Interaction type is not recognized");
    }
    std::unique_ptr<Interaction>& ic = interactions_.back();
    ic->setGroup(group);
    ic->setIndex(bond_id);
    ic->setMoleculeId(molecule_id);
    // Update the interaction groups
    std::cout << "Updateing interaction groups " << std::endl;
    std::map<std::string, int>::iterator iter;
    iter = interaction_groups_.find(ic->getGroup());
    if (iter != interaction_groups_.end()) {
      ic->setGroupId((*iter).second);
    } else {
      int group_size = interaction_groups_.size();
      interaction_groups_[ic->getGroup()] = group_size;
      ic->setGroupId(group_size);
    }
    interactions_by_group_[ic->getGroup()].push_back(ic.get());

    return interactions_.back().get();
  }
};

}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_CSGTOPOLOGY_H
