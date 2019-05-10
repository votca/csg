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
#ifndef VOTCA_CSG_CSGTOPOLOGY_H
#define VOTCA_CSG_CSGTOPOLOGY_H

#include <votca/tools/structureparameters.h>
#include <votca/tools/types.h>

#include "bead.h"
#include "molecule.h"
#include "templatetopology.h"

namespace votca {
namespace csg {

class CSG_Topology : public TemplateTopology<Bead, Molecule> {
 public:
  CSG_Topology(){};
  ~CSG_Topology(){};

  /// Define rule of 3
  /// Copy constructor
  CSG_Topology(const CSG_Topology& top) { this->Copy(top); }

  CSG_Topology& operator=(const CSG_Topology& top) {
    this->Copy(top);
    return *this;
  }

  Molecule& CreateMolecule(tools::StructureParameters& params) {
    std::string molecule_type =
        params.get<std::string>(tools::StructureParameter::MoleculeType);
    int molecule_id = params.get<int>(tools::StructureParameter::MoleculeId);
    CreateMolecule(molecule_id, molecule_type);
    // std::cout << &(molecules_.back()) << std::endl;
    return *(molecules_map_[molecule_id]);
  }

  Molecule& CreateMolecule(int id, std::string molecule_type) {
    if (!type_container_.MoleculeTypeExist(molecule_type)) {
      type_container_.AddMoleculeType(molecule_type);
    }

    assert(!molecules_map_.count(id) &&
           "molecule with the provided id already exists within the topology!");

    //    Molecule molecule = Molecule(id, molecule_type);
    // size_t index = molecules_.size();
    //    molecules_.push_back(molecule);
    // molecules_.resize(index+1);
    size_t index = molecules_.size();
    molecules_.push_back(Molecule(id, molecule_type));
    molecules_map_[id] = &(molecules_[index]);
    return *(molecules_map_[id]);
  }

  Bead& CreateBead(tools::StructureParameters& params) {
    tools::byte_t symmetry =
        params.get<tools::byte_t>(tools::StructureParameter::Symmetry);
    std::string bead_type =
        params.get<std::string>(tools::StructureParameter::BeadType);
    int bead_id = params.get<int>(tools::StructureParameter::BeadId);
    int molecule_id = params.get<int>(tools::StructureParameter::MoleculeId);
    int residue_id = params.get<int>(tools::StructureParameter::ResidueId);
    std::string residue_type =
        params.get<std::string>(tools::StructureParameter::ResidueType);
    std::string element =
        params.get<std::string>(tools::StructureParameter::Element);
    double mass = params.get<double>(tools::StructureParameter::Mass);
    double charge = params.get<double>(tools::StructureParameter::Charge);

    return CreateBead(symmetry, bead_type, bead_id, molecule_id, residue_id,
                      residue_type, element, mass, charge);
  }
  Bead& CreateBead(tools::byte_t symmetry, std::string bead_type, int bead_id,
                   int molecule_id, int residue_id, std::string residue_type,
                   std::string element_symbol, double mass, double charge) {

    assert(bead_id >= 0 && "Bead id is invalid");
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

    beads_[bead_id] = bead;
    return beads_.at(bead_id);
  }

  Interaction* CreateInteraction(InteractionType type, std::string group,
                                 int bond_id, int molecule_id,
                                 std::vector<int> bead_ids) {
    assert(beads_.size() > 0 &&
           "Cannot create interactions before beads have been initialized");

    for (const int& bead_id : bead_ids) {
      assert(beads_.count(bead_id) &&
             "Cannot add interaction as there are no beads, create the beads "
             "before the interactions.");
    }

    if (type == InteractionType::bond) {
      interactions_.push_back(std::unique_ptr<IBond>(new IBond(bead_ids)));
    } else if (type == InteractionType::angle) {
      interactions_.push_back(std::unique_ptr<IAngle>(new IAngle(bead_ids)));
    } else if (type == InteractionType::dihedral) {
      interactions_.push_back(
          std::unique_ptr<IDihedral>(new IDihedral(bead_ids)));
    } else {
      assert(!"Interaction type is not recognized");
    }
    std::unique_ptr<Interaction>& ic = interactions_.back();
    ic->setGroup(group);
    ic->setIndex(bond_id);
    ic->setMoleculeId(molecule_id);
    // Update the interaction groups
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

#endif  // VOTCA_CSG_CSGTOPOLOGY_H
