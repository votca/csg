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

#include "basetopology.h"
#include "bead.h"
#include "molecule.h"

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class CSG_Topology : public Topology<Bead, Molecule> {
 public:
  Molecule *Topology::CreateMolecule(std::string name) {
    Molecule *mol = new Molecule(this, molecules_.size(), name);
    molecules_.push_back(mol);
    return mol;
  }

  Bead *Topology::CreateBead(byte_t symmetry, std::string name,
                             std::string type, int residue_number,
                             std::string residue_name,
                             std::string molecule_name, double m, double q) {

    int molecule_type_id;

    if (molecule_name_and_type_id_.count(molecule_name) == 0) {
      molecule_type_id =
          static_cast<int>(molecule_name_and_type_id_.size()) + 1;
      molecule_name_and_type_id_[molecule_name] = molecule_type_id;
    } else {
      molecule_type_id = molecule_name_and_type_id_[molecule_name];
    }

    std::pair<int, std::string> element{residue_number, residue_name};
    molecule_id_and_residue_id_and_name_[molecule_type_id].insert(element);

    Bead *bead = new Bead(this, beads_.size(), type, symmetry, name,
                          residue_number, residue_name, molecule_type_id, m, q);

    beads_.push_back(bead);
    return bead;
  }
};

}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_CSGTOPOLOGY_H
