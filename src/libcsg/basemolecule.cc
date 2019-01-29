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

#include <cassert>
#include <votca/csg/basemolecule.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace csg {

void BaseMolecule::AddBead(BaseBead* bead) {
  assert(!beads_.count(bead->getId()) &&
         "Cannot add a bead to the basemolecule"
         " when it has been previously added.");

  BeadStructure::AddBead(bead);
  bead_name_and_ids_[bead->getName()].insert(bead->getId());
  bead->setMolecule(this);
}

const string BaseMolecule::getBeadName(int id) const {
  assert(beads_.count(id) &&
         "Cannot get bead name for bead id because "
         "is is not stored in the base molecule.");
  return beads_.at(id)->getName();
}

const string& BaseMolecule::getBeadType(const int& id) const {
  assert(beads_.count(id) &&
         "Cannot get bead type with id beacuse "
         "bead is not stored in base molecule.");
  return beads_.at(id)->getType();
}

const vec& BaseMolecule::getBeadPosition(const int& id) const {
  assert(beads_.count(id) &&
         "Cannot get bead position with id because "
         "bead is not stored in the base molecule.");
  return beads_.at(id)->getPos();
}

unordered_set<int> BaseMolecule::getBeadIdsByName(const string& name) const {
  assert(bead_name_and_ids_.count(name) &&
         "BaseMolecule does not contain any "
         "beads with name ");
  return bead_name_and_ids_.at(name);
}

}  // namespace csg
}  // namespace votca
