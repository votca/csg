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

/*
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <regex>
#include <stddef.h>
#include <stdexcept>
#include <string>
#include <unordered_set>  // IWYU pragma: keep

#include <votca/csg/boundarycondition.h>
#include <votca/csg/interaction.h>
#include <votca/csg/molecule.h>
#include <votca/csg/openbox.h>
//#include <votca/csg/topology.h>
#include <votca/tools/matrix.h>
#include <votca/tools/rangeparser.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

using namespace std;
bool is_digits(const std::string &str) {
  return str.find_first_not_of("0123456789") == std::string::npos;
}*/
/*
Topology::~Topology() {
  Cleanup();
  if (_bc) delete (_bc);
  _bc = nullptr;
}

void Topology::Cleanup() {
  // cleanup beads
  {
    BeadContainer::iterator i;
    for (i = _beads.begin(); i < _beads.end(); ++i) delete *i;
    _beads.clear();
  }
  // cleanup molecules
  {
    MoleculeContainer::iterator i;
    for (i = _molecules.begin(); i < _molecules.end(); ++i) delete *i;
    _molecules.clear();
  }
  // cleanup residues
  {
    molecule_id_and_residue_id_and_name_.clear();
    molecule_name_and_type_id_.clear();
  }
  // cleanup interactions
  {
    InteractionContainer::iterator i;
    for (i = _interactions.begin(); i < _interactions.end(); ++i) delete (*i);
    _interactions.clear();
  }
  // cleanup _bc object
  if (_bc) delete (_bc);
  _bc = new OpenBox();
}*/

/// \todo implement checking, only used in xml topology reader
/*void Topology::CreateMoleculesByRange(string name, int first, int nbeads,
                                      int nmolecules) {
  Molecule *mol = CreateMolecule(name);
  int beadcount = 0;

  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    // xml numbering starts with 1
    if (--first > 0) continue;
    // This is not 100% correct, but let's assume for now that the resnr do
    // increase
    mol->AddBead((*bead));

    if (++beadcount == nbeads) {
      if (--nmolecules <= 0) break;
      mol = CreateMolecule(name);
      beadcount = 0;
    }
  }
}*/
/*
void Topology::CopyTopologyData(Topology *top) {
  BeadContainer::iterator it_bead;
  MoleculeContainer::iterator it_mol;

  _bc->setBox(top->getBox());
  _time = top->_time;
  _step = top->_step;

  // cleanup old data
  Cleanup();

  // Copy residue info
  setResidueIdsAndNames(top->getResidueIdsAndNames());
  setMoleculeNamesAndIds(top->getMoleculeNamesAndIds());

  // create all beads
  for (it_bead = top->_beads.begin(); it_bead != top->_beads.end(); ++it_bead) {
    Bead *bi = *it_bead;
    string type = bi->getType();
    string molecule_name = top->getMolecule(bi->getMoleculeId())->getName();
    CreateBead<Bead>(bi->getSymmetry(), bi->getName(), type,
                     bi->getResidueNumber(), bi->getResidueName(),
                     molecule_name, bi->getMass(), bi->getQ());
  }

  // copy all molecules
  for (it_mol = top->_molecules.begin(); it_mol != top->_molecules.end();
       ++it_mol) {
    Molecule *mi = CreateMolecule((*it_mol)->getName());
    vector<int> bead_ids = (*it_mol)->getBeadIds();
    for (const int &bead_id : bead_ids) {
      mi->AddBead(_beads[bead_id]);
    }
  }
}
*/
/*int Topology::getBeadTypeId(string type) const {
  assert(beadtypes_.count(type));
  return beadtypes_.at(type);
}*/
/*
void Topology::RenameMolecules(string range, string name) {
  RangeParser rp;
  RangeParser::iterator i;

  rp.Parse(range);
  for (i = rp.begin(); i != rp.end(); ++i) {
    if ((unsigned int)*i > _molecules.size()) {
      throw runtime_error(
          string("RenameMolecules: num molecules smaller than"));
    }
    getMolecule(*i - 1)->setName(name);
  }
}*/
/*
void Topology::RenameBeadType(string name, string newname) {
  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    string type = (*bead)->getType();
    if (wildcmp(name.c_str(), type.c_str())) {
      (*bead)->setType(newname);
    }
  }
}*/
/*
void Topology::SetBeadTypeMass(string name, double value) {
  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    string type = (*bead)->getType();
    if (wildcmp(name.c_str(), type.c_str())) {
      (*bead)->setMass(value);
    }
  }
}
*/
/*void Topology::CheckMoleculeNaming(void) {
  map<string, int> nbeads;

  for (MoleculeContainer::iterator iter = _molecules.begin();
       iter != _molecules.end(); ++iter) {
    map<string, int>::iterator entry = nbeads.find((*iter)->getName());
    if (entry != nbeads.end()) {
      if (entry->second != static_cast<int>((*iter)->BeadCount()))
        throw runtime_error(
            "There are molecules which have the same name but different number "
            "of bead please check the section manual topology handling in the "
            "votca manual");
      continue;
    }
    nbeads[(*iter)->getName()] = (*iter)->BeadCount();
  }
}*/
/*
void Topology::AddBondedInteraction(Interaction *ic) {
  map<string, int>::iterator iter;
  iter = _interaction_groups.find(ic->getGroup());
  if (iter != _interaction_groups.end()) {
    ic->setGroupId((*iter).second);
  } else {
    int i = _interaction_groups.size();
    _interaction_groups[ic->getGroup()] = i;
    ic->setGroupId(i);
  }
  _interactions.push_back(ic);
  _interactions_by_group[ic->getGroup()].push_back(ic);
}*/
/*
std::list<Interaction *> Topology::InteractionsInGroup(const string &group) {
  map<string, list<Interaction *> >::iterator iter;
  iter = _interactions_by_group.find(group);
  if (iter == _interactions_by_group.end()) return list<Interaction *>();
  return iter->second;
}
*/
/*bool Topology::BeadTypeExist(string type) const {
  return beadtypes_.count(type);
}

void Topology::RegisterBeadType(string type) {
  unordered_set<int> ids;
  for (pair<const string, int> type_and_id : beadtypes_) {
    ids.insert(type_and_id.second);
  }

  int id = 0;
  // If the type is also a number use it as the id as well provided it is not
  // already taken
  if (is_digits(type)) {
    id = boost::lexical_cast<int>(type);
    assert(!ids.count(id) &&
           "The type passed in is a number and has already"
           " been registered. It is likely that you are passing in numbers as "
           "bead types as well as strings, choose one or the other do not mix "
           "between using numbers and strings ");
  }

  while (ids.count(id)) {
    ++id;
  }
  beadtypes_[type] = id;
}
*/
/*vec Topology::BCShortestConnection(const vec &r_i, const vec &r_j) const {
  return _bc->BCShortestConnection(r_i, r_j);
}

vec Topology::getDist(int bead1, int bead2) const {
  return BCShortestConnection(getBead(bead1)->getPos(),
                              getBead(bead2)->getPos());
}
*/
/*double Topology::BoxVolume() { return _bc->BoxVolume(); }

void Topology::RebuildExclusions() { _exclusions.CreateExclusions(this); }

BoundaryCondition::eBoxtype Topology::autoDetectBoxType(const matrix &box) {
  // set the box type to OpenBox in case "box" is the zero matrix,
  // to OrthorhombicBox in case "box" is a diagonal matrix,
  // or to TriclinicBox otherwise
  if (box.get(0, 0) == 0 && box.get(0, 1) == 0 && box.get(0, 2) == 0 &&
      box.get(1, 0) == 0 && box.get(1, 1) == 0 && box.get(1, 2) == 0 &&
      box.get(2, 0) == 0 && box.get(2, 1) == 0 && box.get(2, 2) == 0) {
    // cout << "box open\n";
    return BoundaryCondition::typeOpen;
  } else if (box.get(0, 1) == 0 && box.get(0, 2) == 0 && box.get(1, 0) == 0 &&
             box.get(1, 2) == 0 && box.get(2, 0) == 0 && box.get(2, 1) == 0) {
    // cout << "box orth\n";
    return BoundaryCondition::typeOrthorhombic;
  } else {
    // cout << "box tric\n";
    return BoundaryCondition::typeTriclinic;
  }
  return BoundaryCondition::typeOpen;
}*/
/*
double Topology::ShortestBoxSize() {
  vec _box_a = getBox().getCol(0);
  vec _box_b = getBox().getCol(1);
  vec _box_c = getBox().getCol(2);

  // create plane normals
  vec _norm_a = _box_b ^ _box_c;
  vec _norm_b = _box_c ^ _box_a;
  vec _norm_c = _box_a ^ _box_b;

  _norm_a.normalize();
  _norm_b.normalize();
  _norm_c.normalize();

  double la = _box_a * _norm_a;
  double lb = _box_b * _norm_b;
  double lc = _box_c * _norm_c;

  return min(la, min(lb, lc));
}

}  // namespace csg
} */ // namespace votca
