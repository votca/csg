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

#ifndef _VOTCA_CSG_MOLECULE_H
#define _VOTCA_CSG_MOLECULE_H

#include "basemolecule.h"
#include "bead.h"
#include "topologyitem.h"
#include <assert.h>
#include <map>
#include <string>
#include <vector>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class Interaction;

namespace molecule_constants {
const std::string molecule_name_unassigned = "unassigned";
}
/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule : public TopologyItem, public BaseMolecule<Bead> {
 public:
  /// get the molecule ID
  // int getId() const { return _id; }

  /// get the name of the molecule
  // const std::string &getName() const { return _name; }

  /// set the name of the molecule
  // void setName(const std::string &name) {  _name=name; }

  /// Add a bead to the molecule
  // void AddBead(Bead *bead, const std::string &name);
  /// get the id of a bead in the molecule
  // Bead *getBead(int bead) { return _beads[bead]; }
  // int getBeadId(int bead) { return _beads[bead]->getId(); }
  // Switched to getBeadIdsByName
  // int getBeadIdByName(const std::string &name);

  /// get the number of beads in the molecule
  // int BeadCount() const { return _beads.size(); }

  /// find a bead by it's name
  // Switched again to getBeadIdsByName
  // int getBeadByName(const std::string &name);
  // std::string getBeadName(int bead) {return _bead_names[bead]; }

  std::unordered_set<int> getBeadIdsByLabel(const std::string &label);

  /// Add an interaction to the molecule
  void AddInteraction(Interaction *ic) { _interactions.push_back(ic); }

  std::vector<Interaction *> Interactions() { return _interactions; }

  template <typename T>
  void setUserData(T *userdata) {
    _userdata = (void *)userdata;
  }

  template <typename T>
  T *getUserData() {
    return (T *)_userdata;
  }

 private:
  // maps a name to a bead id
  // std::map<std::string, int> _beadmap;
  std::vector<Interaction *> _interactions;

  // id of the molecules
  // int _id;

  // name of the molecule
  // std::string _name;
  // the beads in the molecule
  // std::vector<Bead *> _beads;
  // std::vector<std::string> _bead_names;

  void *_userdata;

  /// constructor
  Molecule(Topology *parent, int id, std::string name) : TopologyItem(parent) {
    id_.setId(id);
    name_.setName(name);
  }

  friend class Topology;
};

inline std::unordered_set<int> Molecule::getBeadIdsByLabel(
    const std::string &label) {
  std::unordered_set<int> bead_ids;
  for (const std::pair<const int, Bead *> &id_and_bead : beads_) {
    std::cout << "Lable of bead " << id_and_bead.second->getLabel()
              << std::endl;
    if (label.compare(id_and_bead.second->getLabel()) == 0) {
      bead_ids.insert(id_and_bead.first);
    }
  }
  return bead_ids;
}

/*
inline int Molecule::getBeadIdByName(const std::string &name)
{
    int i = getBeadByName(name);
    if(i<0)
        return i;
    return _beads[i]->getId();
}*/

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_MOLECULE_H */
