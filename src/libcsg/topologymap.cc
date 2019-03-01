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

#include <votca/csg/topologymap.h>

namespace votca {
namespace csg {

void TopologyMap::LoadMap(std::string file_name) {

  Property *options;

  load_property_from_xml(options, file_name);

  Property *options_maps = options.get("cg_molecule.maps");

  list<Property *> molecule_maps = options_maps.Select("map");
  // Cycle maps
  for (list<Property *>::iterator map = molecule_maps.begin();
       molecule_map != molecule_maps.end(); ++molecule_map) {
    maps_[(*molecule_map)->get("name").as<string>()] = *molecule_map;
  }

  if ((unsigned int)cg_mol.BeadCount() != beads_.size()) {
    throw runtime_error(
        "number of beads for cg molecule and mapping definition do "
        "not match, check your molecule naming.");
  }

  AtomisticToCGMolecaleMapper *map =
      new AtomisticToCGMolecaleMapper(atomistic_mol, cg_mol);
  for (vector<CGBeadInfo *>::iterator beaddef = beads_.begin();
       beaddef != beads_.end(); ++beaddef) {

    unordered_set<int> bead_ids = cg_mol.getBeadIdsByType((*beaddef)->type_);
    assert(bead_ids.size() == 1 &&
           "There should only be one bead, if you want a more unique specifier "
           "the beads globally unique id should be used.");

    int bead_id = *bead_ids.begin();
    if (bead_id < 0) {
      throw runtime_error(string("mapping error: reference molecule " +
                                 (*beaddef)->type_ + " does not exist"));
    }

    Property *mdef = getMapByName((*beaddef)->mapping_);
    if (!mdef) {
      throw runtime_error(
          string("mapping " + (*beaddef)->mapping_ + " not found"));
    }

    map->CreateBeadMap((*beaddef)->symmetry_, boundaries, &atomistic_mol,
                       cg_mol.getBead(bead_id), ((*beaddef)->options_), mdef);
  }
  topology_map.AddMoleculeMap(molecule_map);
}

TopologyMap::~TopologyMap() {
  MapContainer::iterator i;

  for (i = _maps.begin(); i != _maps.end(); ++i) delete *i;
  _maps.clear();
}

void TopologyMap::Apply() {
  MapContainer::iterator iter;

  _out->setStep(_in->getStep());
  _out->setTime(_in->getTime());
  _out->setBox(_in->getBox());

  for (iter = _maps.begin(); iter != _maps.end(); ++iter) (*iter)->Apply();
}

}  // namespace csg
}  // namespace votca
