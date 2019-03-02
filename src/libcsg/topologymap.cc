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

  string cg_molecule_type = options.get("cg_molecule.name").as<string>();
  string atomistic_molecule_type =
      options.get("cg_molecule.ident").as<string>();

  AtomAndCGMoleculeTypes molecule_types(atomistic_molecule_type,
                                        cg_molecule_type);

  if (molecule_names_and_maps_.count(molecule_types)) {
    throw invalid_argument(
        "Cannot create map for specified molecules as it "
        "already exists.");
  }

  unordered_map<string, BeadMapInfo> bead_maps_info = ParseBeads_(options);
  ParseMaps_(options, bead_maps_info);

  name_and_molecule_maps_.at(molecule_types) =
      AtomisticToCGMolecaleMapper(atomistic_mol, cg_mol);

  list<Property *> bead_maps_prop = options_maps.Select("map");
  // Cycle maps
  for (list<Property *>::iterator bead_map_prop = bead_maps_prop.begin();
       bead_map_prop != bead_maps_prop.end(); ++bead_map_prop) {
    // maps_[(*bead_map_prop)->get("name").as<string>()] = *bead_map_prop;
  }
  /*
    if ((unsigned int)cg_mol.BeadCount() != beads_.size()) {
      throw runtime_error(
          "number of beads for cg molecule and mapping definition do "
          "not match, check your molecule naming.");
    }
  */
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

unordered_map<string, BeadMapInfo> TopologyMap::ParseBeads_(
    Property &options_in) {
  Property options = options.in.get("cg_beads");
  list<Property *> beads = options.Select("cg_bead");

  unordered_map<string, BeadMapInfo> map_type_and_info;
  for (list<Property *>::iterator bead_iter = beads.begin();
       bead_iter != beads.end(); ++bead_iter) {

    Property *p = *bead_iter;
    BeadMapInfo map_info;
    // map_info.cg_mol_type_ = p->get("name").as<string>();
    map_info.cg_bead_type_ = p->get("type").as<string>();

    // get the beads
    vector<string> subbeads;
    string bead_string(p->get("beads").value());
    Tokenizer tok_beads(bead_string, " \n\t");
    tok_beads.ToVector(subbeads);

    map_info.atomistic_subbeads_ = subbeads;

    map_info.cg_symmetry_ = 1;
    if (p->exists("symmetry")) {
      map_info.symmetry_ = p->get("symmetry").as<int>();
    }

    string map_type = p->get("mapping").as<string>();
    if (map_type_and_info.count(map_type)) {
      throw runtime_error(string("map name ") + map_type +
                          " not unique in mapping");
    }
    map_type_and_info[map_type] = map_info;
  }
  return map_type_and_info;
}

TopologyMap::ParseMaps_(Property &options_in,
                        unordered_map<string, BeadMapInfo> &bead_maps_info) {

  Property maps_prop = options.get("cg_molecule.maps");
  list<Property *> all_maps = maps_prop.Select("map");

  for (list<Property *>::iterator bead_map_iter = all_maps.begin();
       bead_map_iter != all_maps.end(); ++bead_map_iter) {

    Property *p = *bead_map_iter;
    string map_type = p->get("name").as<string>();
    // Ensure that the map is known
    if (bead_maps_info.count(map_type) == 0) {
      throw invalid_argument(
          "Map name " + map_type +
          " is not known because there are no cg_beads using it.");
    }
    // get vector of weights
    Tokenizer tok_weights(opts_map_->get("weights").value(), " \n\t");
    tok_weights.ConvertToVector<double>(weights);

    if (weights.size() != bead_maps_info.cg_subbeads_.size()) {
      string error_msg =
          accumulate(bead_maps_info[map_type].cg_subbeads_.begin(),
                     bead_maps_info[map_type].cg_subbeads_.end(), string(" "));
      throw invalid_argument("The beads in the cg_bead are " + error_msg +
                             " beads the weights are " + tok_weights);
    }

    // get vector of d values used for non-spherical beads
    if (opts_map_->exists("d")) {
      vector<double> d;
      Tokenizer tok_d(opts_map_->get("d").value(), " \n\t");
      tok_d.ConvertToVector(d);
      if (d.size() != bead_maps_info.cg_subbeads_.size()) {
        string error_msg = accumulate(
            bead_maps_info[map_type].cg_subbeads_.begin(),
            bead_maps_info[map_type].cg_subbeads_.end(), string(" "));
        throw invalid_argument("The beads in the cg_bead are " + error_msg +
                               " beads the d values are " + tok_d);
      }
      bead_maps_info[map_type].subbead_d_ = d;
    }

    bead_maps_info[map_type].subbead_weights_ = weights;
  }
}

TopologyMap::~TopologyMap() {
  //  MapContainer::iterator i;

  // for (i = _maps.begin(); i != _maps.end(); ++i) delete *i;
  // _maps.clear();
}

void TopologyMap::Apply(AtomCGConverter converter) {
  MapContainer::iterator iter;

  _out->setStep(_in->getStep());
  _out->setTime(_in->getTime());
  _out->setBox(_in->getBox());

  // Get the cg_beads
  vector<int> cg_molecule_ids = cg_top_->getMoleculeIds();
  for (int &cg_molecule_id : cg_molecule_ids) {
    Molecule *cg_molecule = cg_top_->getMolecule(cg_molecule_id);
    string cg_molecule_type = cg_molecule->getType();
    vector<int> cg_bead_ids = cg_molecule->getBeadIds();
    unordered_map<int, string> bead_ids_and_names =
        converter.MapAtomicBeadIdsToAtomicBeadNames(cg_molecule_type,
                                                    cg_bead_ids);
    map<string, Bead *> name_and_atomic_bead;
    for (pair<int, string> id_and_bead_name : bead_ids_and_names) {
      name_and_atomic_bead[id_and_bead_name.second] =
          atomistic_top_->getBead(id_and_bead_name.first);
    }

    for (int &cg_bead_id : cg_bead_ids) {
      Bead *cg_bead = cg_top_->getBead(cg_bead_id);
    }
  }
  // for (iter = _maps.begin(); iter != _maps.end(); ++iter) (*iter)->Apply();
}

}  // namespace csg
}  // namespace votca
