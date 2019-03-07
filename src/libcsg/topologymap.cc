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
#include <numeric>
#include <votca/csg/topologymap.h>
#include <votca/tools/property.h>
#include "../../include/votca/csg/map.h"

using namespace std;
namespace votca {
namespace csg {

void TopologyMap::LoadMap(const std::string& file_name) {

  Property options;

  load_property_from_xml(options, file_name);

  string cg_molecule_type = options.get("cg_molecule.name").as<string>();
  string atomistic_molecule_type =
      options.get("cg_molecule.ident").as<string>();

  // Be sure that the atomistic and cg topologies have the molecule names that
  // are to be mapped from and too
  if(cg_top_->MoleculeTypeExist(cg_molecule_type)==false){
    throw runtime_error("cg molecule type does not exist within the coarsegrained topology");
  }
  if(atomistic_top_->MoleculeTypeExist(atomistic_molecule_type)==false){
    throw runtime_error("atomistic molecule type does not exist within the atomic topology");
  }

  AtomAndCGMoleculeTypes molecule_types(atomistic_molecule_type,
                                        cg_molecule_type);

  if (molecule_names_and_maps_.count(molecule_types)) {
    throw invalid_argument(
        "Cannot create map for specified molecules as it "
        "already exists.");
  }

  unordered_map<string, BeadMapInfo> bead_maps_info = ParseBeads_(options);
  ParseMaps_(options, bead_maps_info);

  molecule_names_and_maps_.at(molecule_types) =
      AtomisticToCGMoleculeMapper(atomistic_molecule_type, cg_molecule_type);

  //list<Property *> bead_maps_prop = options.Select("map");
  // Cycle maps
  //for (list<Property *>::iterator bead_map_prop = bead_maps_prop.begin();
  //     bead_map_prop != bead_maps_prop.end(); ++bead_map_prop) {
    // maps_[(*bead_map_prop)->get("name").as<string>()] = *bead_map_prop;
 // }
  /*
    if ((unsigned int)cg_mol.BeadCount() != beads_.size()) {
      throw runtime_error(
          "number of beads for cg molecule and mapping definition do "
          "not match, check your molecule naming.");
    }
  */
//  for (vector<CGBeadInfo *>::iterator beaddef = beads_.begin();
 //      beaddef != beads_.end(); ++beaddef) { 
 
 // Consistency check: Ensure that the coarse grained beads exist in the coarse grained topology before initializing the bead maps
    for( pair<const string, BeadMapInfo> & bead_and_map : bead_maps_info ){
    //unordered_set<int> bead_ids = cg_mol.getBeadIdsByType((*beaddef)->type_);
    //assert(bead_ids.size() == 1 &&
     //      "There should only be one bead, if you want a more unique specifier "
      //     "the beads globally unique id should be used.");

    //int cg_bead_id = *bead_ids.begin();
    string cg_bead_type = bead_and_map.second.cg_bead_type_;


    //if (bead_id < 0) {
     // throw runtime_error(string("mapping error: atomic reference molecule type " +
      //                           (*beaddef)->type_ + " does not exist"));
   // }
    if(cg_top_->BeadTypeExist(cg_bead_type)==false){
      throw runtime_error("Cannot create map to the bead, as the cg bead type does not exist in the coarse grained topology");
    }
  /*  Property *mdef = getMapByName((*beaddef)->mapping_);
    if (!mdef) {
      throw runtime_error(
          string("mapping " + (*beaddef)->mapping_ + " not found"));
    }*/

//    map->CreateBeadMap((*beaddef)->symmetry_, boundaries, &atomistic_mol,
 //                      cg_mol.getBead(bead_id), ((*beaddef)->options_), mdef);
  }
    molecule_names_and_maps_.at(molecule_types).Initialize( bead_maps_info );
  //topology_map.AddMoleculeMap(molecule_map);
}

unordered_map<string, BeadMapInfo> TopologyMap::ParseBeads_(
    Property &options_in) {
  Property options = options_in.get("cg_beads");
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

    map_info.atomic_subbeads_ = subbeads;

    map_info.cg_symmetry_ = 1;
    if (p->exists("symmetry")) {
      map_info.cg_symmetry_ = p->get("symmetry").as<int>();
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

void TopologyMap::ParseMaps_(Property &options_in,
                        unordered_map<string, BeadMapInfo> &bead_maps_info) {

  Property maps_prop = options_in.get("cg_molecule.maps");
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
    vector<double> weights;
    Tokenizer tok_weights(p->get("weights").value(), " \n\t");
    tok_weights.ConvertToVector<double>(weights);

    if (weights.size() != bead_maps_info[map_type].atomic_subbeads_.size()) {
      string error_msg =
          accumulate(bead_maps_info[map_type].atomic_subbeads_.begin(),
                     bead_maps_info[map_type].atomic_subbeads_.end(), string(" "));
      throw invalid_argument(string("The beads in the cg_bead are ") + error_msg +
                             string(" beads the weights are ") + p->get("weights").value());
    }

    // get vector of d values used for non-spherical beads
    if (p->exists("d")) {
      vector<double> d;
      Tokenizer tok_d(p->get("d").value(), " \n\t");
      tok_d.ConvertToVector(d);
      if (d.size() != bead_maps_info[map_type].atomic_subbeads_.size()) {
        string error_msg = accumulate(
            bead_maps_info[map_type].atomic_subbeads_.begin(),
            bead_maps_info[map_type].atomic_subbeads_.end(), string(" "));
        throw invalid_argument(string("The beads in the cg_bead are ") + error_msg +
                               string(" beads the d values are ") +p->get("d").value());
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
 // MapContainer::iterator iter;

  cg_top_->setStep(atomistic_top_->getStep());
  cg_top_->setTime(atomistic_top_->getTime());
  cg_top_->setBox(atomistic_top_->getBox());

  // Get the cg molecules
  vector<int> cg_molecule_ids = cg_top_->getMoleculeIds();
  for (int &cg_molecule_id : cg_molecule_ids) {
    Molecule *cg_molecule = cg_top_->getMolecule(cg_molecule_id);
    string cg_molecule_type = cg_molecule->getType();
    int molecule_id = cg_molecule->getId();
    vector<int> cg_bead_ids = cg_molecule->getBeadIds();
    unordered_map<int, string> bead_ids_and_names =
        converter.MapAtomicBeadIdsToAtomicBeadNames(cg_molecule_type,
                                                    cg_bead_ids);
    map<string, Bead *> name_and_atomic_bead;
    for (pair<int, string> id_and_bead_name : bead_ids_and_names) {
      name_and_atomic_bead[id_and_bead_name.second] =
          atomistic_top_->getBead(id_and_bead_name.first);
    }

    string atomistic_molecule_type = atomistic_top_->getMolecule(molecule_id)->getType();

    AtomAndCGMoleculeTypes atom_and_cg_type;
    atom_and_cg_type.first = atomistic_molecule_type;
    atom_and_cg_type.second = cg_molecule_type;

    molecule_names_and_maps_[atom_and_cg_type].Apply(converter);
    //for (int &cg_bead_id : cg_bead_ids) {
     // Bead *cg_bead = cg_top_->getBead(cg_bead_id);
    //}
  }
  // for (iter = _maps.begin(); iter != _maps.end(); ++iter) (*iter)->Apply();
}

}  // namespace csg
}  // namespace votca
