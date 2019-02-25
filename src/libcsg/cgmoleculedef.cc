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

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <stddef.h>
#include <stdexcept>
#include <string>

#include "../../include/votca/csg/csgtopology.h"
#include <votca/csg/bead.h>
#include <votca/csg/cgmoleculedef.h>
#include <votca/csg/interaction.h>
#include <votca/csg/map.h>

#include <votca/tools/constants.h>
#include <votca/tools/property.h>
#include <votca/tools/tokenizer.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {
class Molecule;
class Residue;
}  // namespace csg
}  // namespace votca

namespace votca {
namespace csg {

using boost::lexical_cast;

CGMoleculeDef::~CGMoleculeDef() {
  {
    vector<beaddef_t *>::iterator i;
    for (i = beads_.begin(); i != beads_.end(); ++i) {
      delete *i;
    }
    beads_.clear();
  }
}

void CGMoleculeDef::Load(string filename) {
  load_property_from_xml(options_, filename);
  // parse xml tree
  cg_molecule_type_ = options_.get("cg_molecule.name").as<string>();
  ident_ = options_.get("cg_molecule.ident").as<string>();

  ParseTopology(options_.get("cg_molecule.topology"));
  ParseMapping(options_.get("cg_molecule.maps"));
}

void CGMoleculeDef::ParseTopology(Property &options) {
  ParseBeads(options.get("cg_beads"));
  if (options.exists("cg_bonded")) ParseBonded(options.get("cg_bonded"));
}

void CGMoleculeDef::ParseBeads(Property &options) {
  list<Property *> beads = options.Select("cg_bead");

  for (list<Property *>::iterator iter = beads.begin(); iter != beads.end();
       ++iter) {
    Property *p = *iter;
    beaddef_t *beaddef = new beaddef_t;
    beaddef->options_ = p;

    beaddef->cg_name_ = p->get("name").as<string>();
    beaddef->type_ = p->get("type").as<string>();
    beaddef->mapping_ = p->get("mapping").as<string>();
    if (p->exists("symmetry")) {
      beaddef->symmetry_ = p->get("symmetry").as<int>();
    } else {
      beaddef->symmetry_ = 1;
    }
    if (beads_by_cg_name_.find(beaddef->cg_name_) != beads_by_cg_name_.end()) {
      throw std::runtime_error(string("bead name ") + beaddef->cg_name_ +
                               " not unique in mapping");
    }
    beads_.push_back(beaddef);
    beads_by_cg_name_[beaddef->cg_name_] = beaddef;
  }
}

void CGMoleculeDef::ParseBonded(Property &options) {
  bonded_ = options.Select("*");
}

void CGMoleculeDef::ParseMapping(Property &options) {
  list<Property *> maps = options.Select("map");

  for (list<Property *>::iterator iter = maps.begin(); iter != maps.end();
       ++iter) {
    maps_[(*iter)->get("name").as<string>()] = *iter;
  }
}

Molecule *CGMoleculeDef::CreateMolecule(CSG_Topology &top) {
  int molecule_id = top.MoleculeCount();
  Molecule *minfo = top.CreateMolecule(molecule_id, cg_molecule_type_);

  // create the atoms
  // vector<beaddef_t *>::iterator beaddef_
  // for (beaddef = beads_.begin(); beaddef != beads_.end(); ++beaddef) {
  for (beaddef_t *beaddef : beads_) {
    Bead *bead;

    string bead_type = beaddef->type_;
    bead =
        top.CreateBead(beaddef->symmetry_, bead_type, top.BeadCount(),
                       molecule_id, topology_constants::unassigned_residue_id,
                       topology_constants::unassigned_residue_type,
                       topology_constants::unassigned_element, 0.0, 0.0);
    minfo->AddBead(bead);
  }

  // create the bonds
  list<Property *>::iterator ibnd;
  map<string, string> had_iagroup;

  for (ibnd = bonded_.begin(); ibnd != bonded_.end(); ++ibnd) {
    vector<int> atoms;
    string iagroup = (*ibnd)->get("name").as<string>();

    if (had_iagroup[iagroup] == "yes") {
      throw runtime_error(
          string("double occurence of interactions with name ") + iagroup);
    }
    had_iagroup[iagroup] = "yes";

    Tokenizer tok((*ibnd)->get("beads").value(), " \n\t");
    for (Tokenizer::iterator atom = tok.begin(); atom != tok.end(); ++atom) {
      unordered_set<int> bead_ids = minfo->getBeadIdsByType(*atom);
      assert(bead_ids.size() == 1 &&
             "There is more than one bead with that type if you want a unique "
             "identifier you should probably just use the beads global unique "
             "id.");

      int bead_id = *bead_ids.begin();
      if (bead_id < 0) {
        throw runtime_error(
            string("error while trying to create bonded interaction, bead " +
                   *atom + " not found"));
      }
      atoms.push_back(bead_id);
    }

    int NrBeads = 1;
    if ((*ibnd)->name() == "bond") {
      NrBeads = 2;
    } else if ((*ibnd)->name() == "angle") {
      NrBeads = 3;
    } else if ((*ibnd)->name() == "dihedral") {
      NrBeads = 4;
    }

    if ((atoms.size() % NrBeads) != 0) {
      throw runtime_error("Number of atoms in interaction '" +
                          (*ibnd)->get("name").as<string>() +
                          "' is not a multiple of " +
                          lexical_cast<string>(NrBeads) + "! Missing beads?");
    }
    int index = 0;
    while (!atoms.empty()) {
      Interaction *ic;

      if ((*ibnd)->name() == "bond") {
        ic = top.CreateInteraction(Interaction::interaction_type::bond, iagroup,
                                   index, minfo->getId(), atoms);
      } else if ((*ibnd)->name() == "angle") {
        ic = top.CreateInteraction(Interaction::interaction_type::angle,
                                   iagroup, index, minfo->getId(), atoms);
      } else if ((*ibnd)->name() == "dihedral") {
        ic = top.CreateInteraction(Interaction::interaction_type::dihedral,
                                   iagroup, index, minfo->getId(), atoms);
      } else {
        throw runtime_error("unknown bonded type in map: " + (*ibnd)->name());
      }

      minfo->AddInteraction(ic);
      index++;
    }
  }
  return minfo;
}

Map *CGMoleculeDef::CreateMap(const CSG_Topology *topology,
                              const Molecule &mol_in, Molecule &mol_out) {
  if ((unsigned int)mol_out.BeadCount() != beads_.size()) {
    throw runtime_error(
        "number of beads for cg molecule and mapping definition do "
        "not match, check your molecule naming.");
  }

  Map *map = new Map(mol_in, mol_out);
  for (vector<beaddef_t *>::iterator beaddef = beads_.begin();
       beaddef != beads_.end(); ++beaddef) {

    unordered_set<int> bead_ids = mol_out.getBeadIdsByType((*beaddef)->type_);
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

    BeadMap *bmap;
    switch ((*beaddef)->symmetry_) {
      case 1:
        bmap = new Map_Sphere();
        break;
      case 3:
        bmap = new Map_Ellipsoid();
        break;
      default:
        throw runtime_error(string("unknown symmetry in bead definition!"));
    }
    ////////////////////////////////////////////////////
    bmap->Initialize(topology, &mol_in, mol_out.getBead(bead_id),
                     ((*beaddef)->options_), mdef);
    map->AddBeadMap(bmap);
  }
  return map;
}

CGMoleculeDef::beaddef_t *CGMoleculeDef::getBeadByCGName(
    const string &cg_name) {
  map<string, beaddef_t *>::iterator cg_name_and_beaddef =
      beads_by_cg_name_.find(cg_name);
  if (cg_name_and_beaddef == beads_by_cg_name_.end()) {
    std::cout << "cannot find: <" << cg_name << "> in " << cg_molecule_type_
              << "\n";
    return nullptr;
  }
  return (*cg_name_and_beaddef).second;
}

Property *CGMoleculeDef::getMapByName(const string &mapping_name) {
  map<string, Property *>::iterator iter = maps_.find(mapping_name);
  if (iter == maps_.end()) {
    std::cout << "cannot find map " << mapping_name << "\n";
    return NULL;
  }
  return (*iter).second;
}

}  // namespace csg
}  // namespace votca
