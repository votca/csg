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
#include <boost/tokenizer.hpp>
#include <iostream>
#include <votca/csg/bead.h>
#include <votca/csg/cgmoleculedef.h>
#include <votca/csg/interaction.h>
#include <votca/csg/topology.h>
#include <votca/tools/tokenizer.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

using boost::lexical_cast;
CGMoleculeDef::~CGMoleculeDef() {
  {
    vector<beaddef_t *>::iterator i;
    for (i = _beads.begin(); i != _beads.end(); ++i) delete *i;
    _beads.clear();
  }
}

void CGMoleculeDef::Load(string filename) {
  load_property_from_xml(_options, filename);
  // parse xml tree
  _name = _options.get("cg_molecule.name").as<string>();
  _ident = _options.get("cg_molecule.ident").as<string>();

  ParseTopology(_options.get("cg_molecule.topology"));
  ParseMapping(_options.get("cg_molecule.maps"));
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
    beaddef->_options = p;

    beaddef->residue_number_ = bead_constants::residue_number_unassigned;
    beaddef->_name = p->get("name").as<string>();
    beaddef->_type = p->get("type").as<string>();
    beaddef->_mapping = p->get("mapping").as<string>();
    if (p->exists("symmetry"))
      beaddef->_symmetry = p->get("symmetry").as<int>();
    else
      beaddef->_symmetry = 1;

    if (_beads_by_name.find(beaddef->_name) != _beads_by_name.end())
      throw std::runtime_error(string("bead name ") + beaddef->_name +
                               " not unique in mapping");
    _beads.push_back(beaddef);
    _beads_by_name[beaddef->_name] = beaddef;
  }
}

void CGMoleculeDef::ParseBonded(Property &options) {
  _bonded = options.Select("*");
}

void CGMoleculeDef::ParseMapping(Property &options) {
  list<Property *> maps = options.Select("map");

  for (list<Property *>::iterator iter = maps.begin(); iter != maps.end();
       ++iter)
    _maps[(*iter)->get("name").as<string>()] = *iter;
}

Molecule *CGMoleculeDef::CreateMolecule(Topology &top) {
  // add the residue names
  // Residue *res = top.CreateResidue(_name);
  Molecule *minfo = top.CreateMolecule(_name);

  // create the atoms
  vector<beaddef_t *>::iterator iter;
  for (iter = _beads.begin(); iter != _beads.end(); ++iter) {
    Bead *bead;

    string type = (*iter)->_type;
    if (!top.BeadTypeExist(type)) {
      top.RegisterBeadType(type);
    }
    bead = top.CreateBead<Bead>((*iter)->_symmetry, (*iter)->_name, type,
                                (*iter)->residue_number_, _name, 0, 0);
    // minfo->AddBead(bead, bead->getName());
    minfo->AddBead(bead);

    bead->setOptions(*(*iter)->_options);
  }

  // create the bonds
  list<Property *>::iterator ibnd;
  map<string, string> had_iagroup;

  for (ibnd = _bonded.begin(); ibnd != _bonded.end(); ++ibnd) {
    list<int> atoms;
    string iagroup = (*ibnd)->get("name").as<string>();

    if (had_iagroup[iagroup] == "yes")
      throw runtime_error(
          string("double occurence of interactions with name ") + iagroup);
    had_iagroup[iagroup] = "yes";

    Tokenizer tok((*ibnd)->get("beads").value(), " \n\t");
    for (Tokenizer::iterator atom = tok.begin(); atom != tok.end(); ++atom) {
      // int i = minfo->getBeadIdByName(*atom);
      unordered_set<int> bead_ids = minfo->getBeadIdsByName(*atom);
      assert(bead_ids.size() == 1 &&
             "There is more than one bead with that name "
             "if you want a unique identifier you should probably just use the "
             "beads global unique id.");
      int bead_id = *bead_ids.begin();
      if (bead_id < 0)
        throw runtime_error(
            string("error while trying to create bonded interaction, "
                   "bead " +
                   *atom + " not found"));

      atoms.push_back(bead_id);
    }

    int NrBeads = 1;
    if ((*ibnd)->name() == "bond")
      NrBeads = 2;
    else if ((*ibnd)->name() == "angle")
      NrBeads = 3;
    else if ((*ibnd)->name() == "dihedral")
      NrBeads = 4;

    if ((atoms.size() % NrBeads) != 0)
      throw runtime_error("Number of atoms in interaction '" +
                          (*ibnd)->get("name").as<string>() +
                          "' is not a multiple of " +
                          lexical_cast<string>(NrBeads) + "! Missing beads?");

    int index = 0;
    while (!atoms.empty()) {
      Interaction *ic;

      if ((*ibnd)->name() == "bond")
        ic = new IBond(atoms);
      else if ((*ibnd)->name() == "angle")
        ic = new IAngle(atoms);
      else if ((*ibnd)->name() == "dihedral")
        ic = new IDihedral(atoms);
      else
        throw runtime_error("unknown bonded type in map: " + (*ibnd)->name());

      ic->setGroup(iagroup);
      ic->setIndex(index);
      ic->setMolecule(minfo->getId());
      top.AddBondedInteraction(ic);
      minfo->AddInteraction(ic);
      index++;
    }
  }
  return minfo;
}

Map *CGMoleculeDef::CreateMap(Molecule &in, Molecule &out) {
  if ((unsigned int)out.BeadCount() != _beads.size()) {
    throw runtime_error(
        "number of beads for cg molecule and mapping definition do "
        "not match, check your molecule naming.");
  }

  Map *map = new Map(in, out);
  for (vector<beaddef_t *>::iterator def = _beads.begin(); def != _beads.end();
       ++def) {

    cout << (*def)->_name << endl;
    // int iout = out.getBeadByName((*def)->_name);
    unordered_set<int> bead_ids = out.getBeadIdsByName((*def)->_name);
    assert(
        bead_ids.size() == 1 &&
        "There should only be one bead, if you want a "
        "more unique specifier the beads globally unique id should be used.");
    int bead_id = *bead_ids.begin();
    if (bead_id < 0)
      throw runtime_error(string("mapping error: reference molecule " +
                                 (*def)->_name + " does not exist"));

    Property *mdef = getMapByName((*def)->_mapping);
    if (!mdef)
      throw runtime_error(string("mapping " + (*def)->_mapping + " not found"));

    /// TODO: change this to factory, do not hardcode!!
    BeadMap *bmap;
    switch ((*def)->_symmetry) {
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

    bmap->Initialize(&in, out.getBead(bead_id), ((*def)->_options), mdef);
    map->AddBeadMap(bmap);
  }
  return map;
}

CGMoleculeDef::beaddef_t *CGMoleculeDef::getBeadByName(const string &name) {
  map<string, beaddef_t *>::iterator iter = _beads_by_name.find(name);
  if (iter == _beads_by_name.end()) {
    std::cout << "cannot find: <" << name << "> in " << _name << "\n";
    return NULL;
  }
  // assert(iter != _beadmap.end());
  // return (*iter).second;
  return (*iter).second;
}

Property *CGMoleculeDef::getMapByName(const string &name) {
  map<string, Property *>::iterator iter = _maps.find(name);
  if (iter == _maps.end()) {
    std::cout << "cannot find map " << name << "\n";
    return NULL;
  }
  // assert(iter != _beadmap.end());
  // return (*iter).second;
  return (*iter).second;
}

}  // namespace csg
}  // namespace votca
