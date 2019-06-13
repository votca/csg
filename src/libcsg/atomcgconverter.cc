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

#include "../../include/votca/csg/atomcgconverter.h"
#include "../../include/votca/csg/bead.h"
#include "../../include/votca/csg/beadmap.h"
#include "../../include/votca/csg/interaction.h"
#include <numeric>
#include <stddef.h>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <votca/tools/constants.h>
#include <votca/tools/property.h>
#include <votca/tools/tokenizer.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/*****************************************************************************
 * Public Facing Methods
 *****************************************************************************/
AtomCGConverter::AtomCGConverter(vector<string> ignore_mol_types) {
  for (string &atom_mol_type : ignore_mol_types) {
    atom_mol_types_to_ignore_.insert(atom_mol_type);
  }
}

const std::string &AtomCGConverter::getCGMoleculeType(
    const string &atom_mol_type) const {
  assert(atom_and_cg_mol_types_.left.count(atom_mol_type) &&
         "atomistic molecule is not known");
  return atom_and_cg_mol_types_.left.at(atom_mol_type);
}

bool AtomCGConverter::AtomisticMoleculeTypeExist(
    const std::string &atom_mol_type) const {
  return atom_and_cg_mol_types_.left.count(atom_mol_type);
}

const std::string &AtomCGConverter::getAtomisticMoleculeType(
    const string &cg_mol_type) const {
  assert(atom_and_cg_mol_types_.right.count(cg_mol_type) &&
         "cg molecule is not known");
  return atom_and_cg_mol_types_.right.at(cg_mol_type);
}

Topology AtomCGConverter::Convert(const Topology &atom_top_in) {

  Topology cg_top_out;
  if (cg_top_out.BoundaryExist()) {
    cg_top_out.CopyBoundaryConditions(atom_top_in);
  }
  cg_top_out.setStep(atom_top_in.getStep());
  cg_top_out.setTime(atom_top_in.getTime());

  for (const Molecule &atom_mol : atom_top_in) {
    string atom_mol_type = atom_mol.getType();
    if (atom_mol_types_to_ignore_.count(atom_mol_type)) {
      continue;
    }

    if (AtomisticMoleculeTypeExist(atom_mol_type) == false) {
      cout << "--------------------------------------\n"
           << "WARNING: unknown molecule \"" << atom_mol.getType()
           << "\" with id " << atom_mol.getId() << " in topology" << endl
           << "molecule will not be mapped to CG representation\n"
           << "Check weather a mapping file for all molecule exists, was "
           << "specified in --cg separated by ; and the ident tag in xml-file "
           << "matches the molecule name\n"
           << "--------------------------------------\n";
      continue;
    }

    CoarseGrainMolecule_(atom_top_in, cg_top_out, atom_mol);
  }
  cg_top_out.RebuildExclusions();

  Update(atom_top_in, cg_top_out);

  return cg_top_out;
}

void AtomCGConverter::Update(const Topology &atom_top_in,
                             Topology &cg_top) const {

  assert(atom_top_in.getBoxType() == cg_top.getBoxType() &&
         "box types of topology in and out differ");

  cg_top.setStep(atom_top_in.getStep());
  cg_top.setTime(atom_top_in.getTime());
  cg_top.setBox(atom_top_in.getBox(), atom_top_in.getBoxType());

  // Cycle the cg molecules

  for (const pair<int, map<int, vector<pair<string, int>>>> &cg_mol_with_info :
       cgmolid_cgbeadid_atombeadnames_and_ids_) {

    int mol_id = cg_mol_with_info.first;
    string cg_mol_type = cg_top.getMolecule(mol_id).getType();
    string atom_mol_type = atom_top_in.getMolecule(mol_id).getType();

    // Call the appropriate molecule mapper
    mol_names_and_maps_.at(atom_mol_type)
        .at(cg_mol_type)
        .UpdateCGMolecule(atom_top_in, cg_top, cg_mol_with_info);
  }
}

Molecule &AtomCGConverter::CoarseGrainMolecule_(const Topology &atom_top,
                                                Topology &cg_top_out,
                                                const Molecule &atom_mol) {
  string atom_mol_type = atom_mol.getType();
  const int mol_id = atom_mol.getId();

  assert(cg_top_out.MoleculeExist(mol_id) == false &&
         "Cannot convert atomistic molecule to cg molecule because the cg "
         "molecule with the specified id already exists");
  string cg_mol_type = atom_and_cg_mol_types_.left.at(atom_mol_type);

  // returns ids of the cg molecule vector< atomic_bead_name, atomic_bead_id >
  //  map<int, vector<pair<string, int>>> cg_beads_to_atom_beads =
  //      CreateMolecule_(cg_mol_type, mol_id, cg_top_out, atom_top);

  Molecule &cg_mol = cg_top_out.CreateMolecule(mol_id, cg_mol_type);

  map<int, vector<pair<string, int>>> cg_beads_to_atom_beads =
      MapAtomicMoleculeToCGMolecule_(atom_top, cg_top_out, mol_id);
  cgmolid_cgbeadid_atombeadnames_and_ids_[mol_id] = cg_beads_to_atom_beads;
  return cg_mol;
}

void AtomCGConverter::LoadMoleculeStencil(const string &filename) {

  Property options;
  load_property_from_xml(options, filename);
  // Grab the type of the coarse grained molecule
  string cg_mol_type = options.get("cg_molecule.name").as<string>();
  // Grab the type of the atomistic molecule
  string atom_mol_type = options.get("cg_molecule.ident").as<string>();
  // Store the types in the bimap
  atom_and_cg_mol_types_.insert(
      boost::bimap<string, string>::value_type(atom_mol_type, cg_mol_type));

  // Create the stencil
  cg_mol_and_stencil_.insert(
      make_pair(cg_mol_type, CGMoleculeStencil(cg_mol_type, atom_mol_type)));

  vector<CGBeadStencil> beads_info =
      ParseBeads_(options.get("cg_molecule.topology"));
  // Update the stencil with the bead info
  cg_mol_and_stencil_.at(cg_mol_type).AddBeadStencils(beads_info);

  // Convert vector to map to be used with ParseBeadMaps
  unordered_map<string, CGBeadStencil> name_and_beads_info;
  for (CGBeadStencil &bead_info : beads_info) {
    name_and_beads_info[bead_info.cg_name_] = bead_info;
  }

  ParseBeadMaps_(options, name_and_beads_info);

  // Update the stencil the relevant interactions
  cg_mol_and_stencil_.at(cg_mol_type)
      .AddInteractionStencil(ParseBonded_(options.get("cg_molecule.topology")));

  // Create a mapper to map from the atom to the cg molecule
  mol_names_and_maps_[atom_mol_type].insert(std::make_pair(
      cg_mol_type, move(AtomToCGMoleculeMapper(atom_mol_type, cg_mol_type))));

  vector<string> order_of_beads =
      cg_mol_and_stencil_.at(cg_mol_type).getCGBeadNames();
  // Initialize the mapper
  mol_names_and_maps_.at(atom_mol_type)
      .at(cg_mol_type)
      .InitializeMoleculeMap(name_and_beads_info, order_of_beads);
}

std::unordered_map<int, string>
    AtomCGConverter::MapAtomicBeadIdsToAtomicBeadNames_(
        const string &cg_or_atom_mol_type, const vector<int> &bead_ids) const {

  assert((atom_and_cg_mol_types_.left.count(cg_or_atom_mol_type) ||
          atom_and_cg_mol_types_.right.count(cg_or_atom_mol_type)) &&
         "Cannot map atomic bead ids to atomic bead names because the molecule "
         "type is not recognized.");

  string cg_mol_type = cg_or_atom_mol_type;
  if (atom_and_cg_mol_types_.left.count(cg_or_atom_mol_type)) {
    cg_mol_type = atom_and_cg_mol_types_.left.at(cg_or_atom_mol_type);
  }

  return cg_mol_and_stencil_.at(cg_mol_type)
      .MapAtomicBeadIdsToAtomicBeadNames(bead_ids);
}

// Must use the cg bead name not the type to do this
vector<string> AtomCGConverter::getAtomicBeadNamesOfCGBead(
    const string &cg_mol_type, const string &cg_bead_name) const {
  assert(cg_mol_and_stencil_.count(cg_mol_type) &&
         "cg molecule type is not known to the atom-to-cg converter");
  return cg_mol_and_stencil_.at(cg_mol_type).getAtomicBeadNames(cg_bead_name);
}
/*****************************************************************************
 * Private Internal Methods
 *****************************************************************************/
vector<CGBeadStencil> AtomCGConverter::ParseBeads_(Property &options_in) {

  Property options = options_in.get("cg_beads");
  vector<Property *> beads = options.Select("cg_bead");

  unordered_set<string> bead_cg_names;
  vector<CGBeadStencil> beads_info;
  for (Property *p : beads) {
    CGBeadStencil bead_info;
    bead_info.cg_name_ = p->get("name").as<string>();
    bead_info.cg_bead_type_ = p->get("type").as<string>();
    bead_info.mapping_ = p->get("mapping").as<string>();
    bead_info.cg_symmetry_ = 1;
    if (p->exists("symmetry")) {
      bead_info.cg_symmetry_ = p->get("symmetry").as<int>();
    }

    // get the beads
    vector<string> subbeads;
    string bead_string(p->get("beads").value());
    Tokenizer tok_beads(bead_string, " \n\t");
    tok_beads.ToVector(subbeads);

    bead_info.atomic_subbeads_ = subbeads;

    if (bead_cg_names.count(bead_info.cg_name_)) {
      throw runtime_error(string("bead name ") + bead_info.cg_name_ +
                          " not unique in mapping");
    }
    bead_cg_names.insert(bead_info.cg_name_);
    beads_info.push_back(bead_info);
  }
  return beads_info;
}

void AtomCGConverter::ParseBeadMaps_(
    Property &options_in,
    unordered_map<string, CGBeadStencil> &bead_maps_info) {

  Property maps_prop = options_in.get("cg_molecule.maps");
  vector<Property *> all_maps = maps_prop.Select("map");

  for (Property *p : all_maps) {
    string map_type = p->get("name").as<string>();

    // get vector of weights
    vector<double> weights;
    Tokenizer tok_weights(p->get("weights").value(), " \n\t");
    tok_weights.ConvertToVector<double>(weights);

    // get vector of d values used for non-spherical beads
    vector<double> d;
    if (p->exists("d")) {
      Tokenizer tok_d(p->get("d").value(), " \n\t");
      tok_d.ConvertToVector(d);
    }

    for (pair<const string, CGBeadStencil> &pr : bead_maps_info) {
      if (pr.second.mapping_.compare(map_type) == 0) {
        pr.second.subbead_d_ = d;
        pr.second.subbead_weights_ = weights;
      }
    }
  }
}

bool AtomCGConverter::CheckThatBeadCountAndInteractionTypeAreConsistent_(
    const string &interaction_type, const size_t &bead_count) const noexcept {

  if (interaction_type == "bond") {
    if (bead_count != 2) {
      return false;
    }
  } else if (interaction_type == "angle") {
    if (bead_count != 3) {
      return false;
    }
  } else if (interaction_type == "dihedral") {
    if (bead_count != 4) {
      return false;
    }
  } else if (interaction_type == "improper") {
    return false;
  } else {
    return false;
  }
  return true;
}

vector<CGInteractionStencil> AtomCGConverter::ParseBonded_(
    Property &options_in) {

  vector<CGInteractionStencil> all_interactions_info;
  if (options_in.exists("cg_bonded")) {
    Property options = options_in.get("cg_bonded");
    std::vector<Property *> bonded = options.Select("*");

    vector<Property *>::iterator bond;
    unordered_set<string> interaction_groups;
    // Parse through the different bond sections
    for (bond = bonded.begin(); bond != bonded.end(); ++bond) {
      string interaction_group = (*bond)->get("name").as<string>();
      if (interaction_groups.count(interaction_group)) {
        throw runtime_error(
            string("double occurence of interactions with name ") +
            interaction_group);
      }
      interaction_groups.insert(interaction_group);

      Tokenizer section((*bond)->get("beads").value(), "\n");
      for (Tokenizer::iterator line = section.begin(); line != section.end();
           ++line) {

        CGInteractionStencil interaction_info;

        Tokenizer tok_atoms(*line, " \t");
        int atom_count = 0;
        for (Tokenizer::iterator atom = tok_atoms.begin();
             atom != tok_atoms.end(); ++atom) {
          interaction_info.bead_names_.push_back(*atom);
          ++atom_count;
        }
        if (interaction_info.bead_names_.size() == 0) continue;
        interaction_info.group_ = interaction_group;
        if (atom_count == 2) {
          interaction_info.type_ = "bond";
        } else if (atom_count == 3) {
          interaction_info.type_ = "angle";
        } else if (atom_count == 4) {
          interaction_info.type_ = "dihedral";
        }

        assert(CheckThatBeadCountAndInteractionTypeAreConsistent_(
                   (*bond)->name(), interaction_info.bead_names_.size()) &&
               "Consistency check failed,");

        all_interactions_info.push_back(interaction_info);
      }
    }
  }
  return all_interactions_info;
}

// This function not only creates the beads associated with a cg molecule
// it returns a map expressing how the cg beads of that molecule are associated
// with the atomic beads of a specific atomic molecule
//
// map int - cg bead id
//     vector< atom name, atom id>
//
// The only way this is possible is to assume that the cg molecules and the
// atomic molecules share the same ids, also assumes that the atomic bead ids
// when sorted should line up with the cg beads when they are read sequentially
// from an .xml file
//
map<int, vector<pair<string, int>>> AtomCGConverter::MapAtomicBeadsToCGBeads_(
    const Topology &atom_top, Topology &cg_top_out, Molecule &cg_mol,
    const CGMoleculeStencil &stencil) const {

  map<string, int> cg_bead_name_and_id;
  for (const CGBeadStencil &bead_info : stencil.getBeadStencil()) {

    string bead_type = bead_info.cg_bead_type_;
    Bead &bead = cg_top_out.CreateBead(
        bead_info.cg_symmetry_, bead_type, cg_top_out.BeadCount(),
        cg_mol.getId(), topology_constants::unassigned_residue_id,
        topology_constants::unassigned_residue_type,
        topology_constants::unassigned_element, 0.0, 0.0);

    cg_bead_name_and_id[bead_info.cg_name_] = bead.getId();
    cg_mol.AddBead(bead);
  }

  // cg_bead_id, vector< atom_name, atom_bead_id >
  const Molecule &atom_mol = atom_top.getMolecule(cg_mol.getId());
  vector<int> atom_bead_ids = atom_mol.getBeadIds();
  sort(atom_bead_ids.begin(), atom_bead_ids.end());
  unordered_map<int, string> atom_ids_and_names =
      MapAtomicBeadIdsToAtomicBeadNames_(atom_mol.getType(), atom_bead_ids);
  vector<int> atom_ids;
  for (auto id_and_name : atom_ids_and_names) {
    atom_ids.push_back(id_and_name.first);
  }
  sort(atom_ids.begin(), atom_ids.end());

  map<int, vector<pair<string, int>>> cg_beads_and_atoms;

  int min_index = 0;
  int max_index = 0;
  for (const CGBeadStencil &bead_info : stencil.getBeadStencil()) {

    vector<pair<string, int>> atom_bead_names_ids;
    max_index += bead_info.atomic_subbeads_.size();
    for (int index = min_index; index < max_index; ++index) {

      assert(
          static_cast<size_t>(index) < atom_ids.size() &&
          "ERROR there is a problem in the conversion of the atomic molecule "
          "to the coarse grained molecule, there is a conflict with the number "
          "of beads needed by the coarse grained description and the number of "
          "atoms that are in the atomic molecule");

      int atom_id = atom_ids.at(index);
      string atom_name = atom_ids_and_names.at(atom_id);
      atom_bead_names_ids.push_back(pair<string, int>(atom_name, atom_id));
    }
    min_index = max_index;

    string cg_bead_name = bead_info.cg_name_;
    int cg_bead_id = cg_bead_name_and_id.at(cg_bead_name);
    cg_beads_and_atoms[cg_bead_id] = atom_bead_names_ids;
  }

  return cg_beads_and_atoms;
}

void AtomCGConverter::CreateInteractions_(
    Topology &cg_top_out, Molecule &cg_mol,
    const CGMoleculeStencil &stencil) const {

  // Convert to a map of the atomic bead names and their ids
  vector<int> cg_bead_ids = cg_mol.getBeadIds();
  unordered_map<int, string> bead_ids_and_names =
      stencil.MapCGBeadIdsToCGBeadNames(cg_bead_ids);
  unordered_map<string, int> bead_name_to_id;
  for (pair<const int, string> &value : bead_ids_and_names) {
    bead_name_to_id[value.second] = value.first;
  }

  for (const CGInteractionStencil &interaction_info :
       stencil.getInteractionStencil()) {
    // Convert atoms to vector of ints using the map
    size_t interaction_id = cg_top_out.InteractionCount();
    vector<int> atoms;
    for (const string &atom_name : interaction_info.bead_names_) {
      atoms.push_back(bead_name_to_id.at(atom_name));
    }

    if (!atoms.empty()) {
      Interaction *ic;
      if (interaction_info.type_ == "bond") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::bond, interaction_info.group_, interaction_id,
            cg_mol.getId(), atoms);
      } else if (interaction_info.type_ == "angle") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::angle, interaction_info.group_, interaction_id,
            cg_mol.getId(), atoms);
      } else if (interaction_info.type_ == "dihedral") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::dihedral, interaction_info.group_, interaction_id,
            cg_mol.getId(), atoms);
      } else {
        throw runtime_error("unknown bonded type in map: " +
                            interaction_info.type_);
      }
      cg_mol.AddInteraction(ic);
    }
  }
}

// Should return the id of the cg bead
// followed by the vector< atom name, atom id>

map<int, vector<pair<string, int>>>
    AtomCGConverter::MapAtomicMoleculeToCGMolecule_(const Topology &atom_top,
                                                    Topology &cg_top_out,
                                                    const int &mol_id) {
  // map<int, vector<pair<string, int>>> AtomCGConverter::CreateMolecule_(
  //    string cg_mol_type, int mol_id, Topology &cg_top_out,
  //    Topology &atom_top) {

  Molecule &cg_mol = cg_top_out.getMolecule(mol_id);
  string cg_mol_type = cg_mol.getType();

  CGMoleculeStencil &cg_mol_stencil = cg_mol_and_stencil_.at(cg_mol_type);

  map<int, vector<pair<std::string, int>>> bead_name_to_id =
      MapAtomicBeadsToCGBeads_(atom_top, cg_top_out, cg_mol, cg_mol_stencil);

  CreateInteractions_(cg_top_out, cg_mol, cg_mol_stencil);

  return bead_name_to_id;
}

}  // namespace csg
}  // namespace votca
