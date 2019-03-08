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
#include <votca/csg/atomcgconverter.h>
#include <votca/csg/bead.h>
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

/*****************************************************************************
 * Public Facing Methods
 *****************************************************************************/
AtomCGConverter::AtomCGConverter(vector<string> ignore_molecule_types){
  for( string & atomistic_mol_type : ignore_molecule_types){
    atomistic_molecule_types_to_ignore_.insert(atomistic_mol_type);
  }
}

const std::string &AtomCGConverter::getCGMoleculeType(
    string atomistic_molecule_type_) const {
  assert(atomic_and_cg_molecule_types.left.count(atomistic_molecule_type_) &&
         "atomistic molecule is not known");
  return atomic_and_cg_molecule_types_.left.at(atomistic_molecule_type_);
}

bool AtomCGConverter::AtomisticMoleculeTypeExist(std::string atomistic_mol_type){
  return atomic_and_cg_molecule_types_.left.count(atomistic_mol_type);
}

const std::string &AtomCGConverter::getAtomisticMoleculeType(
    string cg_molecule_type) const {
  assert(atomic_and_cg_molecule_types.right.count(cg_molecule_type) &&
         "cg molecule is not known");
  return atomic_and_cg_molecule_types_.right.at(cg_molecule_type);
}

void AtomCGConverter::Convert(CSG_Topology & atomic_top_in, CSG_Topology & cg_top_out){

    assert(atomistic_top_in.getBoxType() == cg_top_out.getBoxType() &&            
           "box types of topology in and out differ");

    // Grab all the molecules                                                     
    const unordered_map<int, Molecule> &atomistic_mols =                          
        atomic_top_in.Molecules();                                             
                                                                                 
                                                                                 
    // Cycle through the atomistic molecules                                      
   for (const pair<int, Molecule> &id_and_molecule : atomistic_mols) {           
     const Molecule *atomistic_mol = &(id_and_molecule.second);                  
                                                                                 
     string atomistic_mol_type = atomistic_mol->getType();                       
     if (atomistic_molecular_types_to_ignore_.count(atomistic_mol_type)){ 
       continue;
     }
                                                                                 
     if ( AtomisticMoleculeTypeExist(atomistic_mol_type) == false ){   
       cout << "--------------------------------------\n"                        
         << "WARNING: unknown molecule \"" << atomistic_mol->getType()           
         << "\" with id " << atomistic_mol->getId() << " in topology"            
         << endl                                                                 
         << "molecule will not be mapped to CG representation\n"                 
         << "Check weather a mapping file for all molecule exists, was "         
        << "specified in --cg separated by ; and the ident tag in xml-file "    
         << "matches the molecule name\n"                                        
         << "--------------------------------------\n";                          
       continue;                                                                 
                                                                                 
     }                                                                           
                                                                                 
    ConvertAtomisticMoleculeToCGAndAddToCGTopology_(*atomistic_mol,cg_top_out);       
                                                                                  
    }                                                                             
    cg_top_out.RebuildExclusions();   

    Map(atomic_top_in, cg_top_out);
}

void AtomCGConverter::Map(CSG_Topology & atomic_top, CSG_Topology & cg_top){

  assert(atomistic_top_in.getBoxType() == cg_top_out.getBoxType() &&            
      "box types of topology in and out differ");
  cg_top.setStep(atomic_top.getStep());                                  
  cg_top.setTime(atomic_top.getTime());                                  
  cg_top.setBox(atomic_top.getBox());                                    

  // Get the cg molecules                                                       
  vector<int> cg_molecule_ids = cg_top.getMoleculeIds();                      
  for (int &cg_molecule_id : cg_molecule_ids) {                                 
    Molecule *cg_molecule = cg_top.getMolecule(cg_molecule_id);               
    string cg_mol_type = cg_molecule->getType();                           
    int molecule_id = cg_molecule->getId();                                     
    vector<int> cg_bead_ids = cg_molecule->getBeadIds();                        
    unordered_map<int, string> bead_ids_and_names =                             
      MapAtomicBeadIdsToAtomicBeadNames_(cg_mol_type, cg_bead_ids);

    // All the atomic beads in the molecule
    map<string, Bead *> name_and_atomic_bead;                                   
    for (pair<int, string> id_and_bead_name : bead_ids_and_names) {             
      name_and_atomic_bead[id_and_bead_name.second] =                           
        atomic_top.getBead(id_and_bead_name.first);                      
    }                                                                           

    string atomic_mol_type = atomic_top.getMolecule(molecule_id)->getType();

    // Calling the Molecular mapper
    molecule_names_and_maps_[atomic_mol_type][cg_mol_type].Apply(
        cg_top.getBoundaryCondition(),name_and_atomic_bead); 
  }             
}

void AtomCGConverter::ConvertAtomisticMoleculeToCGAndAddToCGTopology_(                 
     const Molecule & atomistic_mol,CSG_Topology & cg_top_out){                        
                                                                                 
   string atomistic_molecule_type_ = atomistic_mol.getType();                    
   int molecule_id = atomistic_mol.getId();                                      
                                                                                 
   assert(cg_top_out.Molecule.Exist(molecule_id) == false &&                     
       "Cannot convert atomistic molecule to cg molecule because the cg "        
       "molecule with the specified id already exists");                         
   string cg_molecule_type = atomic_and_cg_molecule_types.left.at(atomistic_molecule_type_);     
  
   CreateMolecule_(cg_molecule_type, molecule_id, cg_top_out);                   
   /**************************************************************************
    * It is here that we need to record which atom ids belong to a given cg_id
    **************************************************************************/
 }    

void AtomCGConverter::LoadMoleculeStencil(string filename) {

  Property *options;
  load_property_from_xml(options, filename);
  // Grab the type of the coarse grained molecule
  string cg_molecule_type = options_.get("cg_molecule.name").as<string>();
  // Grab the type of the atomistic molecule
  string atomistic_molecule_type_ =
      options_.get("cg_molecule.ident").as<string>();

  // Store the types in the bimap
  atomic_and_cg_molecule_types.insert(atomic_and_cg_molecule_types::value_type(
      cg_molecule_type, atomistic_molecule_type_));

  // Create the stencil
  cg_molecule_and_stencil[cg_molecule_type] =
      CGStencil(cg_molecule_type, atomistic_molecule_type_);

  vector<CGBeadInfo> beads_info = ParseBeads_(options.get("cg_molecule.topology");
  // Update the stencil with the bead info
  cg_molecule_and_stencil[cg_molecule_type].AddBeadInfo(beads_info));

  // Convert vector to map to be used with ParseMaps
  unordered_map<string,CGBeadInfo> map_and_beads_info;
  for( CGBeadInfo & bead_info : beads_info ){
    map_and_beads_info[bead_info.mapping_] = bead_info;
  }

  ParseMaps_(options,man_and_beads_info);

  // Update the stencil the relevant interactions
  cg_molecule_and_stencil[cg_molecule_type].AddInteractionInfo(
      ParseBonded_(options.get("cg_molecule.topology")));

  // Create a mapper to map from the atom to the cg molecule
  molecule_names_and_maps_[atomistic_molecule_type][cg_molecule_type] =
    AtomisticToCGMoleculeMapper(atomistic_molecule_type,cg_molecule_type);
  
  // Initialize the mapper
  molecule_names_and_maps_[atomistic_molecule_type][cg_molecule_type].Initialize( map_and_beads_info);
  
}

std::unordered_map<int, string>
    AtomCGConverter::MapAtomicBeadIdsToAtomicBeadNames(
        string cg_or_atomic_molecule_type, vector<int> bead_ids) {

  assert(
      (atomic_and_cg_molecule_types_.left.count(cg_or_atomic_molecule_type) ||
       atomic_and_cg_molecule_types_.right.count(cg_or_atomic_molecule_type)) &&
      "Cannot map atomic bead ids to atomic bead names because the molecule "
      "type is not recognized.");

  string cg_molecule_type = cg_or_atomic_molecule_type;
  if (atomic_and_cg_molecule_types_.left.count(cg_or_atomic_molecule_type)) {
    cg_molecule_type =
        atomic_and_cg_molecule_types_.at(cg_or_atomic_molecule_type);
  }

  file_names_.insert(filename);

  return cg_molecule_and_stencil[cg_molecule_type]
      .MapAtomicBeadIdsToAtomicBeadNames(bead_ids);
}

// Must use the cg bead name not the type to do this
vector<string> AtomCGConverter::getAtomicBeadNamesOfCGBead(
    string cg_molecule_type, string cg_bead_name) {
  assert(cg_molecule_and_stencil.count(cg_molecule_type) &&
         "cg molecule type is not known to the atom-to-cg converter");
  cg_molecule_and_stencil[cg_molecule_type].getAtomicBeadNames(cg_bead_name);
}
/*****************************************************************************
 * Private Internal Methods
 *****************************************************************************/
std::unordered_map<string,CGBeadInfo> AtomCGConverter::ParseBeads_(Property &options_in) {

  Property options = options.in.get("cg_beads");
  list<Property *> beads = options.Select("cg_bead");

  unordered_set<string> bead_cg_names;
  vector<CGBeadInfo> beads_info;
  for (list<Property *>::iterator bead_iter = beads.begin();
       bead_iter != beads.end(); ++bead_iter) {
    Property *p = *bead_iter;
    CGBeadInfo bead_info;
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

    if (bead_cg_names.count(bead_info->cg_name_)) {
      throw runtime_error(string("bead name ") + bead_info.cg_name_ +
                          " not unique in mapping");
    }
    bead_cg_names.insert(bead_info->cg_name_);
    beads_info.push_back(bead_info);
  }
  return beads_info;
}

 void AtomToCGConverter::ParseMaps_(Property &options_in,                              
                          unordered_map<string, CGBeadInfo> &bead_maps_info) {   
                                                                                  
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


void AtomCGConverter::CheckThatBeadCountAndInteractionTypeAreConsistent_(
    string interaction_type, size_t bead_count) const {

  if (interaction_type == "bond") {
    if (bead_count != 2) {
      throw runtime_error(
          "Error interaction type is specified as bond but"
          " there are not exactly two atoms associated with the bond");
    }
    return;
  } else if (interaction_type == "angle") {
    if (bead_count != 3) {
      throw runtime_error(
          "Error interaction type is specified as angle "
          "but there are not exactly three atoms associated with the bond");
    }
    return;
  } else if (interaction_type == "dihedral") {
    if (bead_count != 4) {
      throw runtime_error(
          "Error interaction type is specified as dihedral"
          " but there are not exactly four atoms associated with the bond");
    }
    return;
  } else if (interaction_type == "improper") {
    throw runtime_error(
        "Error currently the improper interaction type is not "
        "supported.");
  } else {
    throw runtime_error(string("Error the specified interaction type ") +
                        interaction_type + string(" is not known."));
  }
}

vector<CGInteractionInfo> AtomCGConverter::ParseBonded_(Property &options_in) {

  vector<CGInteractionInfo> all_interactions_info;
  if (options_in.exists("cg_bonded")) {
    Property options = options_in.get("cg_bonded");
    std::list<Property *> bonded = options.Select("*");

    list<Property *>::iterator bond;
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

        CGInteractionInfo interaction_info;
        Tokenizer tok_atoms(*line, " \t");
        for (Tokenizer::iterator atom = tok_atoms.begin();
             atom != tok_atoms.end(); ++atom) {
          interaction_info.atoms_.push_back(*atom);
        }

        CheckThatBeadCountAndInteractionTypeAreConsistent_(
            (*bond)->name(), interaction_info.atoms_.size());

        all_interactions_info.push_back(interaction_info);
      }
    }
  }
  return all_interactions_info;
}
/*
void AtomCGConverter::ParseMapping(Property &options) {
  list<Property *> maps = options.Select("map");

  for (list<Property *>::iterator iter = maps.begin(); iter != maps.end();
       ++iter) {
    maps_[(*iter)->get("name").as<string>()] = *iter;
  }
}*/

unordered_map<string, int> AtomCGConverter::CreateBeads_(
    Molecule *cg_mol, CGStencil stencil, CSG_Topology &cg_top_out) {
  unordered_map<std::string, int> bead_name_to_id;
  for (CGBeadInfo &bead_info : stencil.getBeadInfo()) {
    Bead *bead;

    string bead_type = bead_info.type_;
    bead = cg_top_out.CreateBead(
        bead_info.symmetry_, bead_type, cg_top_out.BeadCount(), cg_mol->getId(),
        topology_constants::unassigned_residue_id,
        topology_constants::unassigned_residue_type,
        topology_constants::unassigned_element, 0.0, 0.0);
    cg_mol->AddBead(bead);
    bead_name_to_id[bead_info.cg_name_] = bead->getId();
  }

  return bead_name_to_id;
}

void AtomCGConverter::CreateInteractions_(
    Molecule *cg_mol, CGStencil stencil, CSG_Topology &cg_top_out,
    unordered_map<std::string, int> bead_name_to_id) {
  for (CGInteractionInfo &interaction_info : stencil.getInteractionInfo()) {
    // Convert atoms to vector of ints using the map
    size_t interaction_id = cg_top_out.InteractionCount();
    vector<int> atoms;
    for (string atom_name : interaction_info.bead_names_) {
      atoms.push_back(bead_name_to_id[atom_name]);
    }

    if (!atoms.empty()) {
      Interaction *ic;
      if (interaction_info.type_ == "bond") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::bond, interaction_info.group_, interaction_id,
            cg_mol->getId(), atoms);
      } else if (interaction_info.type_ == "angle") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::angle, interaction_info.group_, interaction_id,
            cg_mol->getId(), atoms);
      } else if (interaction_info.type_ == "dihedral") {
        ic = cg_top_out.CreateInteraction(
            InteractionType::dihedral, interaction_info.group_, interaction_id,
            cg_mol->getId(), atoms);
      } else {
        throw runtime_error("unknown bonded type in map: " +
                            interaction_info.type_);
      }
      cg_mol->AddInteraction(ic);
    }
  }
}

Molecule *AtomCGConverter::CreateMolecule_(string cg_molecule_type,
                                           int molecule_id,
                                           CSG_Topology &cg_top_out) {

  Molecule *cg_mol = cg_top_out.CreateMolecule(molecule_id, cg_molecule_type);
  CGStencil &cg_mol_stencil = cg_molecule_and_stencil.at(cg_molecule_type);

  unordered_map<std::string, int> bead_name_to_id =
      CreateBeads_(cg_mol, cg_mol_stencil, cg_top_out);
  CreateInteractions_(cg_mol, cg_mol_stencil, cg_top_out, bead_name_to_id);

  return cg_mol;
}

Map *AtomCGConverter::CreateMap(const BoundaryCondition *boundaries,
                                const Molecule &mol_in, Molecule &mol_out) {
  if ((unsigned int)mol_out.BeadCount() != beads_.size()) {
    throw runtime_error(
        "number of beads for cg molecule and mapping definition do "
        "not match, check your molecule naming.");
  }

  Map *map = new Map(mol_in, mol_out);
  for (vector<CGBeadInfo *>::iterator beaddef = beads_.begin();
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

    map->CreateBeadMap((*beaddef)->symmetry_, boundaries, &mol_in,
                       mol_out.getBead(bead_id), ((*beaddef)->options_), mdef);
    /*BeadMap *bmap;
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
    bmap->Initialize(boundaries, &mol_in, mol_out.getBead(bead_id),
                     ((*beaddef)->options_), mdef);
    map->AddBeadMap(bmap);*/
  }
  return map;
}
/*
AtomCGConverter::CGBeadInfo *AtomCGConverter::getBeadByCGName(
    const string &cg_name) {
  map<string, CGBeadInfo *>::iterator cg_name_and_beaddef =
      beads_by_cg_name_.find(cg_name);
  if (cg_name_and_beaddef == beads_by_cg_name_.end()) {
    std::cout << "cannot find: <" << cg_name << "> in " << cg_molecule_type_
              << "\n";
    return nullptr;
  }
  return (*cg_name_and_beaddef).second;
}*/
/*
Property *AtomCGConverter::getMapByName(const string &mapping_name) {
  map<string, Property *>::iterator iter = maps_.find(mapping_name);
  if (iter == maps_.end()) {
    std::cout << "cannot find map " << mapping_name << "\n";
    return NULL;
  }
  return (*iter).second;
}*/

}  // namespace csg
}  // namespace votca
