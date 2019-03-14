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

#include "../../include/votca/csg/boundarycondition.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <numeric>
#include <memory>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/beadmap.h>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/vec.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

vector<string> BeadMap::getAtomicBeadNames(){
  vector<string> bead_names;
  for( pair<string,element_t> pr : matrix_){
    bead_names.push_back(pr.first);
  }
  return bead_names;
}

void Map_Sphere::Initialize(vector<string> subbeads, vector<double> weights,
                            vector<double> ds) {

  assert(subbeads.size() == weights.size() &&
         "subbeads and weights are not matched in bead map sphere.");

  vector<double> fweights;
  // normalize the weights
  double norm = 1. / std::accumulate(weights.begin(), weights.end(), 0.);

  transform(weights.begin(), weights.end(), weights.begin(),
            bind2nd(multiplies<double>(), norm));
  // get the d vector if exists or initialize same as weights
  if (ds.size() > 0) {
    // normalize d coefficients
    norm = 1. / std::accumulate(ds.begin(), ds.end(), 0.);
    transform(ds.begin(), ds.end(), ds.begin(),
              bind2nd(multiplies<double>(), norm));
  } else {
    // initialize force-weights with weights
    ds.resize(weights.size());
    copy(weights.begin(), weights.end(), ds.begin());
  }

  assert(subbeads.size() == ds.size() &&
         "subbeads and ds are not matched in bead map sphere.");

  fweights.resize(weights.size());
  // calculate force weights by d_i/w_i
  for (size_t i = 0; i < weights.size(); ++i) {
    if (weights[i] == 0 && ds[i] != 0) {
      throw runtime_error(
          "A d coefficient is nonzero while weights is zero in mapping ");
    }
    if (weights[i] != 0) {
      fweights[i] = ds[i] / weights[i];
    } else {
      fweights[i] = 0;
    }
  }

  int index = 0;
  for (string &name : subbeads) {
    AddElem(name, weights.at(index), fweights.at(index));
    ++index;
  }
  
}

void Map_Sphere::Apply(const BoundaryCondition *boundaries,
                       map<string, Bead *> atomic_beads, Bead *cg_bead) {

  cout << "Applying sphere" << endl;
  assert(cg_bead->getSymmetry()==1 && "Applying spherical bead map on a non spherical corase grained bead");

  assert(matrix_.size() == atomic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified.");
    Bead *atom = atomic_beads.begin()->second;
  assert(atom->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  
  vec reference_position = atom->getPos();
  string bead_type = atom->getType();
  int bead_id = atom->getId();

  bool bVel, bF;
  bVel = bF = false;
  double sum_of_atomistic_mass = 0;
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_velocity(0., 0., 0.);

  for (pair<const string, element_t> &name_and_element : matrix_) {
    cout << "Bead name mapping " << name_and_element.first << " exists " << atomic_beads.count(name_and_element.first)<< endl;
    cout << " weight " << name_and_element.second.weight_ << endl;
    atom = atomic_beads[name_and_element.first];
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    vec shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());

    cout << "Boundary type " << boundaries->getBoxType() << endl;
    if(boundaries->getBoxType()!=BoundaryCondition::eBoxtype::typeOpen){
      double max_dist = 0.5 * boundaries->getShortestBoxDimension();
      if (abs(shortest_distance_beween_beads) > max_dist) {
        cout << reference_position << " " << atom->getPos() << endl;
        throw std::runtime_error(
            "coarse-grained atom is bigger than half the box \n (atoms " +
            bead_type + " (id " + boost::lexical_cast<string>(bead_id + 1) + ")" +
            ", " + atom->getType() + " (id " +
            boost::lexical_cast<string>(atom->getId() + 1) + ")" +
            +" , molecule " +
            boost::lexical_cast<string>(atom->getMoleculeId() + 1) + ")");
      }
    }
    weighted_sum_of_atomistic_pos +=
        name_and_element.second.weight_ *
        (shortest_distance_beween_beads + reference_position);
    if (atom->HasVel()) {
      weighted_sum_of_atomistic_velocity +=
          name_and_element.second.weight_ * atom->getVel();
      bVel = true;
    }
    if (atom->HasF()) {
      weighted_sum_of_atomistic_forces +=
          name_and_element.second.force_weight_ * atom->getF();
      bF = true;
    }
  }
  cg_bead->setMass(sum_of_atomistic_mass);
  cg_bead->setPos(weighted_sum_of_atomistic_pos);
  if (bVel) cg_bead->setVel(weighted_sum_of_atomistic_velocity);
  if (bF) cg_bead->setF(weighted_sum_of_atomistic_forces);
}

/// Warning the atomistic beads must be a map they cannot be an unordered_map
void Map_Ellipsoid::Apply(const BoundaryCondition *boundaries,
                          std::map<string, Bead *> atomistic_beads,
                          Bead *cg_bead) {


  assert(cg_bead->getSymmetry()==3 && "Applying ellipsoidal bead map on a non ellipsoidal corase grained bead");
  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads in "
         "Map_Ellipsoid");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified in "
         "Map_Ellipsoid.");

  std::map<string, Bead *>::iterator name_and_bead_iter;
  name_and_bead_iter = atomistic_beads.begin();
  assert(name_and_bead_iter->second->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  vec reference_position = name_and_bead_iter->second->getPos();
  ++name_and_bead_iter;

  bool bVel, bF;
  bVel = bF = false;
  double sum_of_atomistic_mass = 0;
  vec sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_vel(0., 0., 0.);
  for (pair<string, element_t> name_and_element : matrix_) {
    const Bead *atom = atomistic_beads[name_and_element.first];
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    vec shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());
    if(boundaries->getBoxType()!=BoundaryCondition::eBoxtype::typeOpen){
      double max_dist = 0.5 * boundaries->getShortestBoxDimension();
      if (abs(shortest_distance_beween_beads) > max_dist) {
        throw std::runtime_error(
            "coarse-grained atom is bigger than half the box");
      }
    }
    assert(atom->HasPos() && "Cannot map to coarse grained bead as atomic beads are missing positions");
    assert(name_and_element.second.weight_ && "Cannot map to coarse grained beads with weights of 0.0");
    weighted_sum_of_atomistic_pos +=
      name_and_element.second.weight_ *
      (shortest_distance_beween_beads + reference_position);
    if (atom->HasVel() == true) {
      weighted_sum_of_atomistic_vel +=
          name_and_element.second.weight_ * atom->getVel();
      bVel = true;
    }
    if (atom->HasF()) {
      weighted_sum_of_atomistic_forces +=
          name_and_element.second.force_weight_ * atom->getF();
      bF = true;
    }
  }

  cg_bead->setMass(sum_of_atomistic_mass);
  cg_bead->setPos(weighted_sum_of_atomistic_pos);
  if (bVel) cg_bead->setVel(weighted_sum_of_atomistic_vel);
  if (bF) cg_bead->setF(weighted_sum_of_atomistic_forces);

  // calculate the tensor of gyration
  matrix tensor_of_gyration(0.);
  double number_of_atom_beads = static_cast<double>(matrix_.size());
  vec average_pos = weighted_sum_of_atomistic_pos / number_of_atom_beads;
  for (pair<string, element_t> name_and_element : matrix_) {
    if (name_and_element.second.weight_ == 0) continue;
    const Bead *atom = atomistic_beads[name_and_element.first];
    vec pos = atom->getPos() - average_pos;

    // Normalize the tensor with 1/number_of_atoms_per_bead
    tensor_of_gyration[0][0] +=
        pos.getX() * pos.getX() / number_of_atom_beads;
    tensor_of_gyration[0][1] +=
        pos.getX() * pos.getY() / number_of_atom_beads;
    tensor_of_gyration[0][2] +=
        pos.getX() * pos.getZ() / number_of_atom_beads;
    tensor_of_gyration[1][1] +=
        pos.getY() * pos.getY() / number_of_atom_beads;
    tensor_of_gyration[1][2] +=
        pos.getY() * pos.getZ() / number_of_atom_beads;
    tensor_of_gyration[2][2] +=
        pos.getZ() * pos.getZ() / number_of_atom_beads;
  }
  tensor_of_gyration[1][0] = tensor_of_gyration[0][1];
  tensor_of_gyration[2][0] = tensor_of_gyration[0][2];
  tensor_of_gyration[2][1] = tensor_of_gyration[1][2];

  // calculate the eigenvectors
  matrix::eigensystem_t es;
  tensor_of_gyration.SolveEigensystem(es);

  // I do not like this, this is arbitrarily using the first three beads
  // to determine the orientation parameters
  vec reference_position2 = name_and_bead_iter->second->getPos();
  ++name_and_bead_iter;
  vec reference_position3 = name_and_bead_iter->second->getPos();

  vec u = es.eigenvecs[0];
  vec v = reference_position2 - reference_position;
  v.normalize();

  cg_bead->setV(v);

  // vec w = matrix_[2].bead_in_->getPos() - matrix_[0].bead_in_->getPos();
  vec w = reference_position3 - reference_position; w.normalize();

  if ((v ^ w) * u < 0) { u = vec(0., 0., 0.) - u;}
  cg_bead->setU(u);

  // write out w
  w = u ^ v;
  w.normalize();
  cg_bead->setW(w);
}

AtomToCGMoleculeMapper::~AtomToCGMoleculeMapper() {
  cg_bead_name_and_maps_.clear();
}

// cg_bead_order contains the order of beads names
void AtomToCGMoleculeMapper::Initialize(
    unordered_map<string, CGBeadStencil> bead_maps_info,vector<string> cg_bead_order) {

  for( string & cg_bead_name : cg_bead_order ){
  //for (pair<const string, CGBeadStencil> & bead_info : bead_maps_info) {
    CGBeadStencil & bead_info = bead_maps_info[cg_bead_name];
    cout << "Symmetry " << static_cast<int>(bead_info.cg_symmetry_) << endl;
    int symmetry = static_cast<int>(bead_info.cg_symmetry_);
    //string cg_bead_name = bead_info.cg_name_;
    switch (symmetry) {
      case 1:
        cg_bead_name_and_maps_.insert(make_pair(cg_bead_name, 
            unique_ptr<Map_Sphere>(new Map_Sphere())));
        break;
      case 3:
        cg_bead_name_and_maps_.insert(make_pair(cg_bead_name,
           unique_ptr<Map_Ellipsoid>(new Map_Ellipsoid())));
        break;
      default:
        throw runtime_error("unknown symmetry in bead definition!");
    }

    //cout << "cycle " << bead_info.first << endl;
    for (const string & subbead : bead_info.atomic_subbeads_){
      cout << subbead << " "; 
    }
    cout << endl;
    bead_type_and_names_[bead_info.cg_bead_type_].push_back(cg_bead_name);
    cg_bead_name_and_maps_.at(cg_bead_name)->Initialize(
        bead_info.atomic_subbeads_, bead_info.subbead_weights_);
  }
}

/// Copy Constructor                                                          
AtomToCGMoleculeMapper::AtomToCGMoleculeMapper(const AtomToCGMoleculeMapper & other){ 
  atom_molecule_type_ = other.atom_molecule_type_;                            
  cg_molecule_type_ = other.cg_molecule_type_;                                 

  for( const std::pair<const std::string,std::unique_ptr<BeadMap>> & pr  : other.cg_bead_name_and_maps_){    
    cg_bead_name_and_maps_.at(pr.first) = pr.second->Clone();               
  }                                                                         
}                                                                          

// Move assignment                                                            
AtomToCGMoleculeMapper & AtomToCGMoleculeMapper::operator=(AtomToCGMoleculeMapper&& other){           
  if(this!=&other){                                                           
    cg_bead_name_and_maps_.clear();                                           
    for( std::pair<const std::string,std::unique_ptr<BeadMap>> & pr  : other.cg_bead_name_and_maps_){
      cg_bead_name_and_maps_.at(pr.first) = std::move(pr.second);             
    }                                                                         
    atom_molecule_type_ = other.atom_molecule_type_;                          
    cg_molecule_type_ = other.cg_molecule_type_;                              
  }                                                                           
  return *this;                                                               
}                                                                             

+                                                                                 
// Copy assignment                                                            
AtomToCGMoleculeMapper & AtomToCGMoleculeMapper::operator=(const AtomToCGMoleculeMapper other){       
  if(this!=&other){                                                           
    cg_bead_name_and_maps_.clear();                                           
    for(const std::pair<const std::string,std::unique_ptr<BeadMap>> & pr  : other.cg_bead_name_and_maps_){
      cg_bead_name_and_maps_.at(pr.first) = pr.second->Clone();               
    }                                                                         
    atom_molecule_type_ = other.atom_molecule_type_;                          
    cg_molecule_type_ = other.cg_molecule_type_;                              
  }                                                                           
  return *this;                                                               
}           


void AtomToCGMoleculeMapper::Apply(
    CSG_Topology &atom_top, 
    CSG_Topology& cg_top, 
    pair<int,map<int,vector<pair<string,int>>>> cgmolid_cgbeadid_atomicbeadids ){

  // First int cg_molecule id
  // map
  //   first int - the cg_bead_id 
  //   vector pair
  //       string - atomic name 
  //       int - global atomic bead id


  // Grab the correct cg molecule
  int cg_mol_id = cgmolid_cgbeadid_atomicbeadids.first;
  Molecule * cg_mol = cg_top.getMolecule(cg_mol_id);
 
  // Ensure that the type molecule type is correct
  assert(cg_mol->getType() == cg_molecule_type_ && "Cannot convert to the molecule  is not of the correct type");
  
  // Cycle throught the cg beads in a cg molecule
  vector<int> bead_ids = cg_mol->getBeadIds();
  sort(bead_ids.begin(),bead_ids.end());
  cout << "Grabbed bead_ids of molecule " << bead_ids.size() << endl;
  // Beads that have already been applied
  unordered_set<string> expired_beads;
  for ( int &cg_bead_id : bead_ids ){ 
    Bead * cg_bead = cg_top.getBead(cg_bead_id);
    // get the cg bead type
    string cg_bead_type = cg_bead->getType();
    cout << "Grabbed bead type " << cg_bead_type << endl;
    vector<string> bead_names = bead_type_and_names_[cg_bead_type];
    string bead_name;
    for( const string & name : bead_names){
      if(expired_beads.count(name)==false)  {
        cout << "Chosen name " << name << endl;
        bead_name = name;
        break;
      }
    } 
    expired_beads.insert(bead_name);

    map<string,Bead *> atomic_names_and_beads;
    cout << "bead id is " << cg_bead_id << endl;
    for ( pair<string,int> & atom_name_id : cgmolid_cgbeadid_atomicbeadids.second[cg_bead_id]){
      atomic_names_and_beads[atom_name_id.first] = atom_top.getBead(atom_name_id.second);
      cout << "Passing in atom names and ids " << atom_name_id.first << endl;
    }
    // Grab the correct map
    cout << "Bead name is " << bead_name << endl;
    assert(cg_bead_name_and_maps_.count(bead_name) && "Map for the coarse grained bead type is not known.");
    cout << "Applying to bead name " << bead_name << " found " << cg_bead_name_and_maps_.count(bead_name) << endl;
    cg_bead_name_and_maps_.at(bead_name)->Apply(cg_top.getBoundaryCondition(),
        atomic_names_and_beads,cg_bead);
  } 
}


}  // namespace csg
}  // namespace votca
