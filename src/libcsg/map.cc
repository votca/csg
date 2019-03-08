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
#include <votca/csg/map.h>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/vec.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

AtomisticToCGMoleculeMapper::~AtomisticToCGMoleculeMapper() {
  bead_type_and_maps_.clear();
}

void AtomisticToCGMoleculeMapper::Initialize(
    unordered_map<string, CGBeadInfo> bead_maps_info) {

  for (pair<const string, CGBeadInfo> & bead_info : bead_maps_info) {
    switch (bead_info.second.cg_symmetry_) {
      case 1:
        bead_type_and_maps_[bead_info.second.cg_bead_type_] =
            unique_ptr<Map_Sphere>(new Map_Sphere());
      case 3:
        bead_type_and_maps_[bead_info.second.cg_bead_type_] =
            unique_ptr<Map_Ellipsoid>(new Map_Ellipsoid());
      default:
        throw runtime_error("unknown symmetry in bead definition!");
    }

    bead_type_and_maps_[bead_info.second.cg_bead_type_]->Initialize(
        bead_info.second.atomic_subbeads_, bead_info.second.subbead_weights_);
  }
}

void AtomisticToCGMoleculeMapper::Apply(
    CSG_Topology &atom_top, 
    CSG_Topology& cg_top, 
    pair<int,map<int,vector<pair<string,int>>>> cg_mol_id_cg_bead_id_atomic_bead_ids ){

  // First int cg_molecule id
  // map
  //   first int - the cg_bead_id 
  //   vector pair
  //       string - atomic name 
  //       int - global atomic bead id


  // Grab the correct cg molecule
  int cg_mol_id = cg_mol_id_cg_bead_id_atomic_bead_ids.first;
  Molecule * cg_mol = cg_top.getMolecule(cg_mol_id);
 
  // Ensure that the type molecule type is correct
  assert(cg_mol->getType() == cg_molecule_type_ && string("Cannot convert to ")
        + string("the molecule with id ") + to_string(cg_mol_id)+
        string(" as it is not of the correct type"));
  
  // Cycle throught the cg beads in a cg molecule
  vector<int> bead_ids = cg_mol->getBeadIds();
  for ( int &cg_bead_id : bead_ids ){ 
    Bead * cg_bead = cg_top.getBead(cg_bead_id);
    // get the cg bead type
    string cg_bead_type = cg_bead->getType();

    map<string,Bead *> atomic_names_and_beads;
    for ( pair<string,int> & atom_name_id : cg_mol_id_cg_bead_id_atomic_bead_ids.second[cg_bead_id]){
      atomic_names_and_beads[atom_name_id.first] = atom_top.getBead(atom_name_id.second);
    }
    // Grab the correct map
    bead_type_and_maps_[cg_bead_type]->Apply(cg_top.getBoundaryCondition(),
        atomic_names_and_beads,cg_bead);
  } 
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
          "A d coefficient is nonzero while weights is zero in mapping " +
          opts_map->get("name").as<string>());
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
  }
  
}

void Map_Sphere::Apply(const BoundaryCondition *boundaries,
                       map<string, Bead *> atomistic_beads, Bead *cg_bead) {

  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified.");
    Bead *atom = atomistic_beads.begin()->second;
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

  double max_dist = 0.5 * boundaries->getShortestBoxDimension();
  for (pair<string, element_t> &name_and_element : matrix_) {
    atom = atomisitic_beads[name_and_element.first];
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    vec shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());

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
  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads in "
         "Map_Ellipsoid");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified in "
         "Map_Ellipsoid.");

  // the following is needed for pbc treatment
  double max_dist = 0.5 * boundaries->getShortestBoxDimension();
  std::map<string, Bead *>::iterator name_and_bead_iter;
  name_and_bead_iter = atomistic_beads.begin();
  assert(name_and_bead_iter->second->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  vec reference_position = name_and_bead_iter->second->getPos();
  ++name_and_bead_iter;

  int n;
  n = 0;
  bool bVel, bF;
  bVel = bF = false;
  vec sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_velocity(0., 0., 0.);
  for (pair<string, element_t> name_and_element : matrix_) {
    const Bead *atom = atomistic_beads[name_and_element.first];
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    vec shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());
    if (abs(shortest_distance_beween_beads) > max_dist) {
      throw std::runtime_error(
          "coarse-grained atom is bigger than half the box");
    }
    weighted_sum_of_atomistic_pos +=
        name_and_element.second.weight_ *
        (shortest_distance_beween_beads + reference_position);
    if (atom->HasVel() == true) {
      weighted_sum_of_atomisitc_vel +=
          name_and_element.second.weight_ * atom->getVel();
      bVel = true;
    }
    if (atom->HasF()) {
      weighted_sum_of_atomisitc_forces +=
          name_and_element.second.force_weight_ * atom->getF();
      bF = true;
    }

    if (name_and_element.second.weight_ > 0 && atom->HasPos()) {
      sum_of_atomisitc_pos += atom->getPos();
      n++;
    }
  }

  cg_bead->setPos(weighted_sum_of_atomistic_pos);
  if (bVel) cg_bead->setVel(weighted_sum_of_atomisitc_vel);
  if (bF) cg_bead->setF(weighted_sum_of_atomisitc_forces);

  // calculate the tensor of gyration
  matrix tensor_of_gyration(0.);
  vec average_pos = sum_of_atomisitc_pos / (double)n;
  double number_of_atomistic_beads = static_cast<double>(matrix_.size());
  for (pair<string, element_t> name_and_element : matrix_) {
    if (name_and_element.second.weight_ == 0) continue;
    const Bead *atom = atomistic_beads[name_and_element.first];
    vec pos = atom->getPos() - average_pos;

    // Normalize the tensor with 1/number_of_atoms_per_bead
    tensor_of_gyration[0][0] +=
        pos.getX() * pos.getX() / number_of_atomisitic_beads;
    tensor_of_gyration[0][1] +=
        pos.getX() * pos.getY() / number_of_atomisitic_beads;
    tensor_of_gyration[0][2] +=
        pos.getX() * pos.getZ() / number_of_atomisitic_beads;
    tensor_of_gyration[1][1] +=
        pos.getY() * pos.getY() / number_of_atomisitic_beads;
    tensor_of_gyration[1][2] +=
        pos.getY() * pos.getZ() / number_of_atomisitic_beads;
    tensor_of_gyration[2][2] +=
        pos.getZ() * pos.getZ() / number_of_atomisitic_beads;
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
  vec w = reference_position3 - reference_position;
  w.normalize();

  if ((v ^ w) * u < 0) u = vec(0., 0., 0.) - u;
  cg_bead->setU(u);

  // write out w
  w = u ^ v;
  w.normalize();
  cg_bead->setW(w);
}

}  // namespace csg
}  // namespace votca
