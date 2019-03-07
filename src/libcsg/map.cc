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
    unordered_map<string, BeadMapInfo> bead_maps_info) {

  for (pair<const string, BeadMapInfo> & bead_map_info : bead_maps_info) {
    switch (bead_map_info.second.cg_symmetry_) {
      case 1:
        bead_type_and_maps_[bead_map_info.second.cg_bead_type_] =
            unique_ptr<Map_Sphere>(new Map_Sphere());
      case 3:
        bead_type_and_maps_[bead_map_info.second.cg_bead_type_] =
            unique_ptr<Map_Ellipsoid>(new Map_Ellipsoid());
      default:
        throw runtime_error("unknown symmetry in bead definition!");
    }

    bead_type_and_maps_[bead_map_info.second.cg_bead_type_]->Initialize(
        bead_map_info.second.atomic_subbeads_, bead_map_info.second.subbead_weights_);
  }
}
/*
BeadMap *AtomisticToCGMoleculeMapper::CreateBeadMap(const byte_t symmetry,
                            const BoundaryCondition *boundaries,
                            const Molecule *mol_in, Bead *bead_out,
                            Property *opts_map, Property *opts_bead) {

  switch (symmetry) {
    case 1:
      bead_maps_.push_back(unique_ptr<Map_Sphere>(new Map_Sphere()));
      break;
    case 3:
      bead_maps_.push_back(unique_ptr<Map_Ellipsoid>(new Map_Ellipsoid()));
      break;
    default:
      throw runtime_error(string("unknown symmetry in bead definition!"));
  }
  bead_maps_.back()->Initialize(boundaries, mol_in, bead_out, opts_map,
                                opts_bead);
  return bead_maps_.back().get();
}
*/
void AtomisticToCGMoleculeMapper::Apply(map<string, Bead *> name_and_atomic_bead) {

  for (pair<string,unique_ptr<BeadMap>> & bead_type_and_map : bead_type_and_maps_) {
    bead_type_and_map.second->Apply(name_and_atomic_bead);
  }
  // for (unique_ptr<BeadMap> &map : bead_maps_) {
  //  map->Apply();
  //}
}
/*
void Map_Sphere::Initialize(vector<string> subbeads, vector<double> weights){
        Initialize(subbeads,weights,vector<double>);
}*/

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
  // Which atomistic beads correspond to which weight and which force weight
  /*
    vector<int> bead_ids = mol_in->getBeadIds();
    // WARNING this assumes that the bead weights are read in the same order
    // as the bead ids are allocated, this method is prone to error because it
    // does not guarantee that the weights and forces line up with the correct
    // beads
    sort(bead_ids.begin(), bead_ids.end());
    int index = 0;
    for (const int &bead_id : bead_ids) {

      // int bead_id = *bead_ids.begin();
      if (bead_id < 0) {
        throw std::runtime_error(string("mapping error: molecule " +
                                        beads[index] + " does not exist"));
      }
      AddElem(mol_in->getBeadConst(bead_id), weights[index], fweights[index]);
      ++index;
    }
  */
}
/*
void Map_Sphere::Initialize(const BoundaryCondition *boundaries,
                            const Molecule *mol_in, Bead *bead_out,
                            Property *opts_bead, Property *opts_map) {

  BeadMap::Initialize(boundaries, mol_in, bead_out, opts_bead, opts_map);

  vector<string> beads;
  vector<double> weights;
  vector<double> fweights;

  // get the beads
  string s(opts_bead_->get("beads").value());
  Tokenizer tok_beads(s, " \n\t");
  tok_beads.ToVector(beads);

  // get vector of weights
  Tokenizer tok_weights(opts_map_->get("weights").value(), " \n\t");
  tok_weights.ConvertToVector<double>(weights);

  // check weather weights and # beads matches
  if (beads.size() != weights.size()) {
    throw runtime_error(
        string("number of subbeads in " + opts_bead->get("name").as<string>() +
               " and number of weights in map " +
               opts_map->get("name").as<string>() + " do not match"));
  }
  // normalize the weights
  double norm = 1. / std::accumulate(weights.begin(), weights.end(), 0.);

  transform(weights.begin(), weights.end(), weights.begin(),
            bind2nd(multiplies<double>(), norm));
  // get the d vector if exists or initialize same as weights
  vector<double> d;
  if (opts_map_->exists("d")) {
    Tokenizer tok_weights(opts_map_->get("d").value(), " \n\t");
    tok_weights.ConvertToVector(d);
    // normalize d coefficients
    norm = 1. / std::accumulate(d.begin(), d.end(), 0.);
    transform(d.begin(), d.end(), d.begin(),
              bind2nd(multiplies<double>(), norm));
  } else {
    // initialize force-weights with weights
    d.resize(weights.size());
    copy(weights.begin(), weights.end(), d.begin());
  }

  // check weather number of d coeffs is correct
  if (beads.size() != d.size()) {
    throw runtime_error(
        string("number of subbeads in " + opts_bead->get("name").as<string>() +
               " and number of d-coefficients in map " +
               opts_map->get("name").as<string>() + " do not match"));
  }

  fweights.resize(weights.size());
  // calculate force weights by d_i/w_i
  for (size_t i = 0; i < weights.size(); ++i) {
    if (weights[i] == 0 && d[i] != 0) {
      throw runtime_error(
          "A d coefficient is nonzero while weights is zero in mapping " +
          opts_map->get("name").as<string>());
    }
    if (weights[i] != 0) {
      fweights[i] = d[i] / weights[i];
    } else {
      fweights[i] = 0;
    }
  }

  vector<int> bead_ids = mol_in->getBeadIds();
  // WARNING this assumes that the bead weights are read in the same order
  // as the bead ids are allocated, this method is prone to error because it
  // does not guarantee that the weights and forces line up with the correct
  // beads
  sort(bead_ids.begin(), bead_ids.end());
  int index = 0;
  for (const int &bead_id : bead_ids) {

    // int bead_id = *bead_ids.begin();
    if (bead_id < 0) {
      throw std::runtime_error(string("mapping error: molecule " +
                                      beads[index] + " does not exist"));
    }
    AddElem(mol_in->getBeadConst(bead_id), weights[index], fweights[index]);
    ++index;
  }
}
*/
void Map_Sphere::Apply(BoundaryCondition *boundaries,
                       map<string, Bead *> atomistic_beads, Bead *cg_bead) {

  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified.");
  //	vector<element_t>::iterator iter;
  // bead_out_->ClearParentBeads();

  // the following is needed for pbc treatment
  //  vec reference_position = vec(0, 0, 0);
  //  string bead_type;
  //  int bead_id = 0;
  // if (matrix_.size() > 0) {

  Bead *atom = atomistic_beads.begin()->second;
  assert(atom->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  //	if (matrix_.front().bead_in_->HasPos()) {
  //	if( atom->HasPos() ) {
  //		reference_position = matrix_.front().bead_in_->getPos();
  //     bead_type = matrix_.front().bead_in_->getType();
  //    bead_id = matrix_.front().bead_in_->getId();
  vec reference_position = atom->getPos();
  string bead_type = atom->getType();
  int bead_id = atom->getId();
  //	}
  //}

  bool bVel, bF;
  bVel = bF = false;
  double sum_of_atomistic_mass = 0;
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_velocity(0., 0., 0.);

  double max_dist = 0.5 * boundaries->getShortestBoxDimension();
  for (pair<string, element_t> &name_and_element : matrix_) {
    //    const Bead *bead = iter->bead_in_;
    atom = atomisitic_beads[name_and_element.first];
    // This is not needed because the ids should be the same between the cg bead
    // and the atomistic bead
    //    bead_out_->AddParentBead(bead->getId());
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    //   if (atom->HasPos()) {
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
    //    bPos = true;
    //  }
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

/// \todo implement this function
/// Warning the atomistic beads must be a map they cannot be an unordered_map
void Map_Ellipsoid::Apply(BoundaryCondition *boundaries,
                          std::map<string, Bead *> atomistic_beads,
                          Bead *cg_bead) {
  // vector<element_t>::iterator iter;
  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads in "
         "Map_Ellipsoid");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified in "
         "Map_Ellipsoid.");
  // vec cg(0., 0., 0.), c(0., 0., 0.), f(0., 0., 0.), vel(0., 0., 0.);

  // the following is needed for pbc treatment
  // BoundaryCondition *top = bead_out_->getParent();
  double max_dist = 0.5 * boundaries->getShortestBoxDimension();
  std::map<string, Bead *>::iterator name_and_bead_iter;
  name_and_bead_iter = atomistic_beads.begin();
  // Bead * atom = name_and_bead_iter->second;
  assert(name_and_bead_iter->second->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  vec reference_position = name_and_bead_iter->second->getPos();
  ++name_and_bead_iter;
  // if (matrix_.size() > 0) {
  //  if (matrix_.front().bead_in_->HasPos()) {
  // }
  //}

  int n;
  n = 0;
  // bead_out_->ClearParentBeads();
  // vec tensor_of_gyration(0.,0.,0.);
  bool bVel, bF;
  bVel = bF = false;
  vec sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_velocity(0., 0., 0.);
  // for (iter = matrix_.begin(); iter != matrix_.end(); ++iter) {
  for (pair<string, element_t> name_and_element : matrix_) {
    const Bead *atom = atomistic_beads[name_and_element.first];
    // bead_out_->AddParentBead(atom->getId());
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    // if (atom->HasPos()) {
    vec shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());
    if (abs(shortest_distance_beween_beads) > max_dist) {
      throw std::runtime_error(
          "coarse-grained atom is bigger than half the box");
    }
    weighted_sum_of_atomistic_pos +=
        name_and_element.second.weight_ *
        (shortest_distance_beween_beads + reference_position);
    // bPos = true;
    //}
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

  // if (bPos) bead_out_->setPos(weighted_sum_of_atomistic_pos);
  cg_bead->setPos(weighted_sum_of_atomistic_pos);
  if (bVel) cg_bead->setVel(weighted_sum_of_atomisitc_vel);
  if (bF) cg_bead->setF(weighted_sum_of_atomisitc_forces);

  // if (!matrix_[0].bead_in_->HasPos()) {
  // if (!matrix_[0].bead_in_->HasPos()) {
  //  cg_bead->setU(vec(1.0, 0, 0));
  //   cg_bead->setV(vec(.0, 1, 0));
  //  cg_bead->setW(vec(.0, 0, 1));
  // return;
  // }

  // calculate the tensor of gyration
  matrix tensor_of_gyration(0.);
  vec average_pos = sum_of_atomisitc_pos / (double)n;
  double number_of_atomistic_beads = static_cast<double>(matrix_.size());
  // for (iter = matrix_.begin(); iter != matrix_.end(); ++iter) {
  for (pair<string, element_t> name_and_element : matrix_) {
    if (name_and_element.second.weight_ == 0) continue;
    const Bead *atom = atomistic_beads[name_and_element.first];
    // const Bead *bead = iter->bead_in_;
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
