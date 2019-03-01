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

Map::~Map() { bead_maps_.clear(); }

BeadMap *Map::CreateBeadMap(const byte_t symmetry,
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

void Map::Apply() {
  for (unique_ptr<BeadMap> &map : bead_maps_) {
    map->Apply();
  }
}

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

void Map_Sphere::Apply() {
  vector<element_t>::iterator iter;
  bool bPos, bVel, bF;
  bPos = bVel = bF = false;
  // bead_out_->ClearParentBeads();

  // the following is needed for pbc treatment
  vec reference_position = vec(0, 0, 0);
  string bead_type;
  int bead_id = 0;
  if (matrix_.size() > 0) {
    if (matrix_.front().bead_in_->HasPos()) {
      reference_position = matrix_.front().bead_in_->getPos();
      bead_type = matrix_.front().bead_in_->getType();
      bead_id = matrix_.front().bead_in_->getId();
    }
  }

  double sum_of_atomistic_mass = 0;
  vec weighted_sum_of_atomistic_pos(0., 0., 0.);
  vec weighted_sum_of_atomistic_forces(0., 0., 0.);
  vec weighted_sum_of_atomistic_velocity(0., 0., 0.);

  double max_dist = 0.5 * boundaries_->getShortestBoxDimension();
  for (iter = matrix_.begin(); iter != matrix_.end(); ++iter) {
    const Bead *bead = iter->bead_in_;
    // This is not needed because the ids should be the same between the cg bead
    // and the atomistic bead
    //    bead_out_->AddParentBead(bead->getId());
    sum_of_atomistic_mass += bead->getMass();
    if (bead->HasPos()) {
      vec shortest_distance_beween_beads =
          boundaries_->BCShortestConnection(reference_position, bead->getPos());

      if (abs(shortest_distance_beween_beads) > max_dist) {
        cout << reference_position << " " << bead->getPos() << endl;
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box \n (atoms " +
            bead_type + " (id " + boost::lexical_cast<string>(bead_id + 1) +
            ")" + ", " + bead->getType() + " (id " +
            boost::lexical_cast<string>(bead->getId() + 1) + ")" +
            +" , molecule " +
            boost::lexical_cast<string>(bead->getMoleculeId() + 1) + ")");
      }
      weighted_sum_of_atomistic_pos +=
          (*iter).weight_ *
          (shortest_distance_beween_beads + reference_position);
      bPos = true;
    }
    if (bead->HasVel()) {
      weighted_sum_of_atomistic_velocity += (*iter).weight_ * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      weighted_sum_of_atomistic_forces += (*iter).force_weight_ * bead->getF();
      bF = true;
    }
  }
  bead_out_->setMass(sum_of_atomistic_mass);
  if (bPos) bead_out_->setPos(weighted_sum_of_atomistic_pos);
  if (bVel) bead_out_->setVel(weighted_sum_of_atomistic_velocity);
  if (bF) bead_out_->setF(weighted_sum_of_atomistic_forces);
}

/// \todo implement this function
void Map_Ellipsoid::Apply() {
  vector<element_t>::iterator iter;
  vec cg(0., 0., 0.), c(0., 0., 0.), f(0., 0., 0.), vel(0., 0., 0.);
  matrix m(0.);
  bool bPos, bVel, bF;
  bPos = bVel = bF = false;

  // the following is needed for pbc treatment
  // BoundaryCondition *top = bead_out_->getParent();
  double max_dist = 0.5 * boundaries_->getShortestBoxDimension();
  vec r0 = vec(0, 0, 0);
  if (matrix_.size() > 0) {
    if (matrix_.front().bead_in_->HasPos()) {
      r0 = matrix_.front().bead_in_->getPos();
    }
  }

  int n;
  n = 0;
  bead_out_->ClearParentBeads();
  for (iter = matrix_.begin(); iter != matrix_.end(); ++iter) {
    const Bead *bead = iter->bead_in_;
    bead_out_->AddParentBead(bead->getId());
    if (bead->HasPos()) {
      vec r = boundaries_->BCShortestConnection(r0, bead->getPos());
      if (abs(r) > max_dist) {
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box");
      }
      cg += (*iter).weight_ * (r + r0);
      bPos = true;
    }
    if (bead->HasVel() == true) {
      vel += (*iter).weight_ * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      f += (*iter).force_weight_ * bead->getF();
      bF = true;
    }

    if ((*iter).weight_ > 0 && bead->HasPos()) {
      c += bead->getPos();
      n++;
    }
  }

  if (bPos) bead_out_->setPos(cg);
  if (bVel) bead_out_->setVel(vel);
  if (bF) bead_out_->setF(f);

  if (!matrix_[0].bead_in_->HasPos()) {
    bead_out_->setU(vec(1.0, 0, 0));
    bead_out_->setV(vec(.0, 1, 0));
    bead_out_->setW(vec(.0, 0, 1));
    return;
  }

  // calculate the tensor of gyration
  c = c / (double)n;
  for (iter = matrix_.begin(); iter != matrix_.end(); ++iter) {
    if ((*iter).weight_ == 0) continue;
    const Bead *bead = iter->bead_in_;
    vec v = bead->getPos() - c;

    // Normalize the tensor with 1/number_of_atoms_per_bead
    m[0][0] += v.getX() * v.getX() / (double)matrix_.size();
    m[0][1] += v.getX() * v.getY() / (double)matrix_.size();
    m[0][2] += v.getX() * v.getZ() / (double)matrix_.size();
    m[1][1] += v.getY() * v.getY() / (double)matrix_.size();
    m[1][2] += v.getY() * v.getZ() / (double)matrix_.size();
    m[2][2] += v.getZ() * v.getZ() / (double)matrix_.size();
  }
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];

  // calculate the eigenvectors
  matrix::eigensystem_t es;
  m.SolveEigensystem(es);

  vec u = es.eigenvecs[0];
  vec v = matrix_[1].bead_in_->getPos() - matrix_[0].bead_in_->getPos();
  v.normalize();

  bead_out_->setV(v);

  vec w = matrix_[2].bead_in_->getPos() - matrix_[0].bead_in_->getPos();
  w.normalize();

  if ((v ^ w) * u < 0) u = vec(0., 0., 0.) - u;
  bead_out_->setU(u);

  // write out w
  w = u ^ v;
  w.normalize();
  bead_out_->setW(w);
}

}  // namespace csg
}  // namespace votca
