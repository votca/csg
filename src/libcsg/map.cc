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
#include <numeric>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/csgtopology.h>
#include <votca/csg/map.h>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/vec.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

Map::~Map() {
  vector<BeadMap *>::iterator iter;
  for (iter = _maps.begin(); iter != _maps.end(); ++iter) {
    delete (*iter);
  }
  _maps.clear();
}

void Map::Apply() {
  vector<BeadMap *>::iterator iter;
  for (iter = _maps.begin(); iter != _maps.end(); ++iter) {
    (*iter)->Apply();
  }
}

void Map_Sphere::Initialize(const CSG_Topology *topology_parent_,
                            const Molecule *mol_in, Bead *bead_out,
                            Property *opts_bead, Property *opts_map) {

  BeadMap::Initialize(topology_parent_, mol_in, bead_out, opts_bead, opts_map);

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

  //  for (size_t i = 0; i < beads.size(); ++i) {
  //    unordered_set<int> bead_ids = mol_in->getBeadIdsByLabel(beads.at(i));
  //    assert(bead_ids.size() == 1 &&
  //           "More than a single bead with the same label, maybe the globally
  //           " "unique bead id should be used instead.");

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
  vec cg(0., 0., 0.), f(0., 0., 0.), vel(0., 0., 0.);
  bool bPos, bVel, bF;
  bPos = bVel = bF = false;
  bead_out_->ClearParentBeads();

  // the following is needed for pbc treatment
  // CSG_Topology *top = bead_out_->getParent();
  double max_dist = 0.5 * topology_parent_->ShortestBoxSize();
  vec r0 = vec(0, 0, 0);
  string name0;
  int id0 = 0;
  if (_matrix.size() > 0) {
    if (_matrix.front().bead_in_->HasPos()) {
      r0 = _matrix.front().bead_in_->getPos();
      name0 = _matrix.front().bead_in_->getType();
      id0 = _matrix.front().bead_in_->getId();
    }
  }

  double M = 0;

  for (iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
    const Bead *bead = iter->bead_in_;
    bead_out_->AddParentBead(bead->getId());
    M += bead->getMass();
    if (bead->HasPos()) {
      vec r = topology_parent_->BCShortestConnection(r0, bead->getPos());
      if (abs(r) > max_dist) {
        cout << r0 << " " << bead->getPos() << endl;
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box \n (atoms " +
            name0 + " (id " + boost::lexical_cast<string>(id0 + 1) + ")" +
            ", " + bead->getType() + " (id " +
            boost::lexical_cast<string>(bead->getId() + 1) + ")" +
            +" , molecule " +
            boost::lexical_cast<string>(bead->getMoleculeId() + 1) + ")");
      }
      cg += (*iter)._weight * (r + r0);
      bPos = true;
    }
    if (bead->HasVel()) {
      vel += (*iter)._weight * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      f += (*iter)._force_weight * bead->getF();
      bF = true;
    }
  }
  bead_out_->setMass(M);
  if (bPos) bead_out_->setPos(cg);
  if (bVel) bead_out_->setVel(vel);
  if (bF) bead_out_->setF(f);
}

/// \todo implement this function
void Map_Ellipsoid::Apply() {
  vector<element_t>::iterator iter;
  vec cg(0., 0., 0.), c(0., 0., 0.), f(0., 0., 0.), vel(0., 0., 0.);
  matrix m(0.);
  bool bPos, bVel, bF;
  bPos = bVel = bF = false;

  // the following is needed for pbc treatment
  // CSG_Topology *top = bead_out_->getParent();
  double max_dist = 0.5 * topology_parent_->ShortestBoxSize();
  vec r0 = vec(0, 0, 0);
  if (_matrix.size() > 0) {
    if (_matrix.front().bead_in_->HasPos()) {
      r0 = _matrix.front().bead_in_->getPos();
    }
  }

  int n;
  n = 0;
  bead_out_->ClearParentBeads();
  for (iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
    const Bead *bead = iter->bead_in_;
    bead_out_->AddParentBead(bead->getId());
    if (bead->HasPos()) {
      vec r = topology_parent_->BCShortestConnection(r0, bead->getPos());
      if (abs(r) > max_dist) {
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box");
      }
      cg += (*iter)._weight * (r + r0);
      bPos = true;
    }
    if (bead->HasVel() == true) {
      vel += (*iter)._weight * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      f += (*iter)._force_weight * bead->getF();
      bF = true;
    }

    if ((*iter)._weight > 0 && bead->HasPos()) {
      c += bead->getPos();
      n++;
    }
  }

  if (bPos) bead_out_->setPos(cg);
  if (bVel) bead_out_->setVel(vel);
  if (bF) bead_out_->setF(f);

  if (!_matrix[0].bead_in_->HasPos()) {
    bead_out_->setU(vec(1.0, 0, 0));
    bead_out_->setV(vec(.0, 1, 0));
    bead_out_->setW(vec(.0, 0, 1));
    return;
  }

  // calculate the tensor of gyration
  c = c / (double)n;
  for (iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
    if ((*iter)._weight == 0) continue;
    const Bead *bead = iter->bead_in_;
    vec v = bead->getPos() - c;

    // Normalize the tensor with 1/number_of_atoms_per_bead
    m[0][0] += v.getX() * v.getX() / (double)_matrix.size();
    m[0][1] += v.getX() * v.getY() / (double)_matrix.size();
    m[0][2] += v.getX() * v.getZ() / (double)_matrix.size();
    m[1][1] += v.getY() * v.getY() / (double)_matrix.size();
    m[1][2] += v.getY() * v.getZ() / (double)_matrix.size();
    m[2][2] += v.getZ() * v.getZ() / (double)_matrix.size();
  }
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];

  // calculate the eigenvectors
  matrix::eigensystem_t es;
  m.SolveEigensystem(es);

  vec u = es.eigenvecs[0];
  vec v = _matrix[1].bead_in_->getPos() - _matrix[0].bead_in_->getPos();
  v.normalize();

  bead_out_->setV(v);

  vec w = _matrix[2].bead_in_->getPos() - _matrix[0].bead_in_->getPos();
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
