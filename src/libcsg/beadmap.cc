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

#include "../../include/votca/csg/beadmap.h"
#include "../../include/votca/csg/boundarycondition.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <votca/csg/bead.h>
#include <votca/tools/tokenizer.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

vector<string> BeadMap::getAtomicBeadNames() const {
  vector<string> bead_names;
  for (pair<string, element_t> pr : matrix_) {
    bead_names.push_back(pr.first);
  }
  return bead_names;
}

void Map_Sphere::InitializeBeadMap(const vector<string> &subbeads,
                                   vector<double> weights, vector<double> ds) {

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
    if (weights.at(i) == 0 && ds.at(i) != 0) {
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
  for (string name : subbeads) {
    AddAtomisticBead(name, weights.at(index), fweights.at(index));
    ++index;
  }
}

void Map_Sphere::AddAtomisticBead(const std::string &atomic_bead_name,
                                  const double &weight,
                                  const double &force_weight) {

  element_t el;
  el.weight_ = weight;
  el.force_weight_ = force_weight;
  matrix_[atomic_bead_name] = el;
}

void Map_Sphere::UpdateCGBead(const BoundaryCondition *boundaries,
                              map<string, const Bead *> atomic_beads,
                              Bead *cg_bead) const {

  assert(cg_bead->getSymmetry() == 1 &&
         "Applying spherical bead map on a non spherical corase grained bead");

  assert(matrix_.size() == atomic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified.");
  const Bead *atom = atomic_beads.begin()->second;
  assert(atom->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");

  Eigen::Vector3d reference_position = atom->getPos();
  string bead_type = atom->getType();
  int bead_id = atom->getId();

  bool bVel, bF;
  bVel = bF = false;
  double sum_of_atomistic_mass = 0;
  Eigen::Vector3d weighted_sum_of_atomistic_pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d weighted_sum_of_atomistic_forces = Eigen::Vector3d::Zero();
  Eigen::Vector3d weighted_sum_of_atomistic_velocity = Eigen::Vector3d::Zero();

  for (const pair<const string, element_t> &name_and_element : matrix_) {
    atom = atomic_beads[name_and_element.first];
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    Eigen::Vector3d shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());

    if (boundaries->getBoxType() != BoundaryCondition::eBoxtype::typeOpen) {
      double max_dist = 0.5 * boundaries->getShortestBoxDimension();
      if ((shortest_distance_beween_beads).norm() > max_dist) {
        throw std::runtime_error(
            "coarse-grained atom is bigger than half the box \n (atoms " +
            bead_type + " (id " + boost::lexical_cast<string>(bead_id + 1) +
            ")" + ", " + atom->getType() + " (id " +
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
void Map_Ellipsoid::UpdateCGBead(const BoundaryCondition *boundaries,
                                 std::map<string, const Bead *> atomistic_beads,
                                 Bead *cg_bead) const {

  assert(
      cg_bead->getSymmetry() == 3 &&
      "Applying ellipsoidal bead map on a non ellipsoidal corase grained bead");
  assert(matrix_.size() == atomistic_beads.size() &&
         "Cannot apply mapping mismatch in the number of atomistic beads in "
         "Map_Ellipsoid");
  assert(matrix_.size() > 0 &&
         "Cannot apply mapping no atomistic beads have been specified in "
         "Map_Ellipsoid.");

  std::map<string, const Bead *>::iterator name_and_bead_iter;
  name_and_bead_iter = atomistic_beads.begin();
  assert(name_and_bead_iter->second->HasPos() &&
         "Cannot apply mapping atomistic beads do not have position.");
  Eigen::Vector3d reference_position = name_and_bead_iter->second->getPos();
  ++name_and_bead_iter;

  bool bVel, bF;
  bVel = bF = false;
  double sum_of_atomistic_mass = 0;
  Eigen::Vector3d sum_of_atomistic_pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d weighted_sum_of_atomistic_pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d weighted_sum_of_atomistic_forces = Eigen::Vector3d::Zero();
  Eigen::Vector3d weighted_sum_of_atomistic_vel = Eigen::Vector3d::Zero();
  for (pair<string, element_t> name_and_element : matrix_) {
    const Bead *atom = atomistic_beads[name_and_element.first];
    sum_of_atomistic_mass += atom->getMass();
    assert(atom->HasPos() &&
           "Cannot apply mapping atomistic beads do not have position.");
    Eigen::Vector3d shortest_distance_beween_beads =
        boundaries->BCShortestConnection(reference_position, atom->getPos());
    if (boundaries->getBoxType() != BoundaryCondition::eBoxtype::typeOpen) {
      double max_dist = 0.5 * boundaries->getShortestBoxDimension();
      if (shortest_distance_beween_beads.norm() > max_dist) {
        throw std::runtime_error(
            "coarse-grained atom is bigger than half the box");
      }
    }
    assert(atom->HasPos() &&
           "Cannot map to coarse grained bead as atomic beads are missing "
           "positions");
    assert(name_and_element.second.weight_ &&
           "Cannot map to coarse grained beads with weights of 0.0");
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
  Eigen::Matrix3d tensor_of_gyration = Eigen::Matrix3d::Zero();
  double number_of_atom_beads = static_cast<double>(matrix_.size());
  Eigen::Vector3d average_pos =
      weighted_sum_of_atomistic_pos / number_of_atom_beads;
  for (pair<string, element_t> name_and_element : matrix_) {
    if (name_and_element.second.weight_ == 0) continue;
    const Bead *atom = atomistic_beads[name_and_element.first];
    Eigen::Vector3d pos = atom->getPos() - average_pos;

    // Normalize the tensor with 1/number_of_atoms_per_bead
    tensor_of_gyration += pos * pos.transpose() / number_of_atom_beads;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig;
  eig.computeDirect(tensor_of_gyration);

  Eigen::Vector3d u = eig.eigenvectors().col(0);
  Eigen::Vector3d reference_position2 = name_and_bead_iter->second->getPos();
  Eigen::Vector3d v = reference_position2 - reference_position;
  v.normalize();
  cg_bead->setV(v);

  ++name_and_bead_iter;
  Eigen::Vector3d reference_position3 = name_and_bead_iter->second->getPos();
  Eigen::Vector3d w = reference_position3 - reference_position;
  w.normalize();

  if (v.cross(w).dot(u) < 0) {
    u = -u;
  }
  cg_bead->setU(u);

  // write out w
  w = u.cross(v);
  w.normalize();
  cg_bead->setW(w);
}

AtomToCGMoleculeMapper::~AtomToCGMoleculeMapper() {
  cg_bead_name_and_maps_.clear();
}

// cg_bead_order contains the order of beads names
void AtomToCGMoleculeMapper::InitializeMoleculeMap(
    const unordered_map<string, CGBeadStencil> &bead_maps_info,
    const vector<string> &cg_bead_order) {

  for (string cg_bead_name : cg_bead_order) {
    const CGBeadStencil &bead_info = bead_maps_info.at(cg_bead_name);
    int symmetry = static_cast<int>(bead_info.cg_symmetry_);
    switch (symmetry) {
      case 1:
        cg_bead_name_and_maps_.insert(
            make_pair(cg_bead_name, unique_ptr<Map_Sphere>(new Map_Sphere())));
        break;
      case 3:
        cg_bead_name_and_maps_.insert(make_pair(
            cg_bead_name, unique_ptr<Map_Ellipsoid>(new Map_Ellipsoid())));
        break;
      default:
        throw runtime_error("unknown symmetry in bead definition!");
    }

    bead_type_and_names_[bead_info.cg_bead_type_].push_back(cg_bead_name);
    cg_bead_name_and_maps_.at(cg_bead_name)
        ->InitializeBeadMap(bead_info.atomic_subbeads_,
                            bead_info.subbead_weights_);
  }
}

/// Copy Constructor
AtomToCGMoleculeMapper::AtomToCGMoleculeMapper(
    const AtomToCGMoleculeMapper &other) {
  for (const std::pair<const std::string, std::unique_ptr<BeadMap>> &pr :
       other.cg_bead_name_and_maps_) {
    cg_bead_name_and_maps_.at(pr.first) = pr.second->Clone();
  }
  atom_molecule_type_ = other.atom_molecule_type_;
  cg_molecule_type_ = other.cg_molecule_type_;
}

/// Move Constructor
AtomToCGMoleculeMapper::AtomToCGMoleculeMapper(
    AtomToCGMoleculeMapper &&other) noexcept {
  for (std::pair<const std::string, std::unique_ptr<BeadMap>> &pr :
       other.cg_bead_name_and_maps_) {
    cg_bead_name_and_maps_.at(pr.first) = std::move(pr.second);
  }
  atom_molecule_type_ = other.atom_molecule_type_;
  cg_molecule_type_ = other.cg_molecule_type_;
}

/// Move assignment
AtomToCGMoleculeMapper &AtomToCGMoleculeMapper::operator=(
    AtomToCGMoleculeMapper &&other) noexcept {
  if (this != &other) {
    cg_bead_name_and_maps_.clear();
    for (std::pair<const std::string, std::unique_ptr<BeadMap>> &pr :
         other.cg_bead_name_and_maps_) {
      cg_bead_name_and_maps_.at(pr.first) = std::move(pr.second);
    }
    atom_molecule_type_ = other.atom_molecule_type_;
    cg_molecule_type_ = other.cg_molecule_type_;
  }
  return *this;
}

/// Copy Assignment
AtomToCGMoleculeMapper &AtomToCGMoleculeMapper::operator=(
    const AtomToCGMoleculeMapper &other) {
  if (this != &other) {
    cg_bead_name_and_maps_.clear();
    for (const std::pair<const std::string, std::unique_ptr<BeadMap>> &pr :
         other.cg_bead_name_and_maps_) {
      cg_bead_name_and_maps_.at(pr.first) = pr.second->Clone();
    }
    atom_molecule_type_ = other.atom_molecule_type_;
    cg_molecule_type_ = other.cg_molecule_type_;
  }
  return *this;
}

void AtomToCGMoleculeMapper::UpdateCGMolecule(
    CSG_Topology &atom_top, CSG_Topology &cg_top,
    CGMolToAtom cgmolid_cgbeadid_atomicbeadids) {

  // First int cg_molecule id
  // map
  //   first int - the cg_bead_id
  //   vector pair
  //       string - atomic name
  //       int - global atomic bead id

  // Grab the correct cg molecule
  int cg_mol_id = cgmolid_cgbeadid_atomicbeadids.first;
  Molecule *cg_mol = cg_top.getMolecule(cg_mol_id);

  // Ensure that the type molecule type is correct
  assert(cg_mol->getType() == cg_molecule_type_ &&
         "Cannot convert to the molecule  is not of the correct type");

  // Cycle throught the cg beads in a cg molecule
  vector<int> cg_bead_ids = cg_mol->getBeadIds();
  sort(cg_bead_ids.begin(), cg_bead_ids.end());
  // Beads that have already been applied
  unordered_set<string> expired_beads;
  for (int &cg_bead_id : cg_bead_ids) {
    Bead *cg_bead = cg_top.getBead(cg_bead_id);
    // get the cg bead type
    string cg_bead_type = cg_bead->getType();
    vector<string> bead_names = bead_type_and_names_[cg_bead_type];
    string bead_name;
    for (const string &name : bead_names) {
      if (expired_beads.count(name) == false) {
        bead_name = name;
        break;
      }
    }
    expired_beads.insert(bead_name);

    map<string, const Bead *> atomic_names_and_beads = getAtomicNamesAndBeads_(
        atom_top, cgmolid_cgbeadid_atomicbeadids.second.at(cg_bead_id));

    // Grab the correct map
    assert(cg_bead_name_and_maps_.count(bead_name) &&
           "Map for the coarse grained bead type is not known.");

    cg_bead_name_and_maps_.at(bead_name)->UpdateCGBead(
        cg_top.getBoundaryCondition(), atomic_names_and_beads, cg_bead);
  }
}

map<string, const Bead *> AtomToCGMoleculeMapper::getAtomicNamesAndBeads_(
    const CSG_Topology &atom_top,
    const vector<pair<string, int>> &atomic_names_and_ids) const {

  map<string, const Bead *> atomic_names_and_beads;
  for (const pair<string, int> &atom_name_id : atomic_names_and_ids) {
    atomic_names_and_beads[atom_name_id.first] =
        atom_top.getBeadConst(atom_name_id.second);
  }
  return atomic_names_and_beads;
}

}  // namespace csg
}  // namespace votca
