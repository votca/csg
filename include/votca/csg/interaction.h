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
#pragma once
#ifndef VOTCA_CSG_INTERACTION_H
#define VOTCA_CSG_INTERACTION_H

#include "boundarycondition.h"

#include <cassert>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

enum InteractionType { unassigned, bond, angle, dihedral };

std::string InteractionTypeToString(const InteractionType interaction_type);
/**
\brief base calss for all interactions

This is the base class for all interactions.

\todo double names/groups right, add molecules!!
*/
class Interaction {
 public:
  Interaction()
      : index_(-1),
        group_(""),
        group_id_(-1),
        interaction_type_(InteractionType::unassigned),
        mol_id_(tools::topology_constants::unassigned_molecule_id){};

  virtual ~Interaction() {}

  virtual std::unique_ptr<Interaction> Clone() const = 0;

  /**
   * @brief Determines the magnitude of the interaction
   *
   * @param bc
   * @param bead_positions provides the beads ids and pointers to their
   * positions
   *
   * @return
   */
  virtual double EvaluateVar(const BoundaryCondition &bc,
                             std::unordered_map<int, const Eigen::Vector3d *>
                                 bead_positions) const = 0;

  InteractionType getType() const { return interaction_type_; }

  void setGroup(const std::string &group) { group_ = group; }
  const std::string &getGroup() const {
    assert(group_.compare("") != 0);
    return group_;
  }

  int getGroupId() const {
    assert(group_id_ != -1 && "Cannot access group id as it has not been set");
    return group_id_;
  }
  void setGroupId(int id) { group_id_ = id; }

  void setIndex(const int &index) { index_ = index; }
  const int &getIndex() const {
    assert(index_ != -1 &&
           "Cannot access interaction index as it has not been set");
    return index_;
  }

  void setMoleculeId(const int &mol_id) { mol_id_ = mol_id; }
  const int &getMolecule() const {
    assert(mol_id_ != tools::topology_constants::unassigned_molecule_id &&
           "Cannot access interaction molecule id as it has not been set");
    return mol_id_;
  }

  virtual Eigen::Vector3d Grad(const BoundaryCondition &bc, int bead_index,
                               std::unordered_map<int, const Eigen::Vector3d *>
                                   bead_positions) const = 0;
  int BeadCount() { return bead_ids_.size(); }

  /**
   * @brief Given the bead_id index in the interaction vector return the id
   *
   * @param[in] bead_index value between 0 < BeadCount()
   *
   * @return the bead_ids id
   */
  int getBeadId(const int &bead_index) const {
    assert(bead_index > -1 &&
           boost::lexical_cast<size_t>(bead_index) < bead_ids_.size() &&
           "Cannot access interaction bead_id id as it has not been set");
    return bead_ids_.at(bead_index);
  }

  std::vector<int> getBeadIds() const { return bead_ids_; }

  /**
   * @brief Label returns a comprehensive descrition of the interaction
   *
   * The returned string is of the following form:
   * molecule id id:group name name:group id id:interaction type type:index
   * index
   *
   * So it might look like this
   * molcule id 1:group name Non-Bonded:group id 2:interaction type angle:index
   * 12
   * @return
   */
  std::string getLabel() const;

 protected:
  int index_;
  std::string group_;
  int group_id_;
  InteractionType interaction_type_;
  int mol_id_;
  std::vector<int> bead_ids_;
};

/**
    \brief bond interaction
*/
class IBond : public Interaction {
 public:
  // Constructors SHOULD ONLY BE CALLED BY Topology Object

  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IBond(*this));
  }

  double EvaluateVar(const BoundaryCondition &bc,
                     std::unordered_map<int, const Eigen::Vector3d *>
                         bead_positions) const override;
  Eigen::Vector3d Grad(const BoundaryCondition &bc, int bead_index,
                       std::unordered_map<int, const Eigen::Vector3d *>
                           bead_positions) const override;

 private:
  IBond(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 2 && "IBond must be called with 2 bead_ids.");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::bond;
  }
  friend class Topology;
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  // Constructors SHOULD ONLY BE CALLED BY Topology Object

  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IAngle(*this));
  }
  double EvaluateVar(const BoundaryCondition &bc,
                     std::unordered_map<int, const Eigen::Vector3d *>
                         bead_positions) const override;
  Eigen::Vector3d Grad(const BoundaryCondition &bc, int bead_index,
                       std::unordered_map<int, const Eigen::Vector3d *>
                           bead_positions) const override;

 private:
  IAngle(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 3 &&
           "Cannot create an IAngle with more or less than 3 bead_ids.");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::angle;
  }
  friend class Topology;
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  // Constructors SHOULD ONLY BE CALLED BY Topology Object
  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IDihedral(*this));
  }

  double EvaluateVar(const BoundaryCondition &bc,
                     std::unordered_map<int, const Eigen::Vector3d *>
                         bead_positions) const override;
  Eigen::Vector3d Grad(const BoundaryCondition &bc, int bead_index,
                       std::unordered_map<int, const Eigen::Vector3d *>
                           bead_positions) const override;

 private:
  IDihedral(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 4 &&
           "Cannot create a Dihedral with more or less than three bead_ids");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::dihedral;
  }
  friend class Topology;
};

inline double IBond::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  return bc
      .BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                            *bead_positions.at(bead_ids_.at(1)))
      .norm();
}

inline Eigen::Vector3d IBond::Grad(
    const BoundaryCondition &bc, int bead_index,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  Eigen::Vector3d r = bc.BCShortestConnection(
      *bead_positions.at(bead_ids_.at(0)), *bead_positions.at(bead_ids_.at(1)));
  r.normalize();
  return (bead_index == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  Eigen::Vector3d v1(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(0))));
  Eigen::Vector3d v2(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(2))));
  return std::acos(v1.dot(v2) / sqrt(v1.squaredNorm() * v2.squaredNorm()));
}

inline Eigen::Vector3d IAngle::Grad(
    const BoundaryCondition &bc, int bead_index,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  Eigen::Vector3d v1(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(0))));
  Eigen::Vector3d v2(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(2))));
  double acos_prime =
      1.0 / (sqrt(1 - std::pow(v1.dot(v2), 2) /
                          (v1.squaredNorm() * v2.squaredNorm())));

  switch (bead_index) {
    case (0):
      return acos_prime *
             (-v2 / (v1.norm() * v2.norm()) +
              (v1.dot(v2) * v1) / (v1.squaredNorm() * v2.squaredNorm()));
      break;
    case (1):
      return acos_prime *
             ((v1 + v2) / (v1.norm() * v2.norm()) -
              (v1.dot(v2)) * (v2.squaredNorm() * v1 + v1.squaredNorm() * v2) /
                  (std::pow(v1.norm(), 3) * std::pow(v2.norm(), 3)));
      break;
    case (2):
      return acos_prime * (-v1 / (v1.norm() * v2.norm())) +
             (v1.dot(v2) * v2 / (v1.norm() * std::pow(v2.norm(), 3)));
      break;
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}

inline double IDihedral::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  Eigen::Vector3d v1(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                              *bead_positions.at(bead_ids_.at(1))));
  Eigen::Vector3d v2(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(2))));
  Eigen::Vector3d v3(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(2)),
                              *bead_positions.at(bead_ids_.at(3))));
  Eigen::Vector3d n1 = v1.cross(v2);  // calculate the normal vector
  Eigen::Vector3d n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  return sign *
         std::acos(n1.dot(n2) / sqrt(n1.squaredNorm() * n2.squaredNorm()));
}

inline Eigen::Vector3d IDihedral::Grad(
    const BoundaryCondition &bc, int bead_index,
    std::unordered_map<int, const Eigen::Vector3d *> bead_positions) const {
  Eigen::Vector3d v1(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                              *bead_positions.at(bead_ids_.at(1))));
  Eigen::Vector3d v2(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                              *bead_positions.at(bead_ids_.at(2))));
  Eigen::Vector3d v3(
      bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(2)),
                              *bead_positions.at(bead_ids_.at(3))));

  Eigen::Vector3d n1 = v1.cross(v2);  // calculate the normal vector
  Eigen::Vector3d n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  Eigen::Vector3d returnvec;

  Eigen::Matrix3d e = Eigen::Matrix3d::Identity();

  double acos_prime =
      sign * (-1.0 / (sqrt(1 - std::pow(n1.dot(n2), 2) /
                                   (n1.squaredNorm() * n2.squaredNorm()))));
  switch (bead_index) {
    case (0): {  //
      for (int i = 0; i < 3; i++) {
        returnvec[i] = n2.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n1.dot(v2.cross(e.col(i))) /
                           (n2.norm() * std::pow(n1.norm(), 3));
      }
      return acos_prime * returnvec;
      break;
    }
    case (1): {
      for (int i = 0; i < 3; i++) {
        returnvec[i] =
            (n1.dot(v3.cross(e.col(i))) +
             n2.dot(e.col(i).cross(v1) + e.col(i).cross(v2))) /
                (n1.norm() * n2.norm()) -
            n1.dot(n2) * ((n1.dot(e.col(i).cross(v1) + e.col(i).cross(v2))) /
                              (n2.norm() * std::pow(n1.norm(), 3)) +
                          n2.dot(v3.cross(e.col(i))) /
                              (n1.norm() * std::pow(n2.norm(), 3)));
      }
      return acos_prime * returnvec;
      break;
    };
    case (2): {
      for (int i = 0; i < 3; i++) {
        returnvec[i] =
            (n1.dot(e.col(i).cross(v2) + e.col(i).cross(v3)) +
             n2.dot(v1.cross(e.col(i)))) /
                (n1.norm() * n2.norm()) -
            n1.dot(n2) * (n1.dot(v1.cross(e.col(i))) /
                              (n2.norm() * std::pow(n1.norm(), 3)) +
                          (n2.dot(e.col(i).cross(v2) + e.col(i).cross(v3))) /
                              (n1.norm() * std::pow(n2.norm(), 3)));
      }
      return acos_prime * returnvec;
      break;
    };
    case (3): {  //
      for (int i = 0; i < 3; i++) {
        returnvec[i] = n1.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n2.dot(v2.cross(e.col(i))) /
                           (n1.norm() * std::pow(n2.norm(), 3));
      }
      return acos_prime * returnvec;
      break;
    };
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_INTERACTION_H
