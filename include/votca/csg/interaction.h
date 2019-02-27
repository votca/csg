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

#ifndef VOTCA_CSG_INTERACTION_H
#define VOTCA_CSG_INTERACTION_H

#include "boundarycondition.h"

#include <cassert>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/tools/vec.h>
namespace TOOLS = votca::tools;

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
        mol_id_(TOOLS::topology_constants::unassigned_molecule_id){};

  virtual ~Interaction() {}

  virtual std::unique_ptr<Interaction> Clone() const = 0;

  virtual double EvaluateVar(
      const BoundaryCondition &bc,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const = 0;

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
    assert(mol_id_ != TOOLS::topology_constants::unassigned_molecule_id &&
           "Cannot access interaction molecule id as it has not been set");
    return mol_id_;
  }

  virtual TOOLS::vec Grad(
      const BoundaryCondition &bc, int bead_id,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const = 0;
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
  /*  IBond(int bead1, int bead2) {
      bead_ids_.resize(2);
      bead_ids_[0] = bead1;
      bead_ids_[1] = bead2;
    }

    IBond(std::list<int> &bead_ids) {
      assert(bead_ids.size() >= 2);
      bead_ids_.resize(2);
      for (int i = 0; i < 2; ++i) {
        bead_ids_[i] = bead_ids.front();
        bead_ids.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object

  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IBond(*this));
  }

  double EvaluateVar(
      const BoundaryCondition &bc,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;
  TOOLS::vec Grad(
      const BoundaryCondition &bc, int bead_id,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;

 private:
  IBond(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 2 && "IBond must be called with 2 bead_ids.");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::bond;
  }
  friend class CSG_Topology;
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  /*  IAngle(int bead1, int bead2, int bead3) {
      bead_ids_.resize(3);
      bead_ids_[0] = bead1;
      bead_ids_[1] = bead2;
      bead_ids_[2] = bead3;
    }*/
  /*  IAngle(std::list<int> &bead_ids) {
      assert(bead_ids.size() >= 3);
      bead_ids_.resize(3);
      for (int i = 0; i < 3; ++i) {
        bead_ids_[i] = bead_ids.front();
        bead_ids.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object

  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IAngle(*this));
  }
  double EvaluateVar(
      const BoundaryCondition &bc,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;
  TOOLS::vec Grad(
      const BoundaryCondition &bc, int bead_id,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;

 private:
  IAngle(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 3 &&
           "Cannot create an IAngle with more or less than 3 bead_ids.");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::angle;
  }
  friend class CSG_Topology;
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  /*  IDihedral(int bead1, int bead2, int bead3, int bead4) {
      bead_ids_.resize(4);
      bead_ids_[0] = bead1;
      bead_ids_[1] = bead2;
      bead_ids_[2] = bead3;
      bead_ids_[3] = bead4;
    }
    IDihedral(std::list<int> &bead_ids) {
      assert(bead_ids.size() >= 4);
      bead_ids_.resize(4);
      for (int i = 0; i < 4; ++i) {
        bead_ids_[i] = bead_ids.front();
        bead_ids.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object
  std::unique_ptr<Interaction> Clone() const override {
    return std::unique_ptr<Interaction>(new IDihedral(*this));
  }

  double EvaluateVar(
      const BoundaryCondition &bc,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;
  TOOLS::vec Grad(
      const BoundaryCondition &bc, int bead_id,
      std::unordered_map<int, const TOOLS::vec *> bead_positions) const;

 private:
  IDihedral(std::vector<int> bead_ids) {
    assert(bead_ids.size() == 4 &&
           "Cannot create a Dihedral with more or less than three bead_ids");
    bead_ids_ = bead_ids;
    interaction_type_ = InteractionType::dihedral;
  }
  friend class CSG_Topology;
};

inline double IBond::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  //  std::cout << "Shortest distance between bead_ids " << bead_ids_[0] << "
  //  and "
  //            << bead_ids_[1] << std::endl;
  return abs(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                                     *bead_positions.at(bead_ids_.at(1))));
}

inline TOOLS::vec IBond::Grad(
    const BoundaryCondition &bc, int bead_id,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  TOOLS::vec r = bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                                         *bead_positions.at(bead_ids_.at(1)));
  r.normalize();
  return (bead_id == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  TOOLS::vec v1(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(0))));
  TOOLS::vec v2(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(2))));
  return acos(v1 * v2 / sqrt((v1 * v1) * (v2 * v2)));
}

inline TOOLS::vec IAngle::Grad(
    const BoundaryCondition &bc, int bead_id,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  TOOLS::vec v1(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(0))));
  TOOLS::vec v2(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(2))));

  double acos_prime =
      1.0 / (sqrt(1 - (v1 * v2) * (v1 * v2) /
                          (abs(v1) * abs(v2) * abs(v1) * abs(v2))));
  switch (bead_id) {
    case (0):
      return acos_prime *
             (-v2 / (abs(v1) * abs(v2)) +
              (v1 * v2) * v1 / (abs(v2) * abs(v1) * abs(v1) * abs(v1)));
      break;
    case (1):
      return acos_prime *
             ((v1 + v2) / (abs(v1) * abs(v2)) -
              (v1 * v2) * ((v2 * v2) * v1 + (v1 * v1) * v2) /
                  (abs(v1) * abs(v1) * abs(v1) * abs(v2) * abs(v2) * abs(v2)));
      break;
    case (2):
      return acos_prime *
             (-v1 / (abs(v1) * abs(v2)) +
              (v1 * v2) * v2 / (abs(v1) * abs(v2) * abs(v2) * abs(v2)));
      break;
  }
  // should never reach this
  assert(false);
  return TOOLS::vec(0, 0, 0);
}
inline double IDihedral::EvaluateVar(
    const BoundaryCondition &bc,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  TOOLS::vec v1(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                                        *bead_positions.at(bead_ids_.at(1))));
  TOOLS::vec v2(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(2))));
  TOOLS::vec v3(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(2)),
                                        *bead_positions.at(bead_ids_.at(3))));
  TOOLS::vec n1, n2;
  n1 = v1 ^ v2;  // calculate the normal vector
  n2 = v2 ^ v3;  // calculate the normal vector
  double sign = (v1 * n2 < 0) ? -1 : 1;
  return sign * acos(n1 * n2 / sqrt((n1 * n1) * (n2 * n2)));
}

inline TOOLS::vec IDihedral::Grad(
    const BoundaryCondition &bc, int bead_id,
    std::unordered_map<int, const TOOLS::vec *> bead_positions) const {
  TOOLS::vec v1(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(0)),
                                        *bead_positions.at(bead_ids_.at(1))));
  TOOLS::vec v2(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(1)),
                                        *bead_positions.at(bead_ids_.at(2))));
  TOOLS::vec v3(bc.BCShortestConnection(*bead_positions.at(bead_ids_.at(2)),
                                        *bead_positions.at(bead_ids_.at(3))));
  TOOLS::vec n1, n2;
  n1 = v1 ^ v2;  // calculate the normal vector
  n2 = v2 ^ v3;  // calculate the normal vector
  double sign = (v1 * n2 < 0) ? -1 : 1;
  TOOLS::vec returnvec;                       // vector to return
  double returnvec0, returnvec1, returnvec2;  // components of the return vector
  TOOLS::vec e0(1, 0, 0);  // unit vector pointing in x-direction
  TOOLS::vec e1(0, 1, 0);  // unit vector pointing in y-direction
  TOOLS::vec e2(0, 0, 1);  // unit vector pointing in z-direction

  double acos_prime =
      (-1.0 / (sqrt(1 - (n1 * n2) * (n1 * n2) /
                            (abs(n1) * abs(n2) * abs(n1) * abs(n2))))) *
      sign;
  switch (bead_id) {
    case (0): {  //
      returnvec0 = acos_prime * ((n2 * (v2 ^ e0)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n1 * (v2 ^ e0))) /
                                     (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
      returnvec1 = acos_prime * ((n2 * (v2 ^ e1)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n1 * (v2 ^ e1))) /
                                     (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
      returnvec2 = acos_prime * ((n2 * (v2 ^ e2)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n1 * (v2 ^ e2))) /
                                     (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
      returnvec.setX(returnvec0);
      returnvec.setY(returnvec1);
      returnvec.setZ(returnvec2);
      return returnvec;
      break;
    }
    case (1): {  //
      returnvec0 =
          acos_prime *
          ((n1 * (v3 ^ e0) + n2 * ((e0 ^ v1) + (e0 ^ v2))) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * ((e0 ^ v1) + (e0 ^ v2))) /
                 (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * (v3 ^ e0)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec1 =
          acos_prime *
          ((n1 * (v3 ^ e1) + n2 * ((e1 ^ v1) + (e1 ^ v2))) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * ((e1 ^ v1) + (e1 ^ v2))) /
                 (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * (v3 ^ e1)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec2 =
          acos_prime *
          ((n1 * (v3 ^ e2) + n2 * ((e2 ^ v1) + (e2 ^ v2))) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * ((e2 ^ v1) + (e2 ^ v2))) /
                 (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * (v3 ^ e2)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec.setX(returnvec0);
      returnvec.setY(returnvec1);
      returnvec.setZ(returnvec2);
      return returnvec;
      break;
    };
    case (2): {  //
      returnvec0 =
          acos_prime *
          ((n1 * ((e0 ^ v2) + (e0 ^ v3)) + n2 * (v1 ^ e0)) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * (v1 ^ e0)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * ((e0 ^ v2) + (e0 ^ v3))) /
                 (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec1 =
          acos_prime *
          ((n1 * ((e1 ^ v2) + (e1 ^ v3)) + n2 * (v1 ^ e1)) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * (v1 ^ e1)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * ((e1 ^ v2) + (e1 ^ v3))) /
                 (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec2 =
          acos_prime *
          ((n1 * ((e2 ^ v2) + (e2 ^ v3)) + n2 * (v1 ^ e2)) /
               (abs(n1) * abs(n2)) -
           ((n1 * n2) *
            ((n1 * (v1 ^ e2)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
             (n2 * ((e2 ^ v2) + (e2 ^ v3))) /
                 (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
      returnvec.setX(returnvec0);
      returnvec.setY(returnvec1);
      returnvec.setZ(returnvec2);
      return returnvec;
      break;
    };
    case (3): {  //
      returnvec0 = acos_prime * ((n1 * (v2 ^ e0)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n2 * (v2 ^ e0))) /
                                     (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
      returnvec1 = acos_prime * ((n1 * (v2 ^ e1)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n2 * (v2 ^ e1))) /
                                     (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
      returnvec2 = acos_prime * ((n1 * (v2 ^ e2)) / (abs(n1) * abs(n2)) -
                                 ((n1 * n2) * (n2 * (v2 ^ e2))) /
                                     (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
      returnvec.setX(returnvec0);
      returnvec.setY(returnvec1);
      returnvec.setZ(returnvec2);
      return returnvec;
      break;
    };
  }
  // should never reach this
  assert(false);
  return TOOLS::vec(0, 0, 0);
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_INTERACTION_H
