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

#ifndef _VOTCA_CSG_INTERACTION_H
#define _VOTCA_CSG_INTERACTION_H

#include "bead.h"
#include "boundarycondition.h"

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

/**
    \brief base calss for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction {
 public:
  Interaction()
      : _index(-1), _group(""), _group_id(-1), _name(""), mol_id_(-1){};

  virtual ~Interaction() {}
  virtual double EvaluateVar(const BoundaryCondition &bc) = 0;

  std::string getName() const { return _name; }

  void setGroup(const std::string &group) {
    _group = group;
    RebuildName();
  }
  const std::string &getGroup() const {
    assert(_group.compare("") != 0);
    return _group;
  }

  int getGroupId() {
    assert(_group_id != -1);
    return _group_id;
  }
  void setGroupId(int id) { _group_id = id; }

  void setIndex(const int &index) {
    _index = index;
    RebuildName();
  }
  const int &getIndex() const {
    assert(_index != -1);
    return _index;
  }

  void setMoleculeId(const int &mol_id) {
    mol_id_ = mol_id;
    RebuildName();
  }
  const int &getMolecule() const {
    assert(mol_id_ != -1);
    return mol_id_;
  }

  virtual TOOLS::vec Grad(const BoundaryCondition &bc, int bead) = 0;
  int BeadCount() { return beads_.size(); }

  /**
   * @brief Given the bead index in the interaction vector return the id
   *
   * @param[in] bead_index value between 0 < BeadCount()
   *
   * @return the beads id
   */
  int getBeadId(const int &bead_index) const {
    assert(bead_index > -1 &&
           boost::lexical_cast<size_t>(bead_index) < beads_.size());
    return beads_[bead_index]->getId();
  }

  std::vector<int> getBeadIds() const {
    std::vector<int> bead_ids;
    for (const BaseBead *bead_ptr : beads_) {
      bead_ids.push_back(bead_ptr->getId());
    }
    return bead_ids;
  }

  enum interaction_type { bond, angle, dihedral };

 protected:
  int _index;
  std::string _group;
  int _group_id;
  std::string _name;
  int mol_id_;
  std::vector<const BaseBead *> beads_;

  void RebuildName();
};

inline void Interaction::RebuildName() {
  std::stringstream s;
  if (mol_id_ != -1) s << "molecule " << mol_id_;
  if (!_group.empty()) {
    s << ":" << _group;
    if (_group_id != -1) {
      s << " " << _group_id;
    }
  }
  if (_index != -1) s << ":index " << _index;
  _name = s.str();
}

/**
    \brief bond interaction
*/
class IBond : public Interaction {
 public:
  /*  IBond(int bead1, int bead2) {
      beads_.resize(2);
      beads_[0] = bead1;
      beads_[1] = bead2;
    }

    IBond(std::list<int> &beads) {
      assert(beads.size() >= 2);
      beads_.resize(2);
      for (int i = 0; i < 2; ++i) {
        beads_[i] = beads.front();
        beads.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object
  double EvaluateVar(const BoundaryCondition &bc);
  TOOLS::vec Grad(const BoundaryCondition &bc, int bead);

 private:
  IBond(std::vector<const BaseBead *> beads) {
    assert(beads.size() != 2 && "IBond must be called with 2 beads.");
    beads_ = beads;
  }
  friend class CSG_Topology;
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  /*  IAngle(int bead1, int bead2, int bead3) {
      beads_.resize(3);
      beads_[0] = bead1;
      beads_[1] = bead2;
      beads_[2] = bead3;
    }*/
  /*  IAngle(std::list<int> &beads) {
      assert(beads.size() >= 3);
      beads_.resize(3);
      for (int i = 0; i < 3; ++i) {
        beads_[i] = beads.front();
        beads.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object

  double EvaluateVar(const BoundaryCondition &bc);
  TOOLS::vec Grad(const BoundaryCondition &bc, int bead);

 private:
  IAngle(std::vector<const BaseBead *> beads) {
    assert(beads.size() == 3 &&
           "Cannot create an IAngle with more or less than 3 beads.");
    beads_ = beads;
  }
  friend class CSG_Topology;
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  /*  IDihedral(int bead1, int bead2, int bead3, int bead4) {
      beads_.resize(4);
      beads_[0] = bead1;
      beads_[1] = bead2;
      beads_[2] = bead3;
      beads_[3] = bead4;
    }
    IDihedral(std::list<int> &beads) {
      assert(beads.size() >= 4);
      beads_.resize(4);
      for (int i = 0; i < 4; ++i) {
        beads_[i] = beads.front();
        beads.pop_front();
      }
    }*/
  // SHOULD ONLY BE CALLED BY Topology Object

  double EvaluateVar(const BoundaryCondition &bc);
  TOOLS::vec Grad(const BoundaryCondition &bc, int bead);

 private:
  IDihedral(std::vector<const BaseBead *> beads) {
    assert(beads.size() == 4 &&
           "Cannot create a Dihedral with more or less than three beads");
    beads_ = beads;
  }
  friend class CSG_Topology;
};

inline double IBond::EvaluateVar(const BoundaryCondition &bc) {
  std::cout << "Shortest distance between beads " << beads_[0] << " and "
            << beads_[1] << std::endl;
  return abs(bc.BCShortestConnection(beads_[0]->getPos(), beads_[1]->getPos()));
}

inline TOOLS::vec IBond::Grad(const BoundaryCondition &bc, int bead) {
  TOOLS::vec r =
      bc.BCShortestConnection(beads_[0]->getPos(), beads_[1]->getPos());
  r.normalize();
  return (bead == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(const BoundaryCondition &bc) {
  TOOLS::vec v1(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[0]->getPos()));
  TOOLS::vec v2(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[2]->getPos()));
  return acos(v1 * v2 / sqrt((v1 * v1) * (v2 * v2)));
}

inline TOOLS::vec IAngle::Grad(const BoundaryCondition &bc, int bead) {
  TOOLS::vec v1(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[0]->getPos()));
  TOOLS::vec v2(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[2]->getPos()));

  double acos_prime =
      1.0 / (sqrt(1 - (v1 * v2) * (v1 * v2) /
                          (abs(v1) * abs(v2) * abs(v1) * abs(v2))));
  switch (bead) {
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

inline double IDihedral::EvaluateVar(const BoundaryCondition &bc) {
  TOOLS::vec v1(
      bc.BCShortestConnection(beads_[0]->getPos(), beads_[1]->getPos()));
  TOOLS::vec v2(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[2]->getPos()));
  TOOLS::vec v3(
      bc.BCShortestConnection(beads_[2]->getPos(), beads_[3]->getPos()));
  TOOLS::vec n1, n2;
  n1 = v1 ^ v2;  // calculate the normal vector
  n2 = v2 ^ v3;  // calculate the normal vector
  double sign = (v1 * n2 < 0) ? -1 : 1;
  return sign * acos(n1 * n2 / sqrt((n1 * n1) * (n2 * n2)));
}

inline TOOLS::vec IDihedral::Grad(const BoundaryCondition &bc, int bead) {
  TOOLS::vec v1(
      bc.BCShortestConnection(beads_[0]->getPos(), beads_[1]->getPos()));
  TOOLS::vec v2(
      bc.BCShortestConnection(beads_[1]->getPos(), beads_[2]->getPos()));
  TOOLS::vec v3(
      bc.BCShortestConnection(beads_[2]->getPos(), beads_[3]->getPos()));
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
  switch (bead) {
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

#endif  // _VOTCA_CSG_INTERACTION_H
