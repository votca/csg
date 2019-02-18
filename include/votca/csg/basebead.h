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

#ifndef _VOTCA_CSG_BASEBEAD_H
#define _VOTCA_CSG_BASEBEAD_H

#include <assert.h>
#include <memory>

#include <votca/tools/identity.h>
#include <votca/tools/name.h>
#include <votca/tools/vec.h>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

/**
 * \brief information about a base bead
 *
 * The Base Bead class describes the core functionality of an atom or a coarse
 * grained bead. It stores information like the id, the name, the mass, the
 * charge and the residue it belongs to and the position
 *
 **/
class BaseBead {
 public:
  /**
   * destructor
   */
  virtual ~BaseBead() {}

  /// Gets the id of the bead
  int getId() const { return id_.getId(); }

  /// Sets the id of the bead
  void setId(int id) { id_.setId(id); }

  /// Sets the molecule the bead is attached too
  void setMoleculeId(int molecule_id) { molecule_id_.setId(molecule_id); }

  /// Gets the molecule pointer the bead is attached too
  int getMoleculeId() const { return molecule_id_.getId(); }

  /**
   * get the bead type
   * \return const string
   */
  virtual const std::string getType() const { return type_.getName(); }

  /**
   * set the bead type
   * \param bead type object
   */
  virtual void setType(std::string type) { type_.setName(type); }

  /**
   * get the mass of the base bead
   * \return - base bead mass
   */
  virtual const double &getMass() const { return mass_; }

  /**
   * set the mass of the base bead
   * \param - base bead mass
   */
  virtual void setMass(const double &m) { mass_ = m; }

  /**
   * set the position of the base bead
   * \param - base bead position
   */
  virtual void setPos(const TOOLS::vec &bead_position);

  /**
   * get the position of the base bead
   * \return base bead position
   */
  virtual const TOOLS::vec &getPos() const;

  /**
   * direct access (read/write) to the position of the base bead
   * \return reference to position
   */
  virtual TOOLS::vec &Pos() {
    assert(bead_position_set_ && "Position is not set.");
    return bead_position_;
  }

  /** does this configuration store positions? */
  bool HasPos() const { return bead_position_set_; }

  /** set has position to true */
  void HasPos(bool true_or_false) { bead_position_set_ = true_or_false; }

 protected:
  BaseBead() : mass_(0.0), bead_position_set_(false){};

  TOOLS::Identity<int> molecule_id_;
  TOOLS::Identity<int> id_;
  TOOLS::Name type_;
  TOOLS::Name element_symbol_;
  double mass_;
  TOOLS::vec bead_position_;

  bool bead_position_set_;
};

inline void BaseBead::setPos(const TOOLS::vec &bead_position) {
  bead_position_set_ = true;
  bead_position_ = bead_position;
}

inline const TOOLS::vec &BaseBead::getPos() const {
  assert(bead_position_set_ &&
         "Cannot get bead position as it has not been set.");
  return bead_position_;
}
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_BASEBEAD_H
