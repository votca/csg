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

#ifndef VOTCA_CSG_BOUNDARYCONDITION_H
#define VOTCA_CSG_BOUNDARYCONDITION_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

/**
 * @brief Class keeps track of how the boundaries of the system are handled
 *
 * There are a total of 3 different boundaries:
 * open - no boundaries
 * orthorhombic - orthorombic boundaries
 * triclinic - triclinic boundaries
 *
 * This class enables the correct treatement of distances beteween topology
 * objects, such that distances accound for the periodic boundaries.
 */
class BoundaryCondition {

 public:
  virtual ~BoundaryCondition(){};

  /**
   * @brief Safe way to allow child classes to be copied
   *
   * The child classes must use the same method and override it with their type
   * for this to work.
   *
   * @return standard pointer to child class
   */
  virtual std::unique_ptr<BoundaryCondition> Clone() const = 0;
  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const tools::matrix &box) { box_ = box; };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const tools::matrix &getBox() const { return box_; };

  /**
   * get the volume of the box
   * \return box volume as double
   */
  virtual double BoxVolume() const;

  /**
   * get shortest connection vector between r_i and r_j with respect to the
   * (periodic) box \return shortest distance vector
   */
  virtual tools::vec BCShortestConnection(const tools::vec &r_i,
                                          const tools::vec &r_j) const {
    throw std::runtime_error("BCShortestConnection is not implemented.");
  }

  double getShortestBoxDimension() const;

  enum eBoxtype { typeAuto = 0, typeTriclinic, typeOrthorhombic, typeOpen };
  virtual eBoxtype getBoxType() const {
    throw std::runtime_error("getBoxType is not implemented.");
  }

 protected:
  tools::matrix box_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BOUNDARYCONDITION_H
