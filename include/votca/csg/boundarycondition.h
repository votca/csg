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

#ifndef _VOTCA_CSG_BOUNDARYCONDITION_H
#define _VOTCA_CSG_BOUNDARYCONDITION_H

#include <stdexcept>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

class BoundaryCondition {

 public:
  virtual ~BoundaryCondition(){};

  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const TOOLS::matrix &box) { box_ = box; };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const TOOLS::matrix &getBox() const { return box_; };

  /**
   * get the volume of the box
   * \return box volume as double
   */
  virtual double BoxVolume() const;

  /**
   * get shortest connection vector between r_i and r_j with respect to the
   * (periodic) box \return shortest distance vector
   */
  virtual TOOLS::vec BCShortestConnection(const TOOLS::vec &r_i,
                                          const TOOLS::vec &r_j) const {
    throw std::runtime_error("BCShortestConnection is not implemented.");
  }

  enum eBoxtype { typeAuto = 0, typeTriclinic, typeOrthorhombic, typeOpen };
  virtual eBoxtype getBoxType() {
    throw std::runtime_error("getBoxType is not implemented.");
  }

 protected:
  TOOLS::matrix box_;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_BOUNDARYCONDITION_H */
