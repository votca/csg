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
#ifndef _VOTCA_CSG_TRICLINICBOX_H
#define _VOTCA_CSG_TRICLINICBOX_H

#include "boundarycondition.h"

namespace votca {
namespace csg {

class TriclinicBox : public BoundaryCondition {

 public:
  virtual std::unique_ptr<BoundaryCondition> Clone() const override {
    return std::unique_ptr<BoundaryCondition>(new TriclinicBox(*this));
  }

  Eigen::Vector3d BCShortestConnection(
      const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const override;

  eBoxtype getBoxType() const override { return typeTriclinic; }

 protected:
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_TRICLINICBOX_H */
