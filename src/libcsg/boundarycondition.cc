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

#include <vector>
#include <votca/csg/boundarycondition.h>

namespace votca {
namespace csg {

using namespace votca::tools;

double BoundaryCondition::BoxVolume() const {
  vec a = box_.getCol(0);
  vec b = box_.getCol(1);
  vec c = box_.getCol(2);
  return (a ^ b) * c;
}

double BoundaryCondition::getShortestBoxDimension() const {
  cout << "Getting box columns" << endl;
  assert(getBoxType() != eBoxtype::typeOpen &&
         "Cannot get the shortest dimension of the box because it is open");
  TOOLS::vec _box_a = box_.getCol(0);
  TOOLS::vec _box_b = box_.getCol(1);
  TOOLS::vec _box_c = box_.getCol(2);

  // create plane normals
  TOOLS::vec _norm_a = _box_b ^ _box_c;
  TOOLS::vec _norm_b = _box_c ^ _box_a;
  TOOLS::vec _norm_c = _box_a ^ _box_b;

  _norm_a.normalize();
  _norm_b.normalize();
  _norm_c.normalize();

  double la = _box_a * _norm_a;
  double lb = _box_b * _norm_b;
  double lc = _box_c * _norm_c;

  return std::min(la, std::min(lb, lc));
}

}  // namespace csg
}  // namespace votca
