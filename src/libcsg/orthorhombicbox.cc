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

#include <votca/csg/orthorhombicbox.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;
vec OrthorhombicBox::BCShortestConnection(const vec &r_i,
                                          const vec &r_j) const {

  std::cout << "vector 1 " << r_i.getX() << " " << r_i.getY() << " "
            << r_i.getZ() << std::endl;
  std::cout << "vector 2 " << r_j.getX() << " " << r_j.getY() << " "
            << r_j.getZ() << std::endl;
  vec r_ij;
  double a = box_.get(0, 0);
  double b = box_.get(1, 1);
  double c = box_.get(2, 2);
  std::cout << "Sides of box " << a << " " << b << " " << c << std::endl;
  r_ij = r_j - r_i;
  r_ij.setZ(r_ij.getZ() - c * round(r_ij.getZ() / c));
  r_ij.setY(r_ij.getY() - b * round(r_ij.getY() / b));
  r_ij.setX(r_ij.getX() - a * round(r_ij.getX() / a));
  return r_ij;
}

}  // namespace csg
}  // namespace votca
