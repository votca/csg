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
#ifndef _VOTCA_CSG_BEADPAIR_H
#define _VOTCA_CSG_BEADPAIR_H

#include <votca/tools/eigen.h>

namespace votca {
namespace csg {

class Bead;
/**
   \brief A particle pair

   This class defines a particle pair. The future plan is, that the Pair class
   can be overloaded and Particle list creates these inherited pairs.

 */

class BeadPair {
 public:
  BeadPair() {}
  BeadPair(Bead *bead1, Bead *bead2, Eigen::Vector3d r)
      : _pair(std::pair<Bead *, Bead *>(bead1, bead2)),
        _r(r),
        _dist(r.norm()) {}

  Bead *first() { return _pair.first; }
  Bead *second() { return _pair.second; }
  /// \brief the vector connecting two beads
  Eigen::Vector3d &r() { return _r; }
  /// \brief the distance of the beads
  double &dist() { return _dist; }

 protected:
  std::pair<Bead *, Bead *> _pair;

  Eigen::Vector3d _r;
  double _dist;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_BEADPAIR_H */
