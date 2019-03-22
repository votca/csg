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

#ifndef VOTCA_CSG_BEADLIST_H
#define VOTCA_CSG_BEADLIST_H

#include <votca/tools/vec.h>

#include "bead.h"
#include "csgtopology.h"
#include <list>
#include <string>

namespace votca {
namespace csg {

/**
    \brief Generate lists of beads

    This class generates a list of beads based on some criteria, currently
    only the bead type.

*/
class BeadList : public std::list<Bead *> {
 public:
  BeadList(){};
  ~BeadList() {}

  /// \brief Select all beads of type <select>
  int Generate(CSG_Topology &top, const std::string &select);
  /// \brief Select all beads of type <select> withn a radius <radius> of
  /// reference vector <ref>
  int GenerateInSphericalSubvolume(CSG_Topology &top, const std::string &select,
                                   tools::vec ref, double radius);

  /// Get the csg topology object
  const CSG_Topology *getCSGTopologyConst() const { return topology_; }
  CSG_Topology *getCSGTopology() const { return topology_; }

 private:
  CSG_Topology *topology_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADLIST_H
