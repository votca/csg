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

#include "../../include/votca/csg/beadlist.h"
#include "../../include/votca/csg/bead.h"
#include "../../include/votca/csg/csgtopology.h"
#include "../../include/votca/csg/molecule.h"
#include <Eigen/Dense>
#include <string>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using namespace votca::tools;
using namespace std;

int BeadList::Generate(CSG_Topology &top, const string &select) {
  topology_ = &top;
  bool selectByName = false;
  string pSelect;  // parsed selection string
  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  vector<int> bead_ids = top.getBeadIds();
  for (int &bead_id : bead_ids) {
    Bead *bead_ptr = top.getBead(bead_id);
    if (!selectByName) {
      if (wildcmp(pSelect.c_str(), bead_ptr->getType().c_str())) {
        push_back(bead_ptr);
      }
    } else {
      throw runtime_error("Bead name is no longer supported");
    }
  }
  return size();
}

int BeadList::GenerateInSphericalSubvolume(CSG_Topology &top,
                                           const string &select,
                                           Eigen::Vector3d ref, double radius) {
  topology_ = &(top);
  bool selectByName = false;
  string pSelect;  // parsed selection string

  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  vector<int> bead_ids = top.getBeadIds();
  for (int &bead_id : bead_ids) {
    Bead *bead_ptr = top.getBead(bead_id);
    if (top.BCShortestConnection(ref, bead_ptr->getPos()).norm() > radius)
      continue;
    if (!selectByName) {
      if (wildcmp(pSelect.c_str(), bead_ptr->getType().c_str())) {
        push_back(bead_ptr);
      }
    } else {
      throw runtime_error("Bead name is no longer supported");
    }
  }
  return size();
}

}  // namespace csg
}  // namespace votca
