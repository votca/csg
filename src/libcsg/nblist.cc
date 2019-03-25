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

#include "../../include/votca/csg/csgtopology.h"
#include <iostream>
#include <votca/csg/nblist.h>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;

NBList::NBList() : _do_exclusions(false), _match_function(0) {
  setPairType<BeadPair>();
  SetMatchFunction(NBList::match_always);
}

NBList::~NBList() {
  // TODO: NBList destructor
  if (_match_function) delete _match_function;
}

void NBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions) {
  BeadList::iterator iter1;
  BeadList::iterator iter2;
  _do_exclusions = do_exclusions;

  if (list1.empty()) return;
  if (list2.empty()) return;

  assert(list1.getCSGTopology() == list2.getCSGTopology());
  CSG_Topology *top = list1.getCSGTopology();

  for (iter1 = list1.begin(); iter1 != list1.end(); ++iter1) {
    if (&list1 == &list2) {
      iter2 = iter1;
      ++iter2;
    } else
      iter2 = list2.begin();

    if (*iter1 == *iter2) continue;

    for (; iter2 != list2.end(); ++iter2) {
      Eigen::Vector3d u = (*iter1)->getPos();
      Eigen::Vector3d v = (*iter2)->getPos();

      Eigen::Vector3d r = top->BCShortestConnection(u, v);
      double d = r.norm();
      if (d < _cutoff) {
        if (_do_exclusions)
          if (top->getExclusionsConst().IsExcluded(*iter1, *iter2)) {
            continue;
          }
        if ((*_match_function)(*iter1, *iter2, r, d))
          if (!FindPair(*iter1, *iter2))
            AddPair(_pair_creator(*iter1, *iter2, r));
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
