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
#include <algorithm>
#include <votca/csg/exclusionlist.h>

namespace votca {
namespace csg {

using namespace std;

void ExclusionList::Clear(void) {
  vector<exclusion_t *>::iterator iter;

  for (iter = _exclusions.begin(); iter != _exclusions.end(); ++iter)
    delete *iter;
  _exclusions.clear();
}
/*
//void ExclusionList::CreateExclusions(Topology *top) {
  InteractionContainer &ic = top->BondedInteractions();
  InteractionContainer::iterator ia;

  for (ia = ic.begin(); ia != ic.end(); ++ia) {
    int beads_in_int = (*ia)->BeadCount();
    list<Bead *> l;

    for (int ibead = 0; ibead < beads_in_int; ibead++) {
      int ii = (*ia)->getBeadId(ibead);
      l.push_back(top->getBead(ii));
    }
    ExcludeList(l);
  }
}*/

bool ExclusionList::IsExcluded(const Bead *bead1, const Bead *bead2) const {
  if (bead1->getMoleculeId() != bead2->getMoleculeId()) return false;
  if (bead2->getId() < bead1->getId()) swap(bead1, bead2);

  for (const pair<Bead *, exclusion_t *> bead_excl : _excl_by_bead) {

    if (bead_excl.first->getId() == bead1->getId()) {
      exclusion_t *excl = _excl_by_bead.at(bead_excl.first);
      if (find(excl->_exclude.begin(), excl->_exclude.end(), bead2) !=
          excl->_exclude.end()) {
        return true;
      }
    }
  }
  return false;
}

bool compareAtomIdiExclusionList(const ExclusionList::exclusion_t *a,
                                 const ExclusionList::exclusion_t *b) {
  return a->_atom->getId() < b->_atom->getId();
}

bool compareAtomIdBeadList(const Bead *a, const Bead *b) {
  return a->getId() < b->getId();
}

std::ostream &operator<<(std::ostream &out, ExclusionList &exl) {
  sort(exl._exclusions.begin(), exl._exclusions.end(),
       compareAtomIdiExclusionList);

  vector<ExclusionList::exclusion_t *>::iterator ex;
  for (ex = exl._exclusions.begin(); ex != exl._exclusions.end(); ++ex) {
    (*ex)->_exclude.sort(compareAtomIdBeadList);
    list<Bead *>::iterator i;
    out << (int)((*ex)->_atom->getId()) + 1;
    for (i = (*ex)->_exclude.begin(); i != (*ex)->_exclude.end(); ++i) {
      out << " " << ((*i)->getId() + 1);
    }
    out << endl;
  }
  return out;
}

}  // namespace csg
}  // namespace votca
