/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _BEADLIST_H
#define	_BEADLIST_H

#include <memory>
#include <string>
#include <list>
#include "topology.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
    \brief Generate lists of beads

    This class generates a list of beads based on some criteria, currently
    only the bead type.

*/

class BeadList
    : public list<std::shared_ptr<Bead>>
{
public:
    BeadList() {};
    ~BeadList() {}
    
    /// \brief Select all beads of type <select>
    int Generate(Topology &top, const string &select);
     /// \brief Select all beads of type <select> withn a radius <radius> of reference vector <ref>
    int GenerateInSphericalSubvolume(Topology &top, const string &select,  vec ref, double radius);
    
    Topology *getTopology() {return _topology; }
    
private:
    Topology *_topology;
    
};

}}

#endif	/* _BEADLIST_H */

