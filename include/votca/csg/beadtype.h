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

#ifndef _BEADTYPE_H
#define	_BEADTYPE_H

#include <memory>
#include <string>
#include "topologyitem.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
    \brief Bead Type informaton

    Each bead has a type. While the bead name should be unique,
    several beads can share the same type.
  */
class  BeadType : public TopologyItem {
public:    
    const int &getId() const { return _id; }
    const string &getName() const { return _name; }
    void setName(const string &name) { _name=name; }
    
private:
    int _id;
    string _name;
    
    BeadType(std::shared_ptr<Topology> parent, int id, const string &name)
    : TopologyItem(parent), _id(id), _name(name) {}
    friend class Topology;
};

}}

#endif	/* _BEADTYPE_H */

