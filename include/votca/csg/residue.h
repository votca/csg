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

#ifndef _residue_H
#define	_residue_H

#include <string>
#include "topologyitem.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;
    
/**
    \brief class for a residue
 
    The Residue class describes a residue. When reading in the atoms, all beads belong to a residue. Later on, the molecules
    can be organized into molecules based on their residue.

*/
class Residue : public TopologyItem
{
public:
   
    /// get the name of the residue
    const string &getName();

    /// get the name of the residue
    const int &getId() const { return _id; }

    private:
    int _id;
    string _name;
private:
        /// constructor
    Residue(std::shared_ptr<Topology> parent, int id, const string &name)
        : TopologyItem(parent), _id(id), _name(name)
    {}
    friend class Topology;
};

inline const string &Residue::getName()
{
    return _name;
}

}}

#endif	/* _residue_H */

