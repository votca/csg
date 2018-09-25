/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_MOLECULE_H
#define	_VOTCA_CSG_MOLECULE_H

#include <vector>
#include <map>
#include <string>
#include <assert.h>

#include "beadstructure.h"
#include "topologyitem.h"

//#include "bead.h"

namespace votca { namespace csg {

class Interaction;
class Bead;

/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule : public TopologyItem,
                 public BeadStructure,
                 public virtual Identity<int>,
                 public virtual Name
{
public:            
   
    const std::string getLabel() const { 
      return "Id "+std::to_string(getId())+":Molecule "+getName();
    }
    /// Add an interaction to the molecule
    /// This is seperate from a Connect Beads method, an interaction does not
    /// guarantee a bond as far as I know. Though I will need to check. 
    void AddInteraction(Interaction *ic) { _interactions.push_back(ic);}

    vector<Interaction *> Interactions() { return _interactions; }

    template<typename T>
    void setUserData(T *userdata) { _userdata = (void*)userdata; }

    template<typename T>
    T *getUserData() { return (T *)_userdata; }
    
private:
  
    vector<Interaction*> _interactions;
     
    void *_userdata;
    
    /// constructor
    Molecule(Topology *parent, int id, string name)
        : Name(name), Identity(id), TopologyItem(parent)
    {}

    friend class Topology;
};

}}

#endif	// _VOTCA_CSG_MOLECULE_H

