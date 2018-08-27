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
   
    /// Add a bead to the molecule, we will disable this
    // 1. because Beadstructure already has a AddBead method
    // 2. It makes no sense to add a different name to the bead
    // this is simply confusing 
    //void AddBead(Bead *bead, const string &name);
    // Use instead:
    // void AddBead(BaseBead * bead);


    /// get the id of a bead in the molecule
    // This also does not make since because the beadstructure stores
    // the beads in a map, using the id of the bead, if you know the id
    // you know the index they are the same. This as it is is just confusing
    //int getBeadId(int bead) { return _beads[bead]->getId(); }

    // This is a useful function to add to the beadstructure as the it is likey
    // to be used by any thing that inherits the beadstructure, with the 
    // exception that it should return a vector of all beads with the same name
    // there is no guarantee that there is not more than one bead with the same
    // name
    //int getBeadIdByName(const string &name);
    // Use instead: 
    // vector<int> getIdsOfBeadsWithName(const string &name);
    
    /// find a bead by it's name, this is unsafe because there is nothing
    // This is basically the same as the above function because the id is the index
    // remove this one
    //int getBeadByName(const string &name);

    // Each bead has a unique Id so this is a safe option
    // However, it makes since to put it in the beastructure class instead
    //string getBeadName(int bead) {return _bead_names[bead]; }

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
    // maps a name to a bead id not needed
    //map<string, int> _beadmap;
  
    vector<Interaction*> _interactions;
     
    // id of the molecules
    // int _id; not needed comes from the inherited Identity class
    
    // name of the molecule
    //string _name; Not needed comes from the inherited Name class
    
    // the beads in the molecule
    // Not needed stored in the beadstructure
    //vector<Bead *> _beads;

    // Not needed the bead names are the names stored in each bead object
    // vector<string> _bead_names;

    void *_userdata;
    
    /// constructor
    Molecule(Topology *parent, int id, string name)
        : TopologyItem(parent), _id(id), _name(name)
    {}

    friend class Topology;
};
/*
inline int Molecule::getBeadIdByName(const string &name)
{
    int i = getBeadByName(name);
    if(i<0)
        return i;
    return _beads[i]->getId();
}*/

}}

#endif	// _VOTCA_CSG_MOLECULE_H

