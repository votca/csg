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

#include <votca/csg/molecule.h>
#include <iostream>

namespace votca { namespace csg {

/*void Molecule::AddBead(Bead *bead, const string &name)
{
    _beads.push_back(bead);
    _bead_names.push_back(name);
    _beadmap[name] = _beads.size()-1;

    bead->_mol = this;
}*/
/*
Bead * Molecule::getBead(int bead){
  if(_beads.size()>=bead || bead<0) throw invalid_argument("bead was not found"); 
  return _beads[bead];
}

int Molecule::getBeadId(int bead) {
  if(_beads.size()>=bead || bead<0) throw invalid_argument("bead was not found"); 
  return _beads[bead]->getId();
}

string Molecule::getBeadName(int bead){
  if(_bead_names.size()>=bead || bead<0) throw invalid_argument("bead was not found"); 
  return _bead_names[bead]; 
}

int Molecule::getBeadByName(const string &name)
{
    map<string, int>::iterator iter = _beadmap.find(name);
    if(iter == _beadmap.end()) {
        string err =  "cannot find: <" + name + "> in " + _name;
        throw invalid_argument(err);
    }
    return _beadmap[name];
}
*/
}}
