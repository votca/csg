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
#include <stdexcept>
#include <votca/csg/molecule.h>
#include <iostream>

namespace votca { namespace csg {

void Molecule::AddBead(Bead *bead, const string &name)
{
    _beads.push_back(bead);
    _bead_names.push_back(name);
    _beadmap[name] = _beads.size()-1;

    bead->_mol = this;
}

Bead * Molecule::getBead(int bead){
  if(bead<0 || bead>=_beads.size()){
    throw invalid_argument("getBead bead out of range of vector");
  }
  return _beads[bead];
}

int Molecule::getBeadId(int bead){
  if(bead<0 || bead>=_beads.size()){
    throw invalid_argument("getBeadId bead out of range of vector");
  }
  return _beads[bead]->getId();
}

int Molecule::getBeadByName(const string &name)
{
    map<string, int>::iterator iter = _beadmap.find(name);
    if(iter == _beadmap.end()) {
      string err = "bead name " + name + " is unrecognized in the molecule";
      throw invalid_argument(err);
    }
    return _beadmap[name];
}

string Molecule::getBeadName(int bead){
  if(bead<0 || bead>=_beads.size()){
    throw invalid_argument("getBeadName bead out of range of vector");
  }
  return _bead_names[bead];

}

}}
