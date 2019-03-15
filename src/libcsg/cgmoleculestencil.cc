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

#include <list>
#include <map>
#include <string>
#include <vector>

#include "../../include/votca/csg/cgmoleculestencil.h"

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

namespace votca {
namespace csg {

  using namespace std;
  
  // This assumes that the vector is ordered according to the creation of the
  // beads in the molecule Such the the first cg_bead points to the atoms in the
  // atomistic molecule with the smallest ids And the second cg_bead in the
  // vector points to the beads in the atomistic molecule with the next largest
  // ids etc...
  void CGMoleculeStencil::AddBeadStencil(const vector<CGBeadStencil>& bead_info) {
    bead_info_ = bead_info;

    for (CGBeadStencil &info : bead_info_) {
      string cg_bead_name = info.cg_name_;
      for (string &atomic_name : info.atomic_subbeads_) {
        cg_and_atom_names_.insert(multi_bimap::value_type(cg_bead_name, atomic_name));
      }
    }
  }

  // Assumes that the bead_ids when sorted line up with the CGBeadStencil vector
  unordered_map<int, string> CGMoleculeStencil::MapAtomicBeadIdsToAtomicBeadNames(
      vector<int> bead_ids) const {
    assert(bead_ids.size() == cg_and_atom_names_.right.size() &&
           "number of bead_ids is not consistent with the number of atomic "
           "beads stored in the atomic molecule.");
    sort(bead_ids.begin(), bead_ids.end());
    unordered_map<int, string> id_and_bead_name;
    int index = 0;
    for (const CGBeadStencil & bead_info : bead_info_) {
      for (string atom_name : bead_info.atomic_subbeads_) {
        id_and_bead_name[bead_ids.at(index)] = atom_name;
        ++index;
      }
    }
    return id_and_bead_name;
  }

  unordered_map<int, string> CGMoleculeStencil::MapCGBeadIdsToCGBeadNames(
      vector<int> bead_ids) const {

    assert(bead_ids.size() == bead_info_.size() && "when mapping cg bead ids to cg bead names you must have the same number of beads in the molcule as in the stencil"); 
    sort(bead_ids.begin(),bead_ids.end());
    unordered_map<int,string> cg_id_and_names;
    int index = 0;
    for( const CGBeadStencil & stencil : bead_info_ ){
      cg_id_and_names[bead_ids.at(index)]=stencil.cg_name_; 
      ++index;
    }
    return cg_id_and_names;
  }

  vector<string> CGMoleculeStencil::getCGBeadNames() const{
    vector<string> cg_bead_names;
    for( const CGBeadStencil & stencil : bead_info_ ){
      cg_bead_names.push_back(stencil.cg_name_); 
    }
    return cg_bead_names;
  }

  vector<string> CGMoleculeStencil::getAtomicBeadNames(string cg_bead_name) const {
    vector<string> atom_names;
    std::pair<multi_bimap::left_const_iterator,
              multi_bimap::left_const_iterator>
        begin_and_end = cg_and_atom_names_.left.equal_range(cg_bead_name);
    for (multi_bimap::left_const_iterator iter = begin_and_end.first;
         iter != begin_and_end.second; ++iter) {
      atom_names.push_back(iter->second);
    }
    return atom_names;
  }

  string CGMoleculeStencil::getCGBeadName(string atom_bead_name) const {
    return cg_and_atom_names_.right.at(atom_bead_name);
  }

  void CGMoleculeStencil::AddInteractionStencil(const  vector<CGInteractionStencil> &interaction_info) {
    interaction_info_ = interaction_info;
  }

  const vector<CGBeadStencil> &CGMoleculeStencil::getBeadStencil() const { return bead_info_; }

  const vector<CGInteractionStencil> &CGMoleculeStencil::getInteractionStencil()const {
    return interaction_info_;
  }

  const std::string &CGMoleculeStencil::getCGMoleculeType() const { return cg_molecule_type_; }

  const std::string &CGMoleculeStencil::getAtomisticMoleculeType() const {
    return atomistic_molecule_type_;
  }

}  // namespace csg
}  // namespace votca

