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

#ifndef VOTCA_CSG_CGMOLECULESTENCIL_H
#define VOTCA_CSG_CGMOLECULESTENCIL_H

#include <unordered_map>
#include <string>
#include <vector>

#include "cgbeadinfo.h"
#include <votca/tools/property.h>
#include <votca/tools/types.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

namespace votca {
namespace csg {
  
typedef boost::bimap<boost::bimaps::set_of<std::string>,boost::bimaps::multiset_of<std::string>> multi_bimap;
/**
 * @brief Definition of coarse grained molecule
 *
 * This class is to define a coarse grained molecule, which includes the
 * topology, mapping, ...
 */
class CGMoleculeStencil {
 public:
  CGMoleculeStencil(std::string cg_molecule_type, std::string atomistic_molecule_type)
      : cg_molecule_type_(cg_molecule_type),
        atomistic_molecule_type_(atomistic_molecule_type){};

  ~CGMoleculeStencil(){};

  // This assumes that the vector is ordered according to the creation of the
  // beads in the molecule Such the the first cg_bead points to the atoms in the
  // atomistic molecule with the smallest ids And the second cg_bead in the
  // vector points to the beads in the atomistic molecule with the next largest
  // ids etc...
  void AddBeadInfo(std::vector<CGBeadInfo> bead_info);

  // Assumes that the bead_ids when sorted line up with the CGBeadInfo vector
  std::unordered_map<int, std::string> MapAtomicBeadIdsToAtomicBeadNames(
      std::vector<int> bead_ids);
  
  std::vector<std::string> getAtomicBeadNames(std::string cg_bead_name);

  std::string getCGBeadName(std::string atom_bead_name);

  void AddInteractionInfo(std::vector<CGInteractionInfo> interaction_info);

  const std::vector<CGBeadInfo> &getBeadInfo();

  const std::vector<CGInteractionInfo> &getInteractionInfo();

  const std::string &getCGMoleculeType() const;

  const std::string &getAtomisticMoleculeType() const;

 private:
  std::string cg_molecule_type_;
  std::string atomistic_molecule_type_;
  multi_bimap cg_and_atom_names_;
  // Must appear in the order it is read in from the .xml file
  std::vector<CGBeadInfo> bead_info_;
  std::vector<CGInteractionInfo> interaction_info_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGMOLECULESTENCIL_H
