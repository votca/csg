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

#ifndef VOTCA_CSG_CGMOLECULEDEF_H
#define VOTCA_CSG_CGMOLECULEDEF_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "exclusionlist.h"
#include "map.h"
#include "molecule.h"
#include <votca/tools/property.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BoundaryCondition;
/**
 * @brief Definition of coarse grained molecule
 *
 * This class is to define a coarse grained molecule, which includes the
 * topology, mapping, ...
 */
class CGMoleculeDef {
 public:
  CGMoleculeDef() {}
  ~CGMoleculeDef();

  /**
   * @brief Creates a coarse grained molecule
   *
   * Based on an xml file loaded using the **Load** method this method will
   * generate a coarse grained molecule.
   *
   * @param[in] top
   *
   * @return
   */
  Molecule *CreateMolecule(CSG_Topology &top);
  Map *CreateMap(const BoundaryCondition *boundaries, const Molecule &in,
                 Molecule &out);

  void Load(std::string filename);

  const std::string &getCGType() { return cg_molecule_type_; }
  const std::string &getAtomisticType() { return atomistic_molecule_type_; }

 private:
  Property options_;

  struct beaddef_t {
    std::string cg_name_;
    std::string type_;
    byte_t symmetry_;
    std::string mapping_;
    //  std::vector<std::string> subbeads_;
    Property *options_;
  };

  // type of the coarse grained molecule
  std::string cg_molecule_type_;
  // name of the molecule to coarse grain
  std::string atomistic_molecule_type_;

  // beads of the cg molecule
  std::vector<beaddef_t *> beads_;
  std::map<std::string, beaddef_t *> beads_by_cg_name_;

  // mapping schemes
  std::map<std::string, Property *> maps_;

  std::list<Property *> bonded_;

  void ParseTopology(Property &options);
  void ParseBeads(Property &options);
  void ParseBonded(Property &options);
  void ParseMapping(Property &options);

  beaddef_t *getBeadByCGName(const std::string &cg_name);
  Property *getMapByName(const std::string &map_name);
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGMOLECULEDEF_H
