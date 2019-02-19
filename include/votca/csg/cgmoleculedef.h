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

#ifndef _VOTCA_CSG_CGMOLECULEDEF_H
#define _VOTCA_CSG_CGMOLECULEDEF_H

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

/**
    \brief definition of a coarse grained molecule

    This class is to define a coarse grained molecule, which includes the
   topology, mapping, ...

    \todo clean up this class, do the bonded interactions right!!!!
    \todo check for consistency of xml file, seperate xml parser and class!!
*/
class CGMoleculeDef {
 public:
  CGMoleculeDef() {}
  ~CGMoleculeDef();

  Molecule *CreateMolecule(CSG_Topology &top);
  Map *CreateMap(const CSG_Topology *topology, const Molecule &in,
                 Molecule &out);

  void Load(std::string filename);

  const std::string &getType() { return type_; }
  const std::string &getIdent() { return ident_; }

 private:
  Property options_;

  struct beaddef_t {
    int id_;
    std::string type_;
    byte_t symmetry_;
    int residue_id_;
    std::string mapping_;
    std::vector<std::string> subbeads_;
    Property *options_;
  };

  // name of the coarse grained molecule
  std::string type_;
  // name of the molecule to coarse grain
  std::string ident_;

  // beads of the cg molecule
  std::vector<beaddef_t *> beads_;
  std::map<std::string, beaddef_t *> beads_by_type_;

  // mapping schemes
  std::map<std::string, Property *> maps_;

  std::list<Property *> bonded_;

  void ParseTopology(Property &options);
  void ParseBeads(Property &options);
  void ParseBonded(Property &options);
  void ParseMapping(Property &options);

  beaddef_t *getBeadByType(const std::string &type);
  Property *getMapByType(const std::string &type);
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_CGMOLECULEDEF_H */
