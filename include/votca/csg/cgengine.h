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

#ifndef VOTCA_CSG_CGENGINE_H
#define VOTCA_CSG_CGENGINE_H

#include "cgmoleculedef.h"
#include "cgobserver.h"
#include "csgtopology.h"
#include "topologymap.h"
#include <boost/program_options.hpp>
#include <list>
#include <map>
#include <votca/tools/datacollection.h>

#include "cgengine.h"
#include "cgmoleculedef.h"
#include "molecule.h"
#include "nematicorder.h"
#include "topologyreader.h"
#include "trajectoryreader.h"
#include "trajectorywriter.h"
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;
/**
 * @brief Coarse graining engine
 *
 * This class manages the coarse graining, at the moment it does the
 * measurement stuff
 */
class CGEngine {
 public:
  CGEngine();
  ~CGEngine();

  /**
      create a coarse grained topolgy based on a given topology
  */
  TopologyMap *CreateCGTopology(CSG_Topology &in, CSG_Topology &out);

  /**
      load molecule type from file
  */
  void RegisterCGMolecules(std::string filename);

  //  AtomisiticToCGMoleculeConverter *getMolecularConverter(std::string name);

  /**
   * \brief Adds molecules that are to be ignored during the mapping process
   * \param molecule_type glob molecule_type for molecule molecule_type
   */
  void AddIgnore(std::string molecule_type) {
    _ignores.push_back(molecule_type);
  }

  /**
   * \brief checks whether molecule is ignored
   * \param ident identifyier of molecule
   * \return true if is ignored
   */
  bool IsIgnored(std::string molecule_type);

 private:
  // std::map<std::string, AtomisiticToCGMoleculeConverter *> _molecule_defs;
  AtomToCGConverter converter;

  std::list<std::string> _ignores;
};
/*
inline AtomisiticToCGMoleculeConverter
*CGEngine::getMolecularConverter(std::string name) { std::map<std::string,
AtomisiticToCGMoleculeConverter *>::iterator iter;

  // if there is only 1 molecule definition, don't care about the name
  if (_molecule_defs.size() == 1 && name == "unnamed") {
    return (*(_molecule_defs.begin())).second;
  }

  iter = _molecule_defs.find(name);
  assert(iter != _molecule_defs.end() && "Molecule definition does not exist");
  return (*iter).second;
}*/

inline bool CGEngine::IsIgnored(std::string molecule_type) {
  for (std::list<std::string>::iterator iter = _ignores.begin();
       iter != _ignores.end(); ++iter) {
    if (wildcmp(iter->c_str(), molecule_type.c_str())) return true;
  }
  return false;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGENGINE_H
