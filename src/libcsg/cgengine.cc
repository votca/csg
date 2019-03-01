/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except top_in compliance with the License.
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

#include <fstream>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using namespace std;

namespace po = boost::program_options;

CGEngine::CGEngine() {}

CGEngine::~CGEngine() {
  /*  map<string, AtomToCGConverter *>::iterator i;
    for (i = _molecule_defs.begin(); i != _molecule_defs.end(); ++i)
      delete (*i).second;
    _molecule_defs.clear();*/
}

/**
    \todo melts with different molecules
*/
AtomToCGConverter *CGEngine::CreateCGTopology(CSG_Topology &atomistic_top_in,
                                              CSG_Topology &cg_top_out) {

  assert(atomistic_top_in.getBoxType() == cg_top_out.getBoxType() &&
         "box types of topology in and out differ");
  // Grab all the molecules
  const unordered_map<int, Molecule> &atomistic_mols =
      atomistic_top_in.Molecules();
  // Setup which topology is being added to (cg_top_out) and which topology is
  // being used to determine what is to be added (atomistic_top_in)
  // TopologyMap *topology_map = new TopologyMap(&atomistic_top_in,
  // &cg_top_out);
  for (const pair<int, Molecule> &id_and_molecule : atomistic_mols) {
    // Grab a molecule
    const Molecule *atomistic_mol = &(id_and_molecule.second);
    if (IsIgnored(atomistic_mol->getType())) continue;
    // AtomToCGConverter *converter =
    // getMolecularConverter(atomistic_mol->getType());
    /*    if (!converter) {
          cout << "--------------------------------------\n"
               << "WARNING: unknown molecule \"" << atomistic_mol->getType()
               << "\" with id " << atomistic_mol->getId() << " in topology" <<
       endl
               << "molecule will not be mapped to CG representation\n"
               << "Check weather a mapping file for all molecule exists, was "
                  "specified in --cg "
               << "separated by ; and the ident tag in xml-file matches the "
                  "molecule name\n"
               << "--------------------------------------\n";
          continue;
        }*/
    converter.ConvertAtomisticMoleculeToCGAndAddToCGTopology(atomistic_mol,
                                                             cg_top_out);

    // AtomisticToCGMolecaleMapper *map =
    // converter->CreateMap(atomistic_top_in.getBoundaryCondition(),
    //                                *atomistic_mol, *mcg);
    // topology_map->AddMoleculeMap(map);
  }
  cg_top_out.RebuildExclusions();
  // return topology_map;
  return topology_map;
}

void CGEngine::RegisterCGMolecules(string filename) {
  Tokenizer tok(filename, ";");
  Tokenizer::iterator iter;

  for (iter = tok.begin(); iter != tok.end(); ++iter) {
    //    AtomToCGConverter *converter = new AtomToCGConverter();
    string file = *iter;
    boost::trim(file);
    converter.LoadConversionStencil(file);
    topology_map.LoadMap(file);
    //_molecule_defs[converter->getAtomisticType()] = converter;
  }
}

}  // namespace csg
}  // namespace votca
