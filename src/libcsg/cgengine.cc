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
  map<string, CGMoleculeDef *>::iterator i;
  for (i = _molecule_defs.begin(); i != _molecule_defs.end(); ++i)
    delete (*i).second;
  _molecule_defs.clear();
}

/**
    \todo melts with different molecules
*/
TopologyMap *CGEngine::CreateCGTopology(CSG_Topology &top_in,
                                        CSG_Topology &top_out) {

  // Grab all the molecules
  const unordered_map<int, Molecule> &mols = top_in.Molecules();
  // Setup which topology is being added to (top_out) and which topology is
  // being used to determine what is to be added (top_in)
  TopologyMap *topology_map = new TopologyMap(&top_in, &top_out);
  for (const pair<int, Molecule> &id_and_molecule : mols) {
    // Grab a molecule
    const Molecule *mol_in = &(id_and_molecule.second);
    if (IsIgnored(mol_in->getType())) continue;
    CGMoleculeDef *mol_def = getMoleculeDef(mol_in->getType());
    if (!mol_def) {
      cout << "--------------------------------------\n"
           << "WARNING: unknown molecule \"" << mol_in->getType()
           << "\" with id " << mol_in->getId() << " in topology" << endl
           << "molecule will not be mapped to CG representation\n"
           << "Check weather a mapping file for all molecule exists, was "
              "specified in --cg "
           << "separated by ; and the ident tag in xml-file matches the "
              "molecule name\n"
           << "--------------------------------------\n";
      continue;
    }
    Molecule *mcg = mol_def->CreateMolecule(top_out);
    Map *map = mol_def->CreateMap(&top_in, *mol_in, *mcg);
    topology_map->AddMoleculeMap(map);
  }
  top_out.RebuildExclusions();
  return topology_map;
}

void CGEngine::LoadMoleculeType(string filename) {
  Tokenizer tok(filename, ";");
  Tokenizer::iterator iter;

  for (iter = tok.begin(); iter != tok.end(); ++iter) {
    CGMoleculeDef *mol_def = new CGMoleculeDef();
    string file = *iter;
    boost::trim(file);
    mol_def->Load(file);
    _molecule_defs[mol_def->getIdent()] = mol_def;
  }
}

}  // namespace csg
}  // namespace votca
