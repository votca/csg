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

// TODO: This code need lots of cleaning up! please do not look at anything in
// here!
//

#include "analysistool.h"
#include "bondedstatistics.h"
#include "stdanalysis.h"
#include "tabulatedpotential.h"
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <votca/csg/csgapplication.h>
#include <votca/tools/rangeparser.h>
#include <votca/tools/tokenizer.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

class CsgBoltzmann : public CsgApplication {
 public:
  string ProgramName() { return "csg_boltzmann"; }
  void HelpText(ostream &out) {
    out << "Performs tasks that are needed for simple boltzmann\n"
           "inversion in an interactive environment.";
  }
  bool DoTrajectory() { return true; }
  bool DoMapping() { return true; }

  void Initialize();
  bool EvaluateOptions();
  void Run();

  void InteractiveMode();
  bool EvaluateTopology(CSG_Topology *top_cg, CSG_Topology *top_atomistic);

 protected:
  ExclusionList *CreateExclusionList(CSG_Topology &top_atomistic,
                                     CSG_Topology &top_cg,
                                     Molecule &mol_atomistic, Molecule &mol_cg);
  BondedStatistics _bs;
};
void CsgBoltzmann::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("Special options")(
      "excl", boost::program_options::value<string>(),
      "write atomistic exclusion list to file");

  AddObserver(&_bs);
}

bool CsgBoltzmann::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  if (OptionsMap().count("excl")) {
    CheckRequired("cg", "excl options needs a mapping file");
  }
  return true;
}

bool CsgBoltzmann::EvaluateTopology(CSG_Topology *top_cg,
                                    CSG_Topology *top_atomistic) {
  if (OptionsMap().count("excl")) {
    ExclusionList *ex;
    if (top_atomistic->MoleculeCount() > 1)
      cout << "WARNING: cannot create exclusion list for topology with"
              "multiple molecules, using only first molecule\n";

    vector<int> molecule_ids = top_atomistic->getMoleculeIds();
    Molecule *mol_atomistic = &top_atomistic->getMolecule(molecule_ids.at(0));
    Molecule *mol_cg = &top_cg->getMolecule(molecule_ids.at(0));
    cout << "Writing exclusion list for atomistic molecule "
         << mol_atomistic->getType() << " in coarse grained representation "
         << mol_cg->getType() << endl;

    ex = CreateExclusionList(*top_atomistic, *top_cg, *mol_atomistic, *mol_cg);

    ofstream fl;
    fl.open(OptionsMap()["excl"].as<string>().c_str());
    fl << "# atomistic: " << mol_atomistic->getType()
       << " cg: " << mol_cg->getType()
       << " cgmap: " << OptionsMap()["cg"].as<string>() << endl;
    fl << *ex;
    fl.close();
    delete ex;
    return false;
  }
  return true;
}

ExclusionList *CsgBoltzmann::CreateExclusionList(CSG_Topology &top_atomistic,
                                                 CSG_Topology &top_cg,
                                                 Molecule &mol_atomistic,
                                                 Molecule &mol_cg) {
  ExclusionList *ex = new ExclusionList();
  // exclude all with all
  {
    list<Bead *> excl_list;
    vector<int> bead_ids = mol_atomistic.getBeadIds();
    for (int &bead_id : bead_ids) {
      excl_list.push_back(mol_atomistic.getBead(bead_id));
    }
    ex->ExcludeList(excl_list);
  }

  // remove exclusions from inside a mapped bead
  vector<int> cg_bead_ids = mol_cg.getBeadIds();
  for (int &bead_id : cg_bead_ids) {
    const vector<int> &parent_beads = mol_cg.getBead(bead_id)->ParentBeads();
    list<Bead *> excl_list;

    for (const int &parent_bead_id : parent_beads) {
      excl_list.push_back(top_atomistic.getBead(parent_bead_id));
    }
    ex->Remove(excl_list);
  }

  // remove exclusion which come from mol_atomistic topology and hence bonds and
  // angles
  for (size_t index_1 = 0; index_1 < cg_bead_ids.size() - 1; ++index_1) {
    for (size_t index_2 = index_1 + 1; index_2 < cg_bead_ids.size();
         ++index_2) {
      int bead_id_1 = cg_bead_ids.at(index_1);
      int bead_id_2 = cg_bead_ids.at(index_2);
      if (top_cg.getExclusions().IsExcluded(mol_cg.getBead(bead_id_1),
                                            mol_cg.getBead(bead_id_2))) {
        const vector<int> &parent_beads_w =
            mol_cg.getBead(bead_id_1)->ParentBeads();
        const vector<int> &parent_beads_v =
            mol_cg.getBead(bead_id_2)->ParentBeads();

        for (const int parent_bead_id_w : parent_beads_w) {
          for (const int parent_bead_id_v : parent_beads_v) {
            ex->RemoveExclusion(top_atomistic.getBead(parent_bead_id_w),
                                top_atomistic.getBead(parent_bead_id_v));
          }
        }
      }
    }
  }
  return ex;
}

void CsgBoltzmann::Run() {
  CsgApplication::Run();
  if (OptionsMap().count("excl")) return;
  InteractiveMode();
}

void CsgBoltzmann::InteractiveMode() {
  std::map<std::string, AnalysisTool *> cmds;
  TabulatedPotential tab;
  StdAnalysis std;
  tab.Register(cmds);
  std.Register(cmds);

  string help_text =
      "Interactive mode, expecting commands:\n"
      "help: show this help\n"
      "q: quit\n"
      "list: list all available bonds\n"
      "vals <file> <selection>: write values to file\n"
      "hist <file> <selection>: create histogram\n"
      "tab <file> <selection>: create tabulated potential\n"
      "autocor <file> <selection>: calculate autocorrelation, only one row "
      "allowed in selection!\n"
      "cor <file> <selection>: calculate correlations, first row is correlated "
      "with all other rows";

  cout << help_text << endl;

  while (1) {
    string line;
    cout << "> ";
    getline(cin, line);

    boost::trim(line);
    vector<string> args;
    Tokenizer tok(line, " \t");
    tok.ToVector(args);

    if (args.size() == 0) continue;

    string cmd = args.front();
    args.erase(args.begin());

    if (cmd == "q") break;

    std::map<string, AnalysisTool *>::iterator tool;
    if (cmd == "help") {
      if (args.size() == 0) {
        cout << help_text << endl;
        continue;
      }
      cmd = args.front();
      args.erase(args.begin());
      tool = cmds.find(cmd);
      if (tool == cmds.end()) {
        cout << "error, no help item found" << endl;
        continue;
      }
      tool->second->Help(cmd, args);
      cout << endl;
      continue;
    }

    tool = cmds.find(cmd);
    if (tool == cmds.end()) {
      cout << "error, command not found" << endl;
      continue;
    }

    tool->second->Command(_bs, cmd, args);
  }
}

int main(int argc, char **argv) {
  CsgBoltzmann app;
  app.Exec(argc, argv);
}
