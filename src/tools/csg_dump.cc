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

#include <iostream>
#include <map>
#include <string>
#include <votca/csg/boundarycondition.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/exclusionlist.h>
#include <votca/csg/topology.h>
#include <votca/tools/matrix.h>

using namespace std;
using namespace votca::csg;

class CsgDumpApp : public CsgApplication {
  string ProgramName() { return "csg_dump"; }
  void HelpText(ostream &out) {
    out << "Print atoms that are read from topology file to help"
           " debugging atom naming.";
  }
  void Initialize() {
    CsgApplication::Initialize();
    AddProgramOptions("Specific options")(
        "excl", "  display exclusion list instead of molecule list");
  }

  bool EvaluateTopology(Topology *top, Topology *top_ref);

  bool DoMapping() { return true; }
  bool DoMappingDefault(void) { return false; }
};

int main(int argc, char **argv) {
  CsgDumpApp app;

  return app.Exec(argc, argv);
}

bool CsgDumpApp::EvaluateTopology(Topology *top, Topology *top_ref) {
  if (!OptionsMap().count("excl")) {
    cout << "Boundary Condition: ";
    if (top->getBoxType() == BoundaryCondition::typeAuto) {
      cout << "auto";
    } else if (top->getBoxType() == BoundaryCondition::typeTriclinic) {
      cout << "triclinic";
    } else if (top->getBoxType() == BoundaryCondition::typeOrthorhombic) {
      cout << "orthorhombic";
    } else if (top->getBoxType() == BoundaryCondition::typeOpen) {
      cout << "open";
    }
    cout << endl;
    if (top->getBoxType() != BoundaryCondition::typeOpen) {
      cout << " Box matix:";
      matrix box = top->getBox();
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          cout << " " << box[i][j];
        }
        cout << endl << "           ";
      }
    }

    cout << "\nList of residues:\n";
    map<int, set<pair<int, string>>> molecule_id_residue_ids_and_names =
        top->getResidueIdsAndNames();
    //    for (pair<const int, string> &id_and_name : residue_ids_and_names) {
    //      cout << id_and_name.first << " name: " << id_and_name.second
    //           << " id: " << id_and_name.first << endl;
    //    }
    for (const pair<const int, set<pair<int, string>>>
             &id_and_res_num_and_res_name : molecule_id_residue_ids_and_names) {
      int molecular_type_id = id_and_res_num_and_res_name.first;
      for (const pair<int, string> &res_num_and_res_name :
           id_and_res_num_and_res_name.second) {
        cout << molecular_type_id - 1
             << " name: " << res_num_and_res_name.second
             << " id: " << res_num_and_res_name.first << endl;
      }
    }

    cout << "\nList of molecules:\n";
    MoleculeContainer::iterator mol;
    for (mol = top->Molecules().begin(); mol != top->Molecules().end(); ++mol) {
      cout << "molecule: " << (*mol)->getId() + 1 << " " << (*mol)->getName()
           << " beads: " << (*mol)->BeadCount() << endl;

      vector<int> bead_ids = (*mol)->getBeadIds();

      // for (int i = 0; i < (*mol)->BeadCount(); ++i) {
      for (const int &bead_id : bead_ids) {
        int resnr = (*mol)->getBead(bead_id)->getResidueNumber();

        cout << bead_id << " Name " << (*mol)->getBeadName(bead_id) << " Type "
             << (*mol)->getBead(bead_id)->getType() << " Mass "
             << (*mol)->getBead(bead_id)->getMass() << " Resnr " << resnr
             << " Resname " << (*mol)->getBead(bead_id)->getResidueName()
             << " Charge " << (*mol)->getBead(bead_id)->getQ() << endl;
      }
    }
  } else {
    cout << "\nList of exclusions:\n" << top->getExclusions();
  }

  return true;
}
