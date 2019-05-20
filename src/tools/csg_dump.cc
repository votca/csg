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

#include "../../include/votca/csg/topology.h"
#include <iostream>
#include <map>
#include <string>
#include <votca/csg/boundarycondition.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/exclusionlist.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;
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
      Eigen::Matrix3d box = top->getBox();
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          cout << " " << box(i, j);
        }
        cout << endl << "           ";
      }
    }

    cout << "\nList of residues:" << endl;
    // Get all the residues by cycling through all the molecules
    vector<int> molecule_ids = top->getMoleculeIds();
    sort(molecule_ids.begin(), molecule_ids.end());
    for (const int &molecule_id : molecule_ids) {
      std::map<int, std::string> residue_ids_and_types =
          top->getResidueIdsAndTypesInMolecule(molecule_id);
      for (const pair<int, string> &id_and_type : residue_ids_and_types) {
        cout << "name: " << id_and_type.second;
        cout << " id: " << id_and_type.first << endl;
      }
    }

    cout << "\nList of molecules:\n";
    for (const int &molecule_id : molecule_ids) {
      Molecule *mol = &top->getMolecule(molecule_id);
      cout << "molecule: " << molecule_id + 1 << " " << mol->getType();
      cout << " beads: " << mol->BeadCount() << endl;

      vector<int> bead_ids = mol->getBeadIds();
      sort(bead_ids.begin(), bead_ids.end());
      for (const int &bead_id : bead_ids) {
        Bead *bead = top->getBead(bead_id);
        cout << bead_id << " Type " << bead->getType() << " Mass "
             << bead->getMass() << " Resnr " << bead->getResidueId()
             << " Resname " << bead->getResidueType() << " Charge "
             << bead->getQ() << endl;
      }
    }

  } else {
    cout << "\nList of exclusions:\n" << top->getExclusions();
  }

  return true;
}
