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

#include "../../include/votca/csg/csgapplication.h"
#include <boost/format.hpp>
#include <fstream>
#include <iostream>

using namespace votca::csg;
using namespace std;
using boost::format;
using namespace votca::tools;

/**
    \brief class for writing dlpoly topology files

    This class encapsulates the dlpoly topology writing functions

*/

class DLPTopolApp : public CsgApplication {
 public:
  string ProgramName() { return "csg_dlptopol"; }
  void HelpText(ostream &out) {
    out << "Create a dlpoly topology template based on an existing (atomistic) "
           "topology and \n"
        << "a mapping xml-file. The created template file needs to be "
           "inspected and amended by the user!\n\n"
        << "Examples:\n"
        << "* csg_dlptopol --top .dlpf --out .dlpf --cg cg-map.xml\n  convert "
           "FIELD to FIELD_CGV using cg-map.xml\n"
        << "* csg_dlptopol --top FA-dlpoly.dlpf --out CG-dlpoly.dlpf --cg "
           "cg-map.xml\n"
        << "* csg_dlptopol --top FA-gromacs.tpr --out FA-dlpoly.dlpf "
           "--no-map\n";
  }

  bool DoMapping(void) { return true; }

  void Initialize(void);
  bool EvaluateOptions(void) {
    CsgApplication::EvaluateOptions();
    CheckRequired("out", "no output topology specified");
    return true;
  }
  bool EvaluateTopology(Topology *top, Topology *top_ref);

 protected:
  void WriteMoleculeAtoms(ostream &out, const Molecule &cg);
  void WriteMoleculeInteractions(ostream &out, const Molecule &cg);
  void WriteVDWInteractions(ostream &out, const Molecule &cg);
  void WriteMolecularType(ostream &out, const Molecule &cg, int nummol);
};

void DLPTopolApp::Initialize(void) {
  CsgApplication::Initialize();
  // AddProgramOptions()
  //("top", boost::program_options::value<string>(),
  //"  input topology in any known format:\n  <name>.dlpf for dlpoly, <name>.tpr
  // for gromacs\n  (convention: '.dlpf'='use FIELD')");
  AddProgramOptions()("out", boost::program_options::value<string>(),
                      "  output topology in dlpoly format");
}

bool DLPTopolApp::EvaluateTopology(Topology *top, Topology *top_ref) {
  // check the file names from options

  string fname = OptionsMap()["top"].as<string>();

  if (fname == ".dlpf") {
    fname = "FIELD";
  }

  cout << "input  file-name: " << fname << endl;

  fname = OptionsMap()["out"].as<string>();

#ifdef DEBUG
  cout << "output file-name given: " << fname << endl;
#endif

  if (fname == ".dlpf") {
    fname = "FIELD_CGV";
  }

#ifdef DEBUG
  cout << "output file-name actual: " << fname << endl;
#else
  cout << "output file-name: " << fname << endl;
#endif

  // do CG mapping

  std::vector<Molecule *> MolecularTypes;

  unordered_set<string> previous_molecule_types;
  unordered_map<string, int> molecules_and_count;
  vector<string> vdw_pairs;

  // find all unique molecular types
  vector<int> molecule_ids = top->getMoleculeIds();
  for (const int &molecule_id : molecule_ids) {
    Molecule *mol = &top->getMolecule(molecule_id);
    // molecules are ignored during the mapping stage
    // i.e. the ignored ones do not enter the CG topology (*top) - ?
    // if( IsIgnored(mol->getName()) ) continue;
    const string molecule_type = mol->getType();
    if (previous_molecule_types.count(molecule_type)) {
      ++molecules_and_count[molecule_type];
      // prv_mol_number++;
      continue;
    }

    molecules_and_count[molecule_type] = 1;
    previous_molecule_types.insert(mol->getType());

    //#ifdef DEBUG
    cout << "'" << mol->getType() << "' added to CG molecular types" << endl;
    //#endif

    MolecularTypes.push_back(mol);

    // collect unique bead pairs over all molecular types found
    vector<int> bead_ids = mol->getBeadIds();
    for (int &bead_id : bead_ids) {
      const string bead_name1 = mol->getBead(bead_id).getType();

      for (unsigned int imt = 0; imt < MolecularTypes.size(); imt++) {

        vector<int> bead_ids2 = MolecularTypes[imt]->getBeadIds();
        for (int &bead_id2 : bead_ids2) {

          const string bead_name2 =
              MolecularTypes[imt]->getBead(bead_id2).getType();

          stringstream ss_bp1, ss_bp2;

          ss_bp1 << format("%8s%8s") % bead_name1 % bead_name2;
          ss_bp2 << format("%8s%8s") % bead_name2 % bead_name1;
          bool is_new_pair = true;

          for (unsigned int ibp = 0; ibp < vdw_pairs.size(); ibp++) {
            if (ss_bp1.str() == vdw_pairs[ibp] ||
                ss_bp2.str() == vdw_pairs[ibp]) {
              is_new_pair = false;
              break;
            }
          }
          if (is_new_pair) {
            vdw_pairs.push_back(ss_bp1.str());
#ifdef DEBUG
            cout << "'" << ss_bp1.str() << "' added to CG vdw pairs" << endl;
#endif
          }
        }
      }
    }
  }

  if (MolecularTypes.size() > 1)
    cout << "WARNING: creation of topology for multiple molecular types "
            "is experimental at this stage\n";

  ofstream fl;
  fl.open(fname.c_str());

  fl << "From VOTCA with love"
     << " # check/amend this file if needed!\n";
  fl << "units kJ\n";
  fl << "molecular types " << MolecularTypes.size() << endl;

  for (unsigned int i = 0; i < MolecularTypes.size(); i++) {
    WriteMolecularType(fl, *(MolecularTypes[i]),
                       molecules_and_count[MolecularTypes[i]->getType()]);
  }

  if (vdw_pairs.size() > 0) {

    fl << "vdw " << vdw_pairs.size() << endl;

    for (unsigned int ibp = 0; ibp < vdw_pairs.size(); ibp++) {
      fl << vdw_pairs[ibp] << " tab   1.00000  0.00000\n";
    }
  }

  fl << "close" << endl;

  cout << "Created template for dlpoly topology - please, check & amend if "
          "needed!"
       << endl;

  fl.close();
  return true;
}

void DLPTopolApp::WriteMoleculeAtoms(ostream &out, const Molecule &cg) {
  out << "atoms " << cg.BeadCount() << endl;
  out << "# name  mass  charge  nrept  ifrozen (optional: ngroup, index, "
         "name/type, type/residue, index/res-ID) \n";
  vector<int> bead_ids = cg.getBeadIds();
  int index = 0;
  for (int &bead_id : bead_ids) {
    const Bead b = cg.getBead(bead_id);

    const string btype = b.getType();

    out << format("%8s  %10f  %10f     1     0     1 %10d  %8s  %10d \n") %
               btype % b.getMass() % b.getQ() % (index + 1) % btype %
               (index + 1);
    ++index;
  }
}

void DLPTopolApp::WriteMoleculeInteractions(ostream &out, const Molecule &cg) {
  vector<Interaction *> ics = cg.Interactions();
  vector<Interaction *>::iterator iter;

  stringstream sout;

  int n_entries = 0;
  int nb = -1;

  for (iter = ics.begin(); iter != ics.end(); ++iter) {
    Interaction *ic = *iter;
    if (nb != ic->BeadCount()) {

      if (sout.str() != "") out << n_entries << endl << sout.str();

      sout.str("");
      n_entries = 0;

      nb = ic->BeadCount();
      switch (nb) {
        case 2:
          out << "bonds ";
          break;
        case 3:
          out << "angles ";
          break;
        case 4:
          out << "dihedrals ";
          break;
        default:
          string err = "cannot handle number of beads in interaction:";
          err += to_string(ic->getMolecule() + 1) + ":" + ic->getGroup();
          err += ":" + to_string(ic->getIndex() + 1);
          throw runtime_error(err);
      }
    }
    n_entries++;
    // to do: is it possible to use bond/angle/dihedral function types for 1:1
    // mapping? (CG map overwrites ic->Group anyway)
    // sout << ic->getInteractionFunc(); // something like that (only for 1:1
    // mapping!)
    sout << " tab ";
    for (int i = 0; i < nb; ++i) sout << ic->getBeadId(i) + 1 << " ";
    sout << "   1.00000  0.00000"
         << " # ";
    sout << to_string(ic->getMolecule() + 1);
    sout << ":" + ic->getGroup();
    sout << ":" + to_string(ic->getIndex() + 1) << endl;
  }
  if (sout.str() != "") out << n_entries << endl << sout.str();
}

void DLPTopolApp::WriteMolecularType(ostream &out, const Molecule &cg,
                                     int nummol) {
  out << cg.getType() << endl;
  out << "nummols " << nummol << endl;

  WriteMoleculeAtoms(out, cg);
  WriteMoleculeInteractions(out, cg);

  out << "finish" << endl;
}

int main(int argc, char **argv) {
  DLPTopolApp app;
  return app.Exec(argc, argv);
}
