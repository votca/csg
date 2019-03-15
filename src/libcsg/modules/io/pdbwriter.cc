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

#include "pdbwriter.h"
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;
void PDBWriter::Open(string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void PDBWriter::Close() { fclose(_out); }

void PDBWriter::Write(CSG_Topology *conf) {
  fprintf(_out, "MODEL     %4d\n", conf->getStep());
  vector<int> bead_ids = conf->getBeadIds();
  for (const int &bead_id : bead_ids) {
    Bead *bi = conf->getBead(bead_id);
    vec r = bi->getPos();
    // truncate strings if necessary
    string residue_type = bi->getResidueType();
    string atom_type = bi->getType();
    if (residue_type.size() > 3) {
      residue_type = residue_type.substr(0, 3);
    }
    if (atom_type.size() > 4) {
      atom_type = atom_type.substr(0, 4);
    }
    string element = "  ";
    if (bi->getElement().size() > 3) {
      throw runtime_error("Unrecognized element type when writing pdb file");
    }
    element = bi->getElement();
    fprintf(_out, "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%24.2s \n",
            (bi->getId() + 1) % 100000,          // atom serial number
            atom_type.c_str(),                   // atom type
            residue_type.c_str(),                // residue type
            " ",                                 // chain identifier 1 char
            bi->getResidueId() + 1,              // residue sequence number
            10 * r.x(), 10 * r.y(), 10 * r.z(),  // nm -> Angs
            element.c_str());                    // Element if it is known
                                                 // we skip the charge

    if (bi->getSymmetry() >= 2) {
      vec ru = 0.1 * bi->getU() + r;

      fprintf(_out, "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f\n",
              bi->getId() + 1,          // atom serial number
              bi->getType().c_str(),    // atom type
              "REU",                    // residue type
              " ",                      // chain identifier 1 char
              bi->getResidueId() + 1,   // residue sequence number
              ru.x(), ru.y(), ru.z());  // we skip the charge
    }
    if (bi->getSymmetry() >= 3) {
      vec rv = 0.1 * bi->getV() + r;
      fprintf(_out, "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f\n",
              bi->getId() + 1,          // atom serial number
              bi->getType().c_str(),    // atom type
              "REV",                    // residue type
              " ",                      // chain identifier 1 char
              bi->getResidueId() + 1,   // residue sequence number
              rv.x(), rv.y(), rv.z());  // we skip the charge
    }
  }
  fprintf(_out, "ENDMDL\n");
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
