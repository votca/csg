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

#include <boost/format.hpp>
#include <stdio.h>
#include <string>
#include <votca/csg/pdbwriter.h>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;
void PDBWriter::Open(string file, bool bAppend) {
  if (bAppend) {
    _out.open(file, std::ios_base::app);
  } else {
    _out.open(file);
  }
}

void PDBWriter::WriteHeader(std::string header) {
  if (header.size() < 10 || header.substr(0, 10) != "HEADER    ") {
    _out << "HEADER    ";
  }
  _out << header;
  if (header.back() != '\n') _out << "\n";
}

void PDBWriter::Close() { _out.close(); }

/*<<<<<<< HEAD
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
=======*/
void PDBWriter::Write(CSG_Topology *conf) {

  _out << boost::format("MODEL     %1$4d\n") % (conf->getStep() + 1)
       << std::flush;
  ;
  WriteContainer<CSG_Topology>(conf, *conf);
  _out << "ENDMDL" << std::endl;
}
//>>>>>>> master
/*
void PDBWriter::writeSymmetry(Bead *bead) {
  if (bead->getSymmetry() > 1) {
    tools::vec r = 10 * bead->getPos();
    boost::format beadfrmt(
        "HETATM%1$5d %2$4s %3$3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f\n");
    tools::vec ru = 0.1 * bead->getU() + r;

<<<<<<< HEAD
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
=======*/
/*   _out << beadfrmt % (bead->getId() + 1) % 100000  // atom serial number
               % bead->getName()                    // atom name
               % "REU"                              // residue name
               % " "                                // chain identifier 1 char
               % (bead->getResnr() + 1)             // residue sequence number
               % ru.x() % ru.y() % ru.z();          // we skip the charge

   if (bead->getSymmetry() > 2) {
     tools::vec rv = 0.1 * bead->getV() + r;
     _out << beadfrmt % (bead->getId() + 1) % 100000  // atom serial number
                 % bead->getName()                    // atom name
                 % "REV"                              // residue name
                 % " "                        // chain identifier 1 char
                 % (bead->getResnr() + 1)     // residue sequence number
                 % rv.x() % rv.y() % rv.z();  // we skip the charge*/
//>>>>>>> master
/*    }
  }
  return;
}*/

/*std::string PDBWriter::getRes(Topology &conf, Bead *bead) {
  if (conf.getResidue(bead->getResnr())) {
    return conf.getResidue(bead->getResnr())->getName();
  } else {
    return "";
  }
}*/
}  // namespace csg
}  // namespace votca
