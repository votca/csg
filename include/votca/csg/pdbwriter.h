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

#ifndef VOTCA_CSG_PDBWRITER_H
#define VOTCA_CSG_PDBWRITER_H

#include <stdio.h>
#include <votca/csg/csgtopology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

/**
 * @brief
 *
 * Note that by default the coordinates stored in the pdb file are in units
 * of Angstroms while the following repositories have default units of:
 *
 * csg: nm
 * xtp: bohr
 *
 * Further note that pdb files store the ids, serical numbers and indices
 * starting at 1, votca uses ids, indices and serial numbers starting at 0.
 */
class PDBWriter : public TrajectoryWriter {
 public:
  void Open(std::string file, bool bAppend = false);
  void Close();

  void RegisteredAt(
      tools::ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(CSG_Topology *conf);

  template <class T>
  void WriteContainer(CSG_Topology *conf, T &container);

  void WriteHeader(std::string header);

 private:
  template <class Atom>
  std::string getType(Atom &atom) {
    std::string atomtype = atom.getType();
    if (atomtype.size() > 4) {
      atomtype = atomtype.substr(0, 4);
    }
    return atomtype;
  }

  template <class Atom>
  std::string getElement(Atom &atom) {
    if (atom.getElement() == tools::topology_constants::unassigned_element) {
      return "";
    }
    return atom.getElement();
  }

  // std::string getType(Bead *bead) { return bead->getType(); }

  /**
   * @brief Get the residue type
   *
   * In the case that the residue type is unknown it will follow the pdb
   * convention and assign it a value of "UNK", see the url for more info:
   *
   * www.wwpdb.org/documentation/file-format-content/format33/sect4.html
   *
   * @tparam Atom
   * @param atom
   *
   * @return
   */
  template <class Atom>
  std::string getResidueType(Atom &atom) {
    std::string restype = atom.getResidueType();
    if (restype == tools::topology_constants::unassigned_residue_type) {
      restype = "UNK";
    } else if (restype.size() > 3) {
      restype = restype.substr(0, 3);
    }
    return restype;
  }
  // std::string getResidueType(Bead *bead){
  //  return bead->getResidueType();
  //}

  template <class Atom>
  int getId(Atom &atom) {
    return atom.getId() + 1;
  }
  // int getId(Bead *bead) { return bead->getId(); }

  template <class T, class Atom>
  int getResId(Atom &atom) {
    return atom.getResidueId() + 1;
  }
  // int getResId(CSG_Topology &conf, Bead *bead) { return bead->getResnr() + 1;
  // }

  /*template <class Atom>
  void writeSymmetry(Atom &atom) {};
  void writeSymmetry(Bead & bead);*/

  template <class Atom>
  Eigen::Vector3d getPos(Atom &atom) {
    return atom.getPos() * tools::conv::bohr2ang;
  }

  Eigen::Vector3d getPos(Bead &bead) {
    return bead.Pos().toEigen() * tools::conv::nm2ang;
  }

  template <class T>
  std::vector<Bead *> getIterable(CSG_Topology &top, T &container) {
    std::vector<Bead *> beads;
    std::vector<int> bead_ids = container.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }

  std::vector<Bead *> getIterable(CSG_Topology &top) {
    std::vector<Bead *> beads;
    std::vector<int> bead_ids = top.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }

  std::ofstream _out;
};

template <class T>
inline void PDBWriter::WriteContainer(CSG_Topology *conf, T &container) {
  boost::format atomfrmt(
      "ATOM  %1$5d %2$-4s %3$-3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f           "
      "           %9$+2s\n");

  for (auto &atom : getIterable(conf, container)) {
    int atomid = getId(*atom);
    std::string resname = getResidueType(*atom);
    int residueid = getResId(*atom);
    std::string atomtype = getType(*atom);
    std::string element = getElement(*atom);
    Eigen::Vector3d r = getPos(*atom);

    _out << atomfrmt % (atomid % 100000)    // atom serial number
                % atomtype % resname % " "  // chain identifier 1 char
                % residueid                 // residue sequence number
                % r.x() % r.y() % r.z() % element;

    // we skip the charge
    // writeSymmetry(*atom);
  }
  _out << std::flush;
}

// Super Hacky does not follow the pdb file convention, essentially is breaking
// up coarse grained beads into arbitraty pieces, if we are going to allow
// this it should be exact and break the pieces back into their constituent
// atoms using the cgatomconverter class. Or you should simply use a differnt
// file format
/*
  void PDBWriter::writeSymmetry(Bead &bead) {
    if (bead.getSymmetry() > 1) {
      tools::vec r = tools::conv::nm2ang * bead.getPos();
      boost::format beadfrmt(
          "HETATM%1$5d %2$4s %3$3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f
  %9$+2s\n"); tools::vec ru = 0.1 * bead.getU() + r;

      _out << beadfrmt % getId(bead) % 100000  // bead serial number
        % bead.getType()                    // bead name
        % bead.getResidueType()                              // residue name
        % " "                               // chain identifier 1 char
        % getResId(bead)             // residue sequence number
        % ru.x() % ru.y() % ru.z()          // we skip the charge
        % getElement(bead);

      if (bead.getSymmetry() > 2) {
        tools::vec rv = 0.1 * bead.getV() + r;
        _out << beadfrmt % getId(bead) % 100000  // bead serial number
          % bead.getType()                    // bead name
          % bead.getResidueType()                  // residue name
          % " "                        // chain identifier 1 char
          % getResId(bead)     // residue sequence number
          % rv.x() % rv.y() % rv.z()  // we skip the charge
          % getElement(bead);
      }
    }
    return;
  }*/
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_PDBWRITER_H
