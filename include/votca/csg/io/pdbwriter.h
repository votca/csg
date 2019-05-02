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

#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <stdio.h>
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
template <typename Bead_T, class Molecule_T, class Topology_T>
class PDBWriter : public TrajectoryWriter {
 public:
  void Write(boost::any conf);

  template <class T>
  void WriteContainer(Topology_T *conf, T &container);

  void WriteHeader(std::string header);

  void WriteBox(const Eigen::Matrix3d &box);

 private:
  template <class T>
  std::string getType(T &atom) {
    std::string type = atom.getType();
    if (type.size() > 4) {
      type = type.substr(0, 4);
    }
    return type;
  }

  template <class T>
  std::string getElement(T &item) {
    if (item.getElement() == tools::topology_constants::unassigned_element) {
      return "";
    }
    return item.getElement();
  }

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
  template <class T>
  std::string getResidueType(T &item) {
    std::string restype = item.getResidueType();
    if (restype == tools::topology_constants::unassigned_residue_type) {
      restype = "UNK";
    } else if (restype.size() > 3) {
      restype = restype.substr(0, 3);
    }
    return restype;
  }

  template <class T>
  int getId(T &item) {
    return item.getId() + 1;
  }

  template <class T>
  int getResId(T &item) {
    return item.getResidueId() + 1;
  }

  template <class T>
  Eigen::Vector3d getPos(T &item) {
    return item.getPos() * tools::conv::bohr2ang;
  }

  Eigen::Vector3d getPos(Bead_T &bead) {
    return bead.Pos() * tools::conv::nm2ang;
  }

  template <class T>
  std::vector<Bead_T *> getIterable(Topology_T &top, T &container) {
    std::vector<Bead_T *> beads;
    std::vector<int> bead_ids = container.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }

  std::vector<Bead_T *> getIterable(Topology_T &top) {
    std::vector<Bead_T *> beads;
    std::vector<int> bead_ids = top.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }
};

template <class Bead_T, class Molecule_T, class Topology_T>
template <class T>
inline void PDBWriter<Bead_T, Molecule_T, Topology_T>::WriteContainer(
    Topology_T *conf, T &container) {

  if (out_.is_open()) {
    boost::format atomfrmt(
        "ATOM  %1$5d %2$-4s %3$-3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f         "
        "  "
        "           %9$+2s\n");

    std::vector<Bead_T *> atoms = getIterable(*conf, container);
    for (Bead_T *atom : atoms) {
      int atomid = getId(*atom);
      std::string resname = getResidueType(*atom);
      int residueid = getResId(*atom);
      std::string atomtype = getType(*atom);
      std::string element = getElement(*atom);
      Eigen::Vector3d r = getPos(*atom);

      out_ << atomfrmt % (atomid % 100000)    // atom serial number
                  % atomtype % resname % " "  // chain identifier 1 char
                  % residueid                 // residue sequence number
                  % r.x() % r.y() % r.z() % element;

      // we skip the charge
      // writeSymmetry(*atom);
    }
    out_ << std::flush;
  } else {
    throw std::runtime_error("Cannot write container to file it is not open.");
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void PDBWriter<Bead_T, Molecule_T, Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using pdb writer write, incorrect topology "
        "type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);
  if (out_.is_open()) {
    out_ << boost::format("MODEL     %1$4d\n") % (conf.getStep() + 1)
         << std::flush;
    WriteContainer<Topology_T>(&conf, conf);
    out_ << "ENDMDL" << std::endl;
  } else {
    throw std::runtime_error(
        "Cannot write topology to file, file is not open.");
  }
}
// Super Hacky does not follow the pdb file convention, essentially is breaking
// up coarse grained beads into arbitraty pieces, if we are going to allow
// this it should be exact and break the pieces back into their constituent
// atoms using the cgatomconverter class. Or you should simply use a differnt
// file format
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_PDBWRITER_H
