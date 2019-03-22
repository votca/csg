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

#ifndef __VOTCA_CSG_XYZWRITER_H
#define __VOTCA_CSG_XYZWRITER_H

#include <stdio.h>
#include <votca/csg/csgtopology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

class XYZWriter : public TrajectoryWriter {
 public:
  void Open(std::string file, bool bAppend = false);
  void Close();

  void RegisteredAt(
      tools::ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(CSG_Topology *conf);

  template <class T>
  void Write(CSG_Topology &top, T &container, std::string header);

 private:
  //  template <class T>
  /*  int getSize(T &container) {
      return getIterable(container).size();
    }*/

  template <class Atom>
  std::string getType(Atom &atom) {
    return atom.getType();
  }

  std::string getType(Bead *bead) { return bead->getType(); }

  template <class Atom>
  Eigen::Vector3d getPos(Atom &atom) {
    return atom.getPos() * tools::conv::bohr2ang;
  }

  Eigen::Vector3d getPos(Bead *bead) {
    return bead->Pos().toEigen() * tools::conv::nm2ang;
  }

  // The CSG Topology object is the only object that stores the beads and its
  // pointers, all other containers only store the bead ids
  template <class T>
  std::vector<Bead *> getIterable(CSG_Topology &top, T &container) {
    std::vector<Bead *> beads;
    std::vector<int> bead_ids = container.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }

  std::ofstream _out;
};

template <class T>
inline void XYZWriter::Write(CSG_Topology &top, T &container,
                             std::string header) {

  std::vector<Bead *> atoms = getIterable(top, container);
  _out << atoms.size() << "\n";
  _out << header << "\n";

  boost::format fmter("%1$s%2$10.5f%3$10.5f%4$10.5f\n");

  for (auto &atom : atoms) {
    Eigen::Vector3d r = getPos(atom);
    // truncate strings if necessary
    std::string atomtype = getType(atom);
    if (atomtype.size() > 3) {
      atomtype = atomtype.substr(0, 3);
    }
    while (atomtype.size() < 3) atomtype = " " + atomtype;

    _out << fmter % atomtype % r.x() % r.y() % r.z();
  }
  _out << std::flush;
}
}  // namespace csg
}  // namespace votca

#endif
