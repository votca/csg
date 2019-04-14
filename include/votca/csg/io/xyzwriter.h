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

#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <stdio.h>
#include <votca/tools/constants.h>
namespace votca {
namespace csg {

template <class Bead_T, class Molecule_T, class Topology_T>
class XYZWriter : public TrajectoryWriter {
 public:
  void RegisteredAt(
      tools::ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(boost::any conf);

  template <class T>
  void Write(Topology_T &top, T &container, std::string header);

 private:
  template <class T>
  std::string getType(T &atom) {
    return atom.getType();
  }

  std::string getType(Bead_T *bead) { return bead->getType(); }

  template <class T>
  Eigen::Vector3d getPos(T &atom) {
    return atom.getPos() * tools::conv::bohr2ang;
  }

  Eigen::Vector3d getPos(Bead_T *bead) {
    return bead->Pos() * tools::conv::nm2ang;
  }

  // The CSG Topology object is the only object that stores the beads and its
  // pointers, all other containers only store the bead ids
  template <typename T>
  std::vector<Bead_T *> getIterable(Topology_T &top, T &container) {
    std::vector<Bead_T *> beads;
    std::vector<int> bead_ids = container.getBeadIds();
    for (int &bead_id : bead_ids) {
      beads.push_back(top.getBead(bead_id));
    }
    return beads;
  }
};

template <class Bead_T, class Molecule_T, class Topology_T>
void XYZWriter<Bead_T, Molecule_T, Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz writer write, incorrect topology "
        "type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  std::string header = (boost::format("frame: %1$d time: %2$f\n") %
                        (top.getStep() + 1) % top.getTime())
                           .str();
  Write(top, top, header);
}

template <class Bead_T, class Molecule_T, class Topology_T>
template <typename T>
inline void XYZWriter<Bead_T, Molecule_T, Topology_T>::Write(
    Topology_T &top, T &container, std::string header) {

  std::vector<Bead_T *> atoms = getIterable(top, container);
  out_ << atoms.size() << "\n";
  out_ << header << "\n";

  boost::format fmter("%1$s%2$10.5f%3$10.5f%4$10.5f\n");

  for (Bead_T *atom : atoms) {
    Eigen::Vector3d r = getPos(atom);
    // truncate strings if necessary
    std::string atomtype = getType(atom);
    if (atomtype.size() > 3) {
      atomtype = atomtype.substr(0, 3);
    }
    while (atomtype.size() < 3) atomtype = " " + atomtype;

    out_ << fmter % atomtype % r.x() % r.y() % r.z();
  }
  out_ << std::flush;
}
}  // namespace csg
}  // namespace votca

#endif
