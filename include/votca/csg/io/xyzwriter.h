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

#pragma once
#ifndef VOTCA_CSG_XYZWRITER_H
#define VOTCA_CSG_XYZWRITER_H

#include "../csgtopology.h"
#include "../molecule.h"
#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <stdio.h>
#include <votca/tools/constants.h>
namespace votca {
namespace csg {

template <class Topology_T>
class XYZWriter : public TrajectoryWriter {
 public:
  void Write(boost::any conf);

  template <class T>
  void WriteContainer(T &container);

  void WriteHeader(std::string header, int number_atoms);

 private:
  template <class T>
  void WriteContainer_(T &container);
  /*
  template <class T>
  std::string getType(T &atom) {
    return atom.getType();
  }

  template <class T>
  Eigen::Vector3d getPos(T &atom) {
    return atom.getPos() * tools::conv::bohr2ang;
  }*/

  // The CSG Topology object is the only object that stores the beads and its
  // pointers, all other containers only store the bead ids
  /*  template <typename T>
    std::vector<typename Topology_T::bead_t *> getIterable(Topology_T &top,
                                                           T &container) {
      std::vector<typename Topology_T::bead_t *> beads;
      std::vector<int> bead_ids = container.getBeadIds();
      for (int &bead_id : bead_ids) {
        beads.push_back(top.getBead(bead_id));
      }
      return beads;
    }*/

  template <class T>
  int getAtomCount(T &container) {
    return static_cast<int>(container.size());
  }

  int getAtomCount(Topology_T &container) {
    return static_cast<int>(container.BeadCount());
  }

  template <class T>
  T &getIterable(T &container) {
    return container;
  }

  std::vector<typename Molecule::bead_t> getIterable(Molecule &container) {
    return container.getBeads();
  }

  template <typename T>
  T *ptr(T *obj) {
    return obj;
  }

  template <typename T>
  T *ptr(T &obj) {
    return &obj;
  }
};

template <class Topology_T>
void XYZWriter<Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz writer write, incorrect topology "
        "type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  std::string header = (boost::format("frame: %1$d time: %2$f\n") %
                        (top.getStep() + 1) % top.getTime())
                           .str();

  WriteHeader(header, getAtomCount(top));
  for (auto &container : top) {
    WriteContainer_(container);
  }
}

template <class Topology_T>
inline void XYZWriter<Topology_T>::WriteHeader(std::string header,
                                               int number_atoms) {
  out_ << number_atoms << "\n";
  out_ << header << "\n";
}

template <class Topology_T>
template <typename T>
inline void XYZWriter<Topology_T>::WriteContainer(T &container) {
  WriteHeader("", getAtomCount(container));
  WriteContainer_(container);
}

template <class Topology_T>
template <typename T>
inline void XYZWriter<Topology_T>::WriteContainer_(T &container) {

  boost::format fmter("%1$s%2$10.5f%3$10.5f%4$10.5f\n");
  std::vector<int> ids = container.getBeadIds();
  for (const int &id : ids) {
    //  for (auto & atom : getIterable(container)) {
    //    auto atom_ptr = ptr(atom);
    auto atom_ptr = container.getBead(id);
    tools::StructureParameters params = atom_ptr->getParameters();
    Eigen::Vector3d r = params.get<Eigen::Vector3d>(
        tools::StructureParameter::Position);  // getPos(atom);
    // truncate strings if necessary
    std::string atomtype = params.get<std::string>(
        tools::StructureParameter::BeadType);  // getType(atom);
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
