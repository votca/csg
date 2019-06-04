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
#ifndef VOTCA_CSG_XYZWRITER_PRIV_H
#define VOTCA_CSG_XYZWRITER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
template <class T>
inline int XYZWriter<Topology_T>::getAtomCount(T &container) {
  return static_cast<int>(container.size());
}

template <class Topology_T>
inline int XYZWriter<Topology_T>::getAtomCount(Topology_T &container) {
  return static_cast<int>(container.BeadCount());
}

template <class Topology_T>
template <class T>
inline T &XYZWriter<Topology_T>::getIterable(T &container) {
  return container;
}

template <class Topology_T>
inline std::vector<typename Molecule::bead_t>
    XYZWriter<Topology_T>::getIterable(Molecule &container) {
  return container.getBeads();
}

template <class Topology_T>
template <typename T>
inline T *XYZWriter<Topology_T>::ptr(T *obj) {
  return obj;
}

template <class Topology_T>
template <typename T>
inline T *XYZWriter<Topology_T>::ptr(T &obj) {
  return &obj;
}

template <class Topology_T>
inline void XYZWriter<Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz writer write, incorrect topology "
        "type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  std::string header = (boost::format("frame: %1$d time: %2$f") %
                        (top.getStep() + 1) % top.getTime())
                           .str();

  WriteHeader(header, getAtomCount(top));
  for (auto &container : top) {
    std::cout << "Writing container" << std::endl;
    WriteContainer_(container);
  }
}

template <class Topology_T>
inline void XYZWriter<Topology_T>::formatType(std::string &atomtype) {
  if (atomtype.size() > 3) {
    atomtype = atomtype.substr(0, 3);
  }
  while (atomtype.size() < 3) atomtype = " " + atomtype;
}

template <class Topology_T>
inline void XYZWriter<Topology_T>::formatPosition(Eigen::Vector3d &position) {
  position *=
      converter_.convert(Topology_T::units::distance_unit, this->distance_unit);
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
  std::cout << "Number of bead ids " << ids.size() << " in molecule "
            << container.getId() << std::endl;
  for (const int &id : ids) {
    std::cout << "getting bead id " << id << std::endl;
    auto &atom_ptr = container.getBead(id);
    tools::StructureParameters params = atom_ptr.getParameters();
    Eigen::Vector3d r =
        params.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);
    std::string atomtype =
        params.get<std::string>(tools::StructureParameter::BeadType);

    formatType(atomtype);
    formatPosition(r);

    out_ << fmter % atomtype % r.x() % r.y() % r.z();
  }
  out_ << std::flush;
}

}  // namespace csg
}  // namespace votca

#endif
