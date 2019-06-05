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
#ifndef VOTCA_CSG_PDBWRITER_PRIV_H
#define VOTCA_CSG_PDBWRITER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatType_(std::string &type) {
  if (type.size() > 4) {
    type = type.substr(0, 4);
  }
}

template <class Topology_T>

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatElement_(std::string &element) {
  if (element.compare(tools::topology_constants::unassigned_element) == 0) {
    element = "";
  }
}

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatResidueType_(std::string &restype) {
  if (restype.compare(tools::topology_constants::unassigned_residue_type) ==
      0) {
    restype = "UNK";
  } else if (restype.size() > 3) {
    restype = restype.substr(0, 3);
  }
}

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatId_(int &id) noexcept {
  ++id;
}

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatResId_(int &resId) noexcept {
  ++resId;
}

template <class Topology_T>
inline void PDBWriter<Topology_T>::formatPos_(Eigen::Vector3d &pos) {
  pos *=
      converter_.convert(Topology_T::units::distance_unit, this->distance_unit);
}

template <class Topology_T>
template <class T>
inline T &PDBWriter<Topology_T>::getIterable(T &container) {
  return container;
}

template <class Topology_T>
inline std::vector<typename Molecule::bead_t>
    PDBWriter<Topology_T>::getIterable(Molecule &container) {
  return container.getBeads();
}

template <class Topology_T>
template <typename T>
inline T *PDBWriter<Topology_T>::ptr(T *obj) {
  return obj;
}

template <class Topology_T>
template <typename T>
inline T *PDBWriter<Topology_T>::ptr(T &obj) {
  return &obj;
}

template <class Topology_T>
template <typename T>
inline void PDBWriter::WriteBox(const T &cont) {
  Eigen::Matrix3d box = cont.getBox();
  box *= converter_.convert(Topology_T::units::distance_unit,
                            cont::units::distance_unit);
  boost::format boxfrmt("CRYST1%1$9.3f%2$9.3f%3$9.3f%4$7.2f%5$7.2f%6$7.2f\n");
  double a = box.col(0).norm();
  double b = box.col(1).norm();
  double c = box.col(2).norm();
  double alpha =
      180 / tools::conv::Pi * std::acos(box.col(1).dot(box.col(2)) / (b * c));
  double beta =
      180 / tools::conv::Pi * std::acos(box.col(0).dot(box.col(2)) / (a * c));
  double gamma =
      180 / tools::conv::Pi * std::acos(box.col(0).dot(box.col(1)) / (a * b));
  _out << boxfrmt % a % b % c % alpha % beta % gamma;
}

template <class Topology_T>
template <class T>
inline void PDBWriter<Topology_T>::WriteContainer(T &container) {

  if (out_.is_open()) {
    boost::format atomfrmt(
        "ATOM  %1$5d %2$-4s %3$-3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f         "
        "  "
        "           %9$+2s\n");

    for (auto &atom : getIterable(container)) {
      auto atom_ptr = ptr(atom);
      tools::StructureParameters params = atom_ptr->getParameters();
      int atomid = params.get<int>(tools::StructureParameter::BeadId);
      std::string resname = tools::topology_constants::unassigned_residue_type;
      if (params.ParameterExist(tools::StructureParameter::ResidueType)) {
        resname =
            params.get<std::string>(tools::StructureParameter::ResidueType);
      }
      int residueid = params.get<int>(tools::StructureParameter::ResidueId);
      std::string atomtype =
          params.get<std::string>(tools::StructureParameter::BeadType);
      std::string element =
          params.get<std::string>(tools::StructureParameter::Element);
      Eigen::Vector3d r =
          params.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);

      formatType_(atomtype);
      formatResidueType_(resname);
      formatId_(atomid);
      formatResId_(residueid);
      formatPos_(r);
      formatElement_(element);

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

template <class Topology_T>
void PDBWriter<Topology_T>::WriteHeader(std::string header) {
  if (header.size() < 10 || header.substr(0, 10) != "HEADER    ") {
    out_ << "HEADER    ";
  }
  out_ << header;
  if (header.back() != '\n') out_ << "\n";
}

template <class Topology_T>
void PDBWriter<Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using pdb writer write, incorrect topology "
        "type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);
  if (out_.is_open()) {

    WriteHeader("Frame: " + std::to_string(conf.getStep()) + " Time: " + std::to_string(conf.getTime() );

    WriteBox(conf);
    out_ << boost::format("MODEL     %1$4d\n") % (conf.getStep() + 1)
         << std::flush;
    for (auto &container : conf) {
      WriteContainer(container);
    }
    out_ << "ENDMDL" << std::endl;
  } else {
    throw std::runtime_error(
        "Cannot write topology to file, file is not open.");
  }
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_PDBWRITER_PRIV_H
