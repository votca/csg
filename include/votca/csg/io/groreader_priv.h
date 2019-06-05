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
#ifndef VOTCA_CSG_GROREADER_PRIV_H
#define VOTCA_CSG_GROREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double GROReader<Topology_T>::formatDistance_(const double &distance) {
  return converter_.convert(distance_unit, Topology_T::units::distance_unit) *
         distance;
}

template <class Topology_T>
double GROReader<Topology_T>::formatVelocity_(const double &velocity) {
  return converter_.convert(velocity_unit, Topology_T::units::velocity_unit) *
         velocity;
}

template <class Topology_T>
int GROReader<Topology_T>::formatId_(const int &id) {
  return id - 1;
}

template <class Topology_T>
std::string GROReader<Topology_T>::formatElement_(
    const std::string &atom_name) {
  std::string element = tools::topology_constants::unassigned_element;
  std::string atom_all_caps = boost::to_upper_copy<std::string>(atom_name);
  if (elements_.isEleFull(atom_all_caps)) {
    element = elements_.getEleShort(atom_all_caps);
  } else if (elements_.isEleShort(atom_name)) {
    element = atom_name;
  }
  return element;
}

template <class Topology_T>
double GROReader<Topology_T>::formatMass_(const std::string &element) {
  double atom_weight = 1.0;
  if (element.compare(tools::topology_constants::unassigned_element) == 0) {
    atom_weight = elements_.getMass(element);
  }
  return atom_weight;
}

template <class Topology_T>
bool GROReader<Topology_T>::ReadTopology(const std::string &file,
                                         boost::any top_any) {
  _topology = true;

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gro reader read topology, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  top.Cleanup();

  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open topology file: " + file);

  NextFrame(top_any);

  _fl.close();

  return true;
}

template <class Topology_T>
bool GROReader<Topology_T>::Open(const std::string &file) {
  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  return true;
}

template <class Topology_T>
void GROReader<Topology_T>::Close() {
  _fl.close();
}

template <class Topology_T>
bool GROReader<Topology_T>::FirstFrame(boost::any top_any) {
  _topology = false;
  NextFrame(top_any);
  return true;
}

template <class Topology_T>
bool GROReader<Topology_T>::NextFrame(boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gro reader next frame, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  std::string tmp;
  getline(_fl, tmp);  // title
  if (_fl.eof()) {
    return !_fl.eof();
  }
  getline(_fl, tmp);  // number atoms
  int natoms = atoi(tmp.c_str());
  if (!_topology && static_cast<size_t>(natoms) != top.BeadCount())
    throw std::runtime_error(
        "number of beads in topology and trajectory differ");

  for (int i = 0; i < natoms; i++) {
    std::string line;
    getline(_fl, line);
    std::string resNum, resName, atName, atNum, x, y, z;
    try {
      resNum = std::string(line, 0, 5);   // %5i
      resName = std::string(line, 5, 5);  //%5s
      atName = std::string(line, 10, 5);  // %5s
      atNum = std::string(line, 15, 5);   // %5i not needed
      x = std::string(line, 20, 8);       // %8.3f
      y = std::string(line, 28, 8);       // %8.3f
      z = std::string(line, 36, 8);       // %8.3f
    } catch (std::out_of_range &err) {
      throw std::runtime_error("Misformated gro file");
    }
    boost::algorithm::trim(atName);
    boost::algorithm::trim(atNum);
    boost::algorithm::trim(resName);
    boost::algorithm::trim(resNum);
    boost::algorithm::trim(x);
    boost::algorithm::trim(y);
    boost::algorithm::trim(z);
    std::string vx, vy, vz;
    bool hasVel = true;
    try {
      vx = std::string(line, 44, 8);  // %8.4f
      vy = std::string(line, 52, 8);  // %8.4f
      vz = std::string(line, 60, 8);  // %8.4f
    } catch (std::out_of_range &err) {
      hasVel = false;
    }

    if (_topology) {
      int residue_number = boost::lexical_cast<int>(resNum);
      if (residue_number < 1)
        throw std::runtime_error(
            "Misformated gro file, residue_number has to be > 0");

      // this is not correct, but still better than no type at all!

      // res -1 as internal number starts with 0
      tools::byte_t symmetry = 1;

      double atom_charge = 0.0;
      int atom_number = boost::lexical_cast<int>(atNum);
      std::string element = formatElement_(atName);
      double atom_weight = formatMass_(element);
      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::CSG_Mass, atom_weight);
      params.set(tools::StructureParameter::CSG_Charge, atom_charge);
      params.set(tools::StructureParameter::Element, element);
      params.set(tools::StructureParameter::BeadId, atom_number);
      params.set(tools::StructureParameter::BeadType, atName);
      params.set(tools::StructureParameter::ResidueId,
                 formatId_(residue_number));
      params.set(tools::StructureParameter::ResidueType, resName);
      params.set(tools::StructureParameter::MoleculeId,
                 tools::topology_constants::unassigned_molecule_id);
      top.CreateBead(params);
    }
    typename Topology_T::bead_t &b = top.getBead(i);

    b.setPos(Eigen::Vector3d(formatDistance_(stod(x)), formatDistance_(stod(y)),
                             formatDistance_(stod(z))));
    if (hasVel) {
      boost::algorithm::trim(vx);
      boost::algorithm::trim(vy);
      boost::algorithm::trim(vz);
      b.setVel(Eigen::Vector3d(formatVelocity_(stod(vx)),
                               formatVelocity_(stod(vy)),
                               formatVelocity_(stod(vz))));
    }
  }

  getline(_fl, tmp);  // read box line
  if (_fl.eof())
    throw std::runtime_error(
        "unexpected end of file in poly file, when boxline");
  tools::Tokenizer tok(tmp, " ");
  std::vector<double> fields;
  tok.ConvertToVector<double>(fields);
  Eigen::Matrix3d box;
  if (fields.size() == 3) {
    box = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; i++) box(i, i) = formatDistance_(fields[i]);
  } else if (fields.size() == 9) {
    box(0, 0) = formatDistance_(fields[0]);
    box(1, 1) = formatDistance_(fields[1]);
    box(2, 2) = formatDistance_(fields[2]);
    box(1, 0) = formatDistance_(fields[3]);
    box(2, 0) = formatDistance_(fields[4]);
    box(0, 1) = formatDistance_(fields[5]);
    box(2, 1) = formatDistance_(fields[6]);
    box(0, 2) = formatDistance_(fields[7]);
    box(1, 2) = formatDistance_(fields[8]);
  } else {
    throw std::runtime_error("Error while reading box (last) line");
  }
  top.setBox(box);

  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GROREADER_PRIV_H
