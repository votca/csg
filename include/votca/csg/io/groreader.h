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

#ifndef VOTCA_CSG_GROREADER_H
#define VOTCA_CSG_GROREADER_H

#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <votca/tools/elements.h>
#include <votca/tools/getline.h>
#include <votca/tools/structureparameters.h>

namespace votca {
namespace csg {

/**
    \brief reader for gro files

    This class provides the TrajectoryReader + Topology reader interface
    for gro files

*/
template <class Topology_T>
class GROReader : public TrajectoryReader, public TopologyReader {
 public:
  GROReader() {}
  ~GROReader() {}

  /// open a topology file
  bool ReadTopology(std::string file, boost::any top);

  /// open a trejectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  void Close();

 private:
  std::ifstream _fl;
  bool _topology;
};

template <class Topology_T>
bool GROReader<Topology_T>::ReadTopology(std::string file, boost::any top_any) {
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

  tools::Elements elements;
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

    typename Topology_T::bead_t *b;
    if (_topology) {
      int residue_number = boost::lexical_cast<int>(resNum);
      if (residue_number < 1)
        throw std::runtime_error(
            "Misformated gro file, residue_number has to be > 0");

      // this is not correct, but still better than no type at all!

      // res -1 as internal number starts with 0
      tools::byte_t symmetry = 1;
      std::string element = tools::topology_constants::unassigned_element;
      std::string atom_all_caps = boost::to_upper_copy<std::string>(atName);
      double atom_weight = 1.0;
      double atom_charge = 0.0;
      int atom_number = boost::lexical_cast<int>(atNum);
      if (elements.isEleFull(atom_all_caps)) {
        element = elements.getEleShort(atom_all_caps);
        atom_weight = elements.getMass(element);
      } else if (elements.isEleShort(atName)) {
        element = atName;
        atom_weight = elements.getMass(element);
      }
      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::Mass, atom_weight);
      params.set(tools::StructureParameter::Charge, atom_charge);
      params.set(tools::StructureParameter::Element, element);
      params.set(tools::StructureParameter::BeadId, atom_number);
      params.set(tools::StructureParameter::BeadType, atName);
      params.set(tools::StructureParameter::ResidueId, residue_number - 1);
      params.set(tools::StructureParameter::ResidueType, resName);
      params.set(tools::StructureParameter::MoleculeId,
                 tools::topology_constants::unassigned_molecule_id);
      b = top.CreateBead(params);
    } else {
      b = top.getBead(i);
    }

    b->setPos(Eigen::Vector3d(stod(x), stod(y), stod(z)));
    if (hasVel) {
      boost::algorithm::trim(vx);
      boost::algorithm::trim(vy);
      boost::algorithm::trim(vz);
      b->setVel(Eigen::Vector3d(stod(vx), stod(vy), stod(vz)));
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
    for (int i = 0; i < 3; i++) box(i, i) = fields[i];
  } else if (fields.size() == 9) {
    box(0, 0) = fields[0];
    box(1, 1) = fields[1];
    box(2, 2) = fields[2];
    box(1, 0) = fields[3];
    box(2, 0) = fields[4];
    box(0, 1) = fields[5];
    box(2, 1) = fields[6];
    box(0, 2) = fields[7];
    box(1, 2) = fields[8];
  } else {
    throw std::runtime_error("Error while reading box (last) line");
  }
  top.setBox(box);

  if (_topology) {
    std::cout
        << "WARNING: topology created from .gro file, masses and charges are "
           "wrong!\n";
  }

  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GROREADER_H
