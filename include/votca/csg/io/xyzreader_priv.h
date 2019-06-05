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
#ifndef VOTCA_CSG_XYZREADER_PRIV_H
#define VOTCA_CSG_XYZREADER_PRIV_H

namespace votca {
namespace csg {

namespace TopConst = tools::topology_constants;

template <class Topology_T>
template <bool topology>
inline void XYZReader<Topology_T>::AddAtom(Topology_T &container,
                                           tools::StructureParameters &params) {

  Eigen::Vector3d posnm =
      params.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);
  if (topology) {
    container.CreateMolecule(params).AddBead(container.CreateBead(params));
  }
  typename Topology_T::bead_t &b =
      container.getBead(params.get<int>(tools::StructureParameter::BeadId));

  b.setPos(posnm);
}

template <class Topology_T>
template <class T>
inline void XYZReader<Topology_T>::ReadFile(T &container) {
  if (!ReadFrame<true, T>(container)) {
    throw std::runtime_error("Reading xyz file '" + _file + "' failed");
  }
}

template <class Topology_T>
bool XYZReader<Topology_T>::FirstFrame(boost::any top) {
  return NextFrame(top);
}

template <class Topology_T>
bool XYZReader<Topology_T>::NextFrame(boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz reader next frame, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);

  bool success = ReadFrame<false>(top);
  return success;
}

template <class Topology_T>
bool XYZReader<Topology_T>::ReadTopology(const std::string &file,
                                         boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz reader read topology, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  top.Cleanup();
  top.setStep(0);

  _file = file;
  _fl.open(file);
  if (!_fl.is_open()) {
    throw std::ios_base::failure("Error on open topology file: " + file);
  }

  ReadFrame<true, Topology_T>(top);

  _fl.close();

  return true;
}

template <class Topology_T>
inline void XYZReader<Topology_T>::formatElement(std::string &element) {
  std::string name_upper_case = boost::to_upper_copy<std::string>(element);
  if (elements_.isEleFull(name_upper_case)) {
    element = elements_.getEleShort(name_upper_case);
  } else if (!(elements_.isEleShort(element))) {
    element = TopConst::unassigned_element;
  }
}

template <class Topology_T>
inline void XYZReader<Topology_T>::formatDistance(double &dist) {
  dist *=
      converter_.convert(this->distance_unit, Topology_T::units::distance_unit);
}

template <class Topology_T>
template <bool topology, class T>
inline bool XYZReader<Topology_T>::ReadFrame(T &container) {

  std::string line;
  getline(_fl, line);
  ++_line;
  if (!_fl.eof()) {
    // read the number of atoms
    tools::Tokenizer tok1(line, " \t");
    std::vector<std::string> line1;
    tok1.ToVector(line1);
    if (line1.size() != 1) {
      throw std::runtime_error(
          "First line of xyz file should contain number "
          "of atoms/beads, nothing else.");
    }
    size_t natoms = boost::lexical_cast<size_t>(line1[0]);
    if (!topology && natoms != container.BeadCount()) {
      throw std::runtime_error(
          "number of beads in topology and trajectory differ");
    }
    // the title line
    getline(_fl, line);
    ++_line;
    // read atoms

    for (int bead_id = 0; bead_id < natoms; ++bead_id) {
      getline(_fl, line);
      ++_line;

      if (_fl.eof()) {
        throw std::runtime_error("unexpected end of file in xyz file");
      }
      std::vector<std::string> fields;
      tools::Tokenizer tok(line, " ");
      tok.ToVector(fields);
      if (fields.size() != 4) {
        throw std::runtime_error("invalide line " +
                                 boost::lexical_cast<std::string>(_line) +
                                 " in xyz file\n" + line);
      }

      std::string element = fields[0];
      std::string atom_type = fields[0];
      double x_pos = boost::lexical_cast<double>(fields[1]);
      double y_pos = boost::lexical_cast<double>(fields[2]);
      double z_pos = boost::lexical_cast<double>(fields[3]);

      formatElement(element);
      formatDistance(x_pos);
      formatDistance(y_pos);
      formatDistance(z_pos);

      Eigen::Vector3d pos = Eigen::Vector3d(x_pos, y_pos, z_pos);
      tools::byte_t symmetry = 1;
      double mass = 0.0;
      double charge = 0.0;
      // Because xyz files do not contain any bond info each atom is assumed
      // to be in its own molecule
      int molecule_id = bead_id;
      int residue_id = TopConst::unassigned_residue_id;
      std::string residue_type = TopConst::unassigned_residue_type;
      std::string molecule_type = TopConst::unassigned_molecule_type;

      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::CSG_Mass, mass);
      params.set(tools::StructureParameter::CSG_Charge, charge);
      params.set(tools::StructureParameter::Element, element);
      params.set(tools::StructureParameter::MoleculeId, molecule_id);
      params.set(tools::StructureParameter::MoleculeType, molecule_type);
      params.set(tools::StructureParameter::ResidueId, residue_id);
      params.set(tools::StructureParameter::ResidueType, residue_type);
      params.set(tools::StructureParameter::BeadId, bead_id);
      params.set(tools::StructureParameter::BeadType, atom_type);
      params.set(tools::StructureParameter::CSG_Position, pos);

      AddAtom<topology>(container, params);
    }
  }
  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca

#endif
