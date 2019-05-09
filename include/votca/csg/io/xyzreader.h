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
#ifndef __VOTCA_CSG_XYZREADER_H
#define __VOTCA_CSG_XYZREADER_H

#include "../csgtopology.h"
#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

#include "../molecule.h"
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/structureparameters.h>
namespace votca {
namespace csg {

/**
    \brief class for reading xyz files

    This class provides the TrajectoryReader + Topology reader interface
    for xyz files

*/
template <class Topology_T>
class XYZReader : public TrajectoryReader, public TopologyReader {
 public:
  XYZReader() {}
  ~XYZReader() {}

  /// open a topology file
  template <bool topology>
  bool ReadTopology(std::string file, boost::any top);

  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  template <class T>
  void ReadFile(T &container) {
    if (!ReadFrame<true, T>(container)) {
      throw std::runtime_error("Reading xyz file '" + _file + "' failed");
    }
  }

 private:
  template <class T>
  size_t getContainerSize(T &container) {
    return container.size();
  }

  size_t getContainerSize(Topology_T &top) { return top.BeadCount(); }

  template <bool topology, class T>
  void AddAtom(T &container, tools::StructureParameters &params) {
    // the typedef returns the type of the objects the container holds
    typedef
        typename std::iterator_traits<decltype(container.begin())>::value_type
            atom;
    Eigen::Vector3d pos2 =
        params.get<Eigen::Vector3d>(tools::StructureParameter::Position) *
        tools::conv::ang2bohr;
    params.set(tools::StructureParameter::Position, pos2);
    // container.push_back(atom(id, name, pos2));
    container.push_back(atom(params));
  }

  template <bool topology>
  void AddAtom(Molecule *molecule, tools::StructureParameters &params) {}

  template <bool topology>
  void AddAtom(Topology_T &container, tools::StructureParameters &params) {

    typename Topology_T::bead_t *b;
    Eigen::Vector3d posnm =
        params.get<Eigen::Vector3d>(tools::StructureParameter::Position) *
        tools::conv::ang2nm;
    if (topology) {
      b = &(container.CreateBead(params));
    } else {
      b = container.getBead(params.get<int>(tools::StructureParameter::BeadId));
    }
    b->setPos(posnm);
  }

  /// For adding to xtp topology assume adding all the atoms to a single segment
  /*template<bool topology>
    void AddAtom(Topology_T & top, tools::StructureParameters params){
      if(top.Segments.size()==0){
        top.AddSegment(tools::topology_constants::unassigned_segment_type);
      }else{
        top.begin()->push_back(top.bead_t(params));
      }
    }*/
  /*  template <bool topology, class T>
    void AddAtom(T &container, std::string bead_type, int bead_id,
                 std::string element, const Eigen::Vector3d &pos) {

      typename Topology_T::bead_t *b;
      Eigen::Vector3d posnm = pos * tools::conv::ang2nm;
      if (topology) {

        tools::byte_t symmetry = 1;
        double mass = 0.0;
        double charge = 0.0;
        std::string element = tools::topology_constants::unassigned_element;
        int molecule_id = tools::topology_constants::unassigned_molecule_id;
        int residue_id = tools::topology_constants::unassigned_residue_id;
        std::string residue_type =
            tools::topology_constants::unassigned_residue_type;

        tools::StructureParameters params;
        params.set(tools::StructureParameter::Symmetry, symmetry);
        params.set(tools::StructureParameter::Mass, mass);
        params.set(tools::StructureParameter::Charge, charge);
        params.set(tools::StructureParameter::Element, element);
        params.set(tools::StructureParameter::MoleculeId, molecule_id);
        params.set(tools::StructureParameter::ResidueId, residue_id);
        params.set(tools::StructureParameter::ResidueType, residue_type);
        params.set(tools::StructureParameter::BeadId, bead_id);
        params.set(tools::StructureParameter::BeadType, bead_type);

        b = container.CreateBead(params);

      } else {
        b = container.getBead(bead_id);
      }
      b->setPos(posnm);
    }*/

  template <bool topology, class T>
  bool ReadFrame(T &container);

  std::ifstream _fl;
  std::string _file;
  int _line;
};

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
template <bool topology>
bool XYZReader<Topology_T>::ReadTopology(std::string file, boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using xyz reader read topology, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  top.Cleanup();

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
    if (!topology && natoms != getContainerSize(container)) {
      throw std::runtime_error(
          "number of beads in topology and trajectory differ");
    }
    // the title line
    getline(_fl, line);
    ++_line;
    tools::Elements elements;
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

      std::string name_upper_case =
          boost::to_upper_copy<std::string>(fields[0]);
      std::string element = tools::topology_constants::unassigned_element;
      if (elements.isEleFull(name_upper_case)) {
        element = elements.getEleShort(name_upper_case);
      } else if (elements.isEleShort(fields[0])) {
        element = fields[0];
      }

      Eigen::Vector3d pos =
          Eigen::Vector3d(boost::lexical_cast<double>(fields[1]),
                          boost::lexical_cast<double>(fields[2]),
                          boost::lexical_cast<double>(fields[3]));

      tools::byte_t symmetry = 1;
      double mass = 0.0;
      double charge = 0.0;
      int molecule_id = tools::topology_constants::unassigned_molecule_id;
      int residue_id = tools::topology_constants::unassigned_residue_id;
      std::string residue_type =
          tools::topology_constants::unassigned_residue_type;

      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::Mass, mass);
      params.set(tools::StructureParameter::Charge, charge);
      params.set(tools::StructureParameter::Element, element);
      params.set(tools::StructureParameter::MoleculeId, molecule_id);
      params.set(tools::StructureParameter::ResidueId, residue_id);
      params.set(tools::StructureParameter::ResidueType, residue_type);
      params.set(tools::StructureParameter::BeadId, bead_id);
      params.set(tools::StructureParameter::BeadType, fields[0]);
      params.set(tools::StructureParameter::Position, pos);

      AddAtom<topology>(container, params);
    }
  }
  return !_fl.eof();
}
}  // namespace csg
}  // namespace votca

#endif
