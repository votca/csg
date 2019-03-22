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

#ifndef __VOTCA_CSG_XYZREADER_H
#define __VOTCA_CSG_XYZREADER_H

#include "csgtopology.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

namespace votca {
namespace csg {

/**
    \brief class for reading xyz files

    This class provides the TrajectoryReader + Topology reader interface
    for xyz files

*/
class XYZReader : public TrajectoryReader, public TopologyReader {
 public:
  XYZReader() {}
  ~XYZReader() {}

  /// open a topology file
  bool ReadTopology(std::string file, CSG_Topology &top);

  /// open a trajectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(CSG_Topology &top);
  /// read in the next frame
  bool NextFrame(CSG_Topology &top);

  template <class T>
  void ReadFile(T &container) {
    if (!ReadFrame<true, T>(container)) {
      throw std::runtime_error("Reading xyz file failed");
    }
  }

  void Close();

 private:
  template <class T>
  int getContainerSize(T &container) {
    return container.size();
  }

  int getContainerSize(CSG_Topology &container) {
    return container.BeadCount();
  }

  template <bool topology, class T>
  void AddAtom(T &container, std::string name, int id, const tools::vec &pos) {
    // the typedef returns the type of the objects the container holds
    typedef
        typename std::iterator_traits<decltype(container.begin())>::value_type
            atom;
    Eigen::Vector3d pos2 = pos.toEigen() * tools::conv::ang2bohr;
    container.push_back(atom(id, name, pos2));
  }

  template <bool topology, class T>

  void AddAtom(CSG_Topology &container, std::string bead_type, int bead_id,
               std::string element, const tools::vec &pos) {

    Bead *b;
    tools::vec posnm = pos * tools::conv::ang2nm;
    if (topology) {

      tools::byte_t symmetry = 1;
      std::string element = tools::topology_constants::unassigned_element;
      int molecule_id = tools::topology_constants::unassigned_molecule_id;
      int residue_id = tools::topology_constants::unassigned_residue_id;
      std::string residue_type =
          tools::topology_constants::unassigned_residue_type;
      double mass = 0.0;
      double charge = 0.0;

      b = container.CreateBead(symmetry, bead_type, bead_id, molecule_id,
                               residue_id, residue_type, element, mass, charge);

    } else {
      b = container.getBead(bead_id);
    }
    b->setPos(posnm);
  }

  template <bool topology, class T>
  bool ReadFrame(T &container);

  std::ifstream _fl;

  int _line;
};

template <bool topology, class T>
inline bool XYZReader::ReadFrame(T &container) {
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
    int natoms = boost::lexical_cast<int>(line1[0]);
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

      tools::vec pos = tools::vec(boost::lexical_cast<double>(fields[1]),
                                  boost::lexical_cast<double>(fields[2]),
                                  boost::lexical_cast<double>(fields[3]));

      AddAtom<topology, T>(container, fields[0], bead_id, element, pos);
    }
  }
  return !_fl.eof();
}
}  // namespace csg
}  // namespace votca

#endif
