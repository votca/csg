/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/getline.h>
#include <votca/tools/unitconverter.h>

// Local VOTCA includes
#include "topologyreader.h"
#include "trajectoryreader.h"

namespace votca {
namespace csg {

/**
    \brief class for reading xyz files

    This class provides the TrajectoryReader + Topology reader interface
    for xyz files

*/
class XYZReader : public TrajectoryReader, public TopologyReader {
 public:
  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;

  XYZReader() = default;
  ~XYZReader() override = default;

  /// open a topology file
  bool ReadTopology(std::string file, Topology &top) override;

  /// open a trajectory file
  bool Open(const std::string &file) override;
  /// read in the first frame
  bool FirstFrame(Topology &top) override;
  /// read in the next frame
  bool NextFrame(Topology &top) override;

  template <class T>
  void ReadFile(T &container) {
    if (!ReadFrame<true, T>(container)) {
      throw std::runtime_error("Reading xyz file '" + file_ + "' failed");
    }
  }

  void Close() override;

 private:
  template <class T>
  Index getContainerSize(T &container) {
    return container.size();
  }

  Index getContainerSize(Topology &container) { return container.BeadCount(); }

  template <bool topology, class T>
  void AddAtom(T &container, std::string name, Index id,
               const Eigen::Vector3d &pos) {
    // the typedef returns the type of the objects the container holds
    using atom =
        typename std::iterator_traits<decltype(container.begin())>::value_type;
    Eigen::Vector3d pos2 = pos * tools::conv::ang2bohr;
    container.push_back(atom(id, name, pos2));
  }

  template <bool topology, class T>
  void AddAtom(Topology &container, std::string name, Index id,
               const Eigen::Vector3d &pos) {
    Bead *b;
    Eigen::Vector3d posnm = pos * tools::conv::ang2nm;
    if (topology) {
      b = container.CreateBead(Bead::spherical,
                               name + boost::lexical_cast<std::string>(id),
                               name, 0, 0, 0);
    } else {
      b = container.getBead(id);
    }
    b->setPos(posnm);
  }

  template <bool topology, class T>
  bool ReadFrame(T &container);

  std::ifstream fl_;
  std::string file_;
  Index line_;
};

template <bool topology, class T>
inline bool XYZReader::ReadFrame(T &container) {
  std::string line;
  tools::getline(fl_, line);
  ++line_;
  if (!fl_.eof()) {
    // read the number of atoms
    std::vector<std::string> line1 = tools::Tokenizer(line, " \t").ToVector();
    if (line1.size() != 1) {
      throw std::runtime_error(
          "First line of xyz file should contain number "
          "of atoms/beads, nothing else.");
    }
    Index natoms = boost::lexical_cast<Index>(line1[0]);
    if (!topology && natoms != getContainerSize(container)) {
      throw std::runtime_error(
          "number of beads in topology and trajectory differ");
    }
    // the title line
    tools::getline(fl_, line);
    ++line_;
    // read atoms
    for (Index i = 0; i < natoms; ++i) {
      tools::getline(fl_, line);
      ++line_;
      if (fl_.eof()) {
        throw std::runtime_error("unexpected end of file in xyz file");
      }
      std::vector<std::string> fields =  tools::Tokenizer(line, " \t").ToVector();
      if (fields.size() != 4) {
        throw std::runtime_error("invalid line " +
                                 boost::lexical_cast<std::string>(line_) +
                                 " in xyz file\n" + line);
      }
      Eigen::Vector3d pos =
          Eigen::Vector3d(boost::lexical_cast<double>(fields[1]),
                          boost::lexical_cast<double>(fields[2]),
                          boost::lexical_cast<double>(fields[3]));

      AddAtom<topology, T>(container, fields[0], i, pos);
    }
  }
  return !fl_.eof();
}
}  // namespace csg
}  // namespace votca

#endif
