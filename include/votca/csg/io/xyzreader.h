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
#ifndef VOTCA_CSG_XYZREADER_H
#define VOTCA_CSG_XYZREADER_H

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
  void ReadFile(T &container);

 private:
  template <class T>
  size_t getContainerSize(T &container);

  size_t getContainerSize(Topology_T &top);

  template <bool topology, class T>
  void AddAtom(T &container, tools::StructureParameters &params);

  template <bool topology>
  void AddAtom(Topology_T &container, tools::StructureParameters &params);

  template <bool topology, class T>
  bool ReadFrame(T &container);

  std::ifstream _fl;
  std::string _file;
  int _line;
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/xyzreader_priv.h"

#endif
