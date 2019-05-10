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
#ifndef VOTCA_CSG_PDBREADER_H
#define VOTCA_CSG_PDBREADER_H

#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include <boost/any.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include <votca/tools/elements.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/unitconverter.h>
namespace votca {
namespace csg {

/**
  brief typename for reading pdb files

  This class provides the Trajectory and Topology reader interface
  for pdb files

*/
template <class Topology_T>
class PDBReader : public csg::TopologyReader, public csg::TrajectoryReader {
 public:
  /// Constuctor
  PDBReader() {}
  /// Destructor
  ~PDBReader() {}
  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  bool ReadTopology(std::string file, boost::any top);

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;

 private:
  std::ifstream _fl;
  bool _topology;
  tools::UnitConverter converter_;

  tools::Elements elements_;
  void formatId(int &id);
  void formatElement(std::string &element, const std::string &atom_type);
  void formatDistance(double &distance);
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/pdbreader_priv.h"

#endif  // VOTCA_CSG_PDBREADER_H
