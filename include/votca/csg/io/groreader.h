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
#include <votca/tools/unitconverter.h>

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
  bool ReadTopology(const std::string &file, boost::any top);

  /// open a trejectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  void Close();

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::nanometers;
  const tools::MassUnit mass_unit = tools::MassUnit::atomic_mass_units;
  const tools::TimeUnit time_unit = tools::TimeUnit::picoseconds;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::EnergyUnit energy_unit = tools::EnergyUnit::kilojoules_per_mole;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::nanometers_per_picosecond;
  const tools::ForceUnit force_unit =
      tools::ForceUnit::kilojoules_per_mole_nanometer;

 private:
  std::ifstream _fl;
  tools::UnitConverter converter_;
  tools::Elements elements_;
  double formatDistance_(const double &distance);
  double formatVelocity_(const double &velocity);
  double formatMass_(const std::string &atom_name);
  int formatId_(const int &id);
  std::string formatElement_(const std::string &element);
  bool _topology;
};

}  // namespace csg
}  // namespace votca

#include "groreader_priv.h"
#endif  // VOTCA_CSG_GROREADER_H
