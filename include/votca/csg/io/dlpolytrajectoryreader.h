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
#ifndef VOTCA_CSG_DLPOLYTRAJECTORYREADER_H
#define VOTCA_CSG_DLPOLYTRAJECTORYREADER_H

#include "dlpolytrajectoryreader.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <boost/filesystem/convenience.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../boundarycondition.h"
#include "../trajectoryreader.h"
#include <votca/tools/constants.h>
#include <votca/tools/getline.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

/**
    \brief class for reading dlpoly trajectory and configuration files

    This class encapsulates the dlpoly trajectory and configuration reading
   function and provides an interface to fill a topology class

*/

template <class Topology_T>
class DLPOLYTrajectoryReader : public TrajectoryReader {
 public:
  /// open original trajectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(boost::any &conf);
  /// read in the next frame
  bool NextFrame(boost::any &conf);
  /// close original trajectory file
  void Close();

  /// set/get the original configuration or trajectory file name:
  /// <name>.dlpc/<name>.dlph (convention: ".dlpc"="CONFIG", ".dlph"="HISTORY")
  void setFname(std::string name) {
    _fname = name;
    return;
  }
  std::string getFname() { return _fname; }

  /// set/check the flag for the read-in file as configuration, i.e. not
  /// trajectory format
  void setIsConfig(bool isConf) {
    _isConfig = isConf;
    return;
  }
  bool getIsConfig() { return _isConfig; }

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;
  const tools::MassUnit mass_unit = tools::MassUnit::atomic_mass_units;
  const tools::TimeUnit time_unit = tools::TimeUnit::picoseconds;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::EnergyUnit energy_unit = tools::EnergyUnit::joules_per_mole;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::angstroms_per_picosecond;
  const tools::ForceUnit force_unit =
      tools::ForceUnit::kilojoules_per_mole_angstrom;

 private:
  tools::UnitConverter converter_;
  double formatDistance_(const double &distance) const;
  Eigen::Vector3d formatDistance_(const Eigen::Vector3d &distances) const;
  Eigen::Vector3d formatVelocity_(const Eigen::Vector3d &velocity) const;
  Eigen::Vector3d formatForce_(const Eigen::Vector3d &force) const;
  double formatTime_(const double &time) const;
  std::ifstream _fl;
  std::string _fname;
  bool _first_frame;
  bool _isConfig;
};

}  // namespace csg
}  // namespace votca

#include "dlpolytrajectoryreader_priv.h"
#endif  // VOTCA_CSG_DLPOLYTRAJECTORYREADER_H
