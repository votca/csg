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
#ifndef VOTCA_CSG_GMXTRAJECTORYREADER_H
#define VOTCA_CSG_GMXTRAJECTORYREADER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../trajectoryreader.h"

#include <boost/any.hpp>
#include <cstdlib>
#include <stdexcept>

#include <gromacs/fileio/oenv.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/utility/programcontext.h>
#include <stdexcept>
#include <string>
#include <votca/tools/unitconverter.h>
// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

/**
    \brief class for reading gromacs trajectory files

    This class provides the TrajectoryReader interface and encapsulates the
   trajectory reading function of gromacs

*/
template <class Topology_T>
class GMXTrajectoryReader : public TrajectoryReader {
 public:
  GMXTrajectoryReader() {}

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
  double formatTime_(const double &time);
  double formatDistance_(const double &distance);
  double formatForce_(const double &force);
  double formatVelocity_(const double &velocity);

  tools::UnitConverter converter_;

  std::string _filename;

  // gmx status used in read_first_frame and _read_next_frame;
  t_trxstatus *_gmx_status;
  /// gmx frame
  t_trxframe _gmx_frame;
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/gmxtrajectoryreader_priv.h"
#endif  // VOTCA_CSG_GMXTRAJECTORYREADER_H
