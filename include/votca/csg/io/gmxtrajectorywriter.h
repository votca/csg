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
#ifndef VOTCA_CSG_GMXTRAJECTORYWRITER_H
#define VOTCA_CSG_GMXTRAJECTORYWRITER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif
#include <boost/any.hpp>

#include "../trajectorywriter.h"
#include <gromacs/fileio/trxio.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <votca/tools/unitconverter.h>
// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

template <class Topology_T>
class GMXTrajectoryWriter : public TrajectoryWriter {
 public:
  GMXTrajectoryWriter() {}

  void Open(std::string file, bool bAppend = false);
  void Close();
  void Write(boost::any conf);

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
  t_trxstatus *_file;

  tools::UnitConverter converter_;
  double formatTime(const double &time);
  double formatDistance(const double &distance);
  double formatForce(const double &force);
  double formatVelocity(const double &velocity);
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/gmxtrajectorywriter_priv.h"
#endif  // VOTCA_CSG_GMXTRAJECTORYWRITER_H
