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
#ifndef VOTCA_CSG_GROWRITER_H
#define VOTCA_CSG_GROWRITER_H

#include "../trajectorywriter.h"
#include "growriter.h"
#include <boost/any.hpp>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

template <class Topology_T>
class GROWriter : public TrajectoryWriter {
 public:
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
  tools::UnitConverter converter_;
  double formatDistance_(const double &distance);
  double formatVelocity_(const double &velocity);
  int formatId_(const int &id);
  FILE *_out = nullptr;
};

}  // namespace csg
}  // namespace votca

#include "growriter_priv.h"
#endif  // VOTCA_CSG_GROWRITER_H
