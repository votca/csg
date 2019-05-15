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
#ifndef VOTCA_CSG_LAMMPSDUMPWRITER_H
#define VOTCA_CSG_LAMMPSDUMPWRITER_H

#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <stdio.h>
#include <string>
#include <votca/tools/objectfactory.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

template <class Topology_T>
class LAMMPSDumpWriter : public TrajectoryWriter {
 public:
  void Open(std::string file, bool bAppend = false);
  void Close();

  //  void RegisteredAt(
  //      tools::ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(boost::any conf);

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;
  const tools::TimeUnit time_unit = tools::TimeUnit::femtoseconds;
  const tools::MassUnit mass_unit = tools::MassUnit::grams_per_mole;
  const tools::EnergyUnit energy_unit =
      tools::EnergyUnit::kilocalories_per_mole;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::ForceUnit force_unit =
      tools::ForceUnit::kilocalories_per_mole_ansgtrom;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::angstroms_per_femtosecond;

 private:
  tools::UnitConverter converter_;
  int formatId_(const int &id);
  double formatDistance_(const double &distance);
  double formatForce_(const double &force);
  double formatVelocity_(const double &velocity);
  FILE *_out;
};

}  // namespace csg
}  // namespace votca

#include "lammpsdumpwriter_priv.h"

#endif  // VOTCA_CSG_LAMMPSDUMPWRITER_H
