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
#ifndef VOTCA_CSG_GMXTOPOLOGYREADER_H
#define VOTCA_CSG_GMXTOPOLOGYREADER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../topologyreader.h"
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/types.h>
#include <votca/tools/unitconverter.h>

#include <gromacs/fileio/tpxio.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/topology/atoms.h>
#include <gromacs/topology/topology.h>

// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

/**
    \brief reader for gromacs topology files

    This class encapsulates the gromacs reading functions and provides an
   interface to fill a topolgy class

*/
template <class Topology_T>
class GMXTopologyReader : public TopologyReader {
 public:
  GMXTopologyReader() {}

  /// read a topology file
  bool ReadTopology(const std::string &file, boost::any top);

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
  tools::Elements elements_;
  tools::UnitConverter converter_;

  std::string formatElement_(const std::string &bead_type);
  double formatDistance_(const double &distance);
  double formatCharge_(const double &charge);
  double formatMass_(const double &mass);
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/gmxtopologyreader_priv.h"
#endif  // VOTCA_CSG_GMXTOPOLOGYREADER_H
