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
#ifndef VOTCA_CSG_LAMMPSDUMPREADER_H
#define VOTCA_CSG_LAMMPSDUMPREADER_H

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

#include <votca/tools/constants.h>
#include <votca/tools/getline.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/unitconverter.h>

#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace votca {
namespace csg {

/**
    \brief class for reading lammps dump files

    This class provides the TrajectoryReader + Topology reader interface
    for lammps dump files

*/
template <class Topology_T>
class LAMMPSDumpReader : public TrajectoryReader, public TopologyReader {
 public:
  LAMMPSDumpReader() {}
  ~LAMMPSDumpReader() {}

  /// open a topology file
  bool ReadTopology(const std::string &file, boost::any top);

  /// open a trejectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  void Close();

  /// Assuming units are using 'units real' lammps command
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
  double formatCharge_(const double &charge);
  double formatMass_(const double &mass);

  void ReadTimestep(Topology_T &top, const std::string &itemline);
  void ReadBox(Topology_T &top, const std::string &itemline);
  void ReadNumAtoms(Topology_T &top, const std::string &itemline);
  void ReadAtoms(Topology_T &top, std::string itemline);

  std::ifstream _fl;
  std::string _fname;
  bool read_topology_data_ = false;
  int number_of_atoms_ = 0;
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/lammpsdumpreader_priv.h"
#endif  // VOTCA_CSG_LAMMPSDUMPREADER_H
