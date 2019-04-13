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

#ifndef VOTCA_CSG_DLPOLYTRAJECTORYWRITER_H
#define VOTCA_CSG_DLPOLYTRAJECTORYWRITER_H

#include "../boundarycondition.h"
#include "../trajectorywriter.h"
#include <boost/any.hpp>
#include <boost/filesystem/convenience.hpp>
#include <iomanip>
#include <string>

namespace votca {
namespace csg {

/**
    \brief class for writing dlpoly trajectory and configuration files

    This class encapsulates the dlpoly trajectory and configuration writing
   function

*/
template <class Bead_T, class Molecule_T, class Topology_T>
class DLPOLYTrajectoryWriter : public TrajectoryWriter {
 public:
  // open transformed trajectory file
  void Open(std::string file, bool bAppend = false);
  // close transformed trajectory file
  void Close();
  // write a frame into transformed trajectory file
  void Write(boost::any conf);

  /// set/get the created configuration or trajectory file name:
  /// <name>.dlpc or <name>.dlph (convention: ".dlpc"="CONFIG_CGV",
  /// ".dlph"="HISTORY_CGV")
  void setFname(std::string name) {
    _fname = name;
    return;
  }
  std::string getFname() { return _fname; }

  /// set/check the flag for the created file as configuration, i.e. not
  /// trajectory format
  void setIsConfig(bool isConf) {
    _isConfig = isConf;
    return;
  }
  bool getIsConfig() { return _isConfig; }

 private:
  std::ofstream _fl;
  std::string _fname;
  bool _isConfig;
};

template <class Bead_T, class Molecule_T, class Topology_T>
void DLPOLYTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Open(
    std::string file, bool bAppend)
// open/create a dlpoly configuration or trajectory file
// NOTE: allowed file naming - <name>.dlpc or <name>.dlph (convention:
// ".dlpc"="CONFIG_CGV", ".dlph"="HISTORY_CGV")
{
  if (bAppend)
    throw std::runtime_error(
        "Error: appending to dlpoly files not implemented");

  boost::filesystem::path filepath(file.c_str());
  std::string out_name = "HISTORY_CGV";

  if (boost::filesystem::extension(filepath).size() == 0) {

    throw std::ios_base::failure("Error on creating dlpoly file '" + file +
                                 "' - extension is expected, .dlph or .dlpc");

  } else if (boost::filesystem::extension(filepath) == ".dlpc") {

    _isConfig = true;
    out_name = "CONFIG_CGV";

  } else if (boost::filesystem::extension(filepath) == ".dlph") {

    _isConfig = false;

  } else {
    throw std::ios_base::failure("Error on creating dlpoly file '" + file +
                                 "' - wrong extension, use .dlph or .dlpc");
  }

  if (boost::filesystem::basename(filepath).size() == 0) {
    if (filepath.parent_path().string().size() == 0) {
      _fname = out_name;
    } else {
      _fname = filepath.parent_path().string() + "/" + out_name;
    }
  } else {
    _fname = file;
  }

  _fl.open(_fname.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on creating dlpoly file '" + _fname +
                                 "'");
}

template <class Bead_T, class Molecule_T, class Topology_T>
void DLPOLYTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Close() {
  _fl.close();
}

template <class Bead_T, class Molecule_T, class Topology_T>
void DLPOLYTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Write(
    boost::any conf_any) {

  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot write topology using dlpolytopology writer, incorrect "
        "topology type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);
  static int nstep = 1;
  static double dstep = 0.0;
  const double scale = 10.0;  // nm -> A factor
  int mavecs = 0;
  int mpbct = 0;
  double energy = 0.0;

  if (conf.HasForce() && conf.HasVel()) {
    mavecs = 2;
  } else if (conf.HasVel()) {
    mavecs = 1;
  }

  if (conf.getBoxType() == BoundaryCondition::typeOrthorhombic) mpbct = 2;
  if (conf.getBoxType() == BoundaryCondition::typeTriclinic) mpbct = 3;

  if (_isConfig) {

    _fl << "From VOTCA with love" << std::endl;
    _fl << std::setw(10) << mavecs << std::setw(10) << mpbct << std::setw(10)
        << conf.BeadCount() << std::setw(20) << energy << std::endl;
    Eigen::Matrix3d m = conf.getBox();
    for (int i = 0; i < 3; i++)
      _fl << std::fixed << std::setprecision(10) << std::setw(20)
          << m(i, 0) * scale << std::setw(20) << m(i, 1) * scale
          << std::setw(20) << m(i, 2) * scale << std::endl;

  } else {

    if (nstep == 1) {
      _fl << "From VOTCA with love" << std::endl;
      _fl << std::setw(10) << mavecs << std::setw(10) << mpbct << std::setw(10)
          << conf.BeadCount() << std::endl;
      dstep = conf.getTime() / (double)(conf.getStep());
    }

    _fl << "timestep" << std::setprecision(9) << std::setw(10) << conf.getStep()
        << std::setw(10) << conf.BeadCount() << std::setw(10) << mavecs
        << std::setw(10) << mpbct;
    _fl << std::setprecision(9) << std::setw(12) << dstep << std::setw(12)
        << conf.getTime() << std::endl;

    Eigen::Matrix3d m = conf.getBox();
    for (int i = 0; i < 3; i++)
      _fl << std::setprecision(12) << std::setw(20) << m(i, 0) * scale
          << std::setw(20) << m(i, 1) * scale << std::setw(20)
          << m(i, 2) * scale << std::endl;
  }

  for (int i = 0; static_cast<size_t>(i) < conf.BeadCount(); i++) {
    Bead_T *bead = conf.getBead(i);

    // AB: DL_POLY needs bead TYPE, not name!

    if (_isConfig) {
      _fl << std::setw(8) << std::left << bead->getType() << std::right
          << std::setw(10) << i + 1 << std::endl;
    } else {
      _fl << std::setw(8) << std::left << bead->getType() << std::right
          << std::setw(10) << i + 1;
      _fl << std::setprecision(6) << std::setw(12) << bead->getMass()
          << std::setw(12) << bead->getQ() << std::setw(12) << "   0.0"
          << std::endl;
    }

    // nm -> Angs
    _fl << resetiosflags(std::ios::fixed) << std::setprecision(12)
        << std::setw(20) << bead->getPos().x() * scale;
    _fl << std::setw(20) << bead->getPos().y() * scale << std::setw(20)
        << bead->getPos().z() * scale << std::endl;

    if (mavecs > 0) {
      if (!bead->HasVel())
        throw std::ios_base::failure(
            "Error: dlpoly frame is supposed to contain velocities, but bead "
            "does not have v-data");

      // nm -> Angs
      _fl << std::setprecision(12) << std::setw(20)
          << bead->getVel().x() * scale << std::setw(20);
      _fl << bead->getVel().y() * scale << std::setw(20)
          << bead->getVel().z() * scale << std::endl;

      if (mavecs > 1) {
        if (!bead->HasF())
          throw std::ios_base::failure(
              "Error: dlpoly frame is supposed to contain forces, but bead "
              "does not have f-data");

        // nm -> Angs
        _fl << std::setprecision(12) << std::setw(20)
            << bead->getF().x() * scale << std::setw(20);
        _fl << bead->getF().y() * scale << std::setw(20)
            << bead->getF().z() * scale << std::endl;
      }
    }
  }
  nstep++;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_DLPOLYTRAJECTORYWRITER_H
