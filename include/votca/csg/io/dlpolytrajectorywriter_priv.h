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
#ifndef VOTCA_CSG_DLPOLYTRAJECTORYWRITER_PRIV_H
#define VOTCA_CSG_DLPOLYTRAJECTORYWRITER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
void DLPOLYTrajectoryWriter<Topology_T>::Open(std::string file, bool bAppend)
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

  out_.open(_fname.c_str());
  if (!out_.is_open())
    throw std::ios_base::failure("Error on creating dlpoly file '" + _fname +
                                 "'");
}

template <class Topology_T>
void DLPOLYTrajectoryWriter<Topology_T>::Close() {
  out_.close();
}

template <class Topology_T>
void DLPOLYTrajectoryWriter<Topology_T>::Write(boost::any conf_any) {

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

    out_ << "From VOTCA with love" << std::endl;
    out_ << std::setw(10) << mavecs << std::setw(10) << mpbct << std::setw(10)
         << conf.BeadCount() << std::setw(20) << energy << std::endl;
    Eigen::Matrix3d m = conf.getBox();
    for (int i = 0; i < 3; i++)
      out_ << std::fixed << std::setprecision(10) << std::setw(20)
           << m(i, 0) * scale << std::setw(20) << m(i, 1) * scale
           << std::setw(20) << m(i, 2) * scale << std::endl;

  } else {

    if (nstep == 1) {
      out_ << "From VOTCA with love" << std::endl;
      out_ << std::setw(10) << mavecs << std::setw(10) << mpbct << std::setw(10)
           << conf.BeadCount() << std::endl;
      dstep = conf.getTime() / (double)(conf.getStep());
    }

    out_ << "timestep" << std::setprecision(9) << std::setw(10)
         << conf.getStep() << std::setw(10) << conf.BeadCount() << std::setw(10)
         << mavecs << std::setw(10) << mpbct;
    out_ << std::setprecision(9) << std::setw(12) << dstep << std::setw(12)
         << conf.getTime() << std::endl;

    Eigen::Matrix3d m = conf.getBox();
    for (int i = 0; i < 3; i++)
      out_ << std::setprecision(12) << std::setw(20) << m(i, 0) * scale
           << std::setw(20) << m(i, 1) * scale << std::setw(20)
           << m(i, 2) * scale << std::endl;
  }

  for (int i = 0; static_cast<size_t>(i) < conf.BeadCount(); i++) {
    typename Topology_T::bead_t &bead = conf.getBead(i);

    // AB: DL_POLY needs bead TYPE, not name!

    if (_isConfig) {
      out_ << std::setw(8) << std::left << bead.getType() << std::right
           << std::setw(10) << i + 1 << std::endl;
    } else {
      out_ << std::setw(8) << std::left << bead.getType() << std::right
           << std::setw(10) << i + 1;
      out_ << std::setprecision(6) << std::setw(12) << bead.getMass()
           << std::setw(12) << bead.getQ() << std::setw(12) << "   0.0"
           << std::endl;
    }

    // nm -> Angs
    out_ << resetiosflags(std::ios::fixed) << std::setprecision(12)
         << std::setw(20) << bead.getPos().x() * scale;
    out_ << std::setw(20) << bead.getPos().y() * scale << std::setw(20)
         << bead.getPos().z() * scale << std::endl;

    if (mavecs > 0) {
      if (!bead.HasVel())
        throw std::ios_base::failure(
            "Error: dlpoly frame is supposed to contain velocities, but bead "
            "does not have v-data");

      // nm -> Angs
      out_ << std::setprecision(12) << std::setw(20)
           << bead.getVel().x() * scale << std::setw(20);
      out_ << bead.getVel().y() * scale << std::setw(20)
           << bead.getVel().z() * scale << std::endl;

      if (mavecs > 1) {
        if (!bead.HasF())
          throw std::ios_base::failure(
              "Error: dlpoly frame is supposed to contain forces, but bead "
              "does not have f-data");

        // nm -> Angs
        out_ << std::setprecision(12) << std::setw(20)
             << bead.getF().x() * scale << std::setw(20);
        out_ << bead.getF().y() * scale << std::setw(20)
             << bead.getF().z() * scale << std::endl;
      }
    }
  }
  nstep++;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_DLPOLYTRAJECTORYWRITER_PRIV_H
