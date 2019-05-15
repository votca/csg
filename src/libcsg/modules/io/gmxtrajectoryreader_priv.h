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
#ifndef VOTCA_CSG_GMXTRAJECTORYREADER_PRIV_H
#define VOTCA_CSG_GMXTRAJECTORYREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double GMXTrajectoryReader<Topology_T>::formatDistance_(
    const double &distance) {
  return converter_.convert(distance_unit, Topology_T::distance_unit) *
         distance;
}

template <class Topology_T>
double GMXTrajectoryReader<Topology_T>::formatTime_(const double &time) {
  return converter_.convert(time_unit, Topology_T::time_unit) * time;
}

template <class Topology_T>
double GMXTrajectoryReader<Topology_T>::formatVelocity_(
    const double &velocity) {
  return converter_.convert(velocity_unit, Topology_T::velocity_unit) *
         velocity;
}

template <class Topology_T>
double GMXTrajectoryReader<Topology_T>::formatForce_(const double &force) {
  return converter_.convert(force_unit, Topology_T::force_unit) * force;
}

template <class Topology_T>
bool GMXTrajectoryReader<Topology_T>::Open(const std::string &file) {
  _filename = file;
  return true;
}

template <class Topology_T>
void GMXTrajectoryReader<Topology_T>::Close() {
  close_trx(_gmx_status);
}

template <class Topology_T>
bool GMXTrajectoryReader<Topology_T>::FirstFrame(boost::any conf_any) {

  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gmxtrajectory reader first frame, "
        "incorrect topology type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);
  gmx_output_env_t *oenv;
  output_env_init(&oenv, gmx::getProgramContext(), time_ps, FALSE, exvgNONE, 0);
  if (!read_first_frame(oenv, &_gmx_status, (char *)_filename.c_str(),
                        &_gmx_frame, TRX_READ_X | TRX_READ_V | TRX_READ_F))
    throw std::runtime_error(std::string("cannot open ") + _filename);
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m(i, j) = formatDistance_(_gmx_frame.box[j][i]);
  conf.setBox(m);
  conf.setTime(formatTime_(_gmx_frame.time));
  conf.setStep(_gmx_frame.step);

  if (_gmx_frame.natoms != (int)conf.BeadCount())
    throw std::runtime_error(
        "number of beads in trajectory do not match topology");

  for (int i = 0; i < _gmx_frame.natoms; i++) {
    Eigen::Vector3d r = {formatDistance_(_gmx_frame.x[i][XX]),
                         formatDistance_(_gmx_frame.x[i][YY]),
                         formatDistance_(_gmx_frame.x[i][ZZ])};
    conf.getBead(i)->setPos(r);
    if (_gmx_frame.bF) {
      Eigen::Vector3d f = {formatForce_(_gmx_frame.f[i][XX]),
                           formatForce_(_gmx_frame.f[i][YY]),
                           formatForce_(_gmx_frame.f[i][ZZ])};
      conf.getBead(i)->setF(f);
    }
    if (_gmx_frame.bV) {
      Eigen::Vector3d v = {formatVelocity_(_gmx_frame.v[i][XX]),
                           formatVelocity_(_gmx_frame.v[i][YY]),
                           formatVelocity_(_gmx_frame.v[i][ZZ])};
      conf.getBead(i)->setVel(v);
    }
  }
  return true;
}

template <class Topology_T>
bool GMXTrajectoryReader<Topology_T>::NextFrame(boost::any conf_any) {

  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gmx trajectory reader next frame, "
        "incorrect topology type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);
  gmx_output_env_t *oenv;
  output_env_init(&oenv, gmx::getProgramContext(), time_ps, FALSE, exvgNONE, 0);
  if (!read_next_frame(oenv, _gmx_status, &_gmx_frame)) return false;
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m(i, j) = formatDistance_(_gmx_frame.box[j][i]);
  conf.setTime(formatTime_(_gmx_frame.time));
  conf.setStep(_gmx_frame.step);
  conf.setBox(m);

  for (int i = 0; i < _gmx_frame.natoms; i++) {
    Eigen::Vector3d r = {formatDistance_(_gmx_frame.x[i][XX]),
                         formatDistance_(_gmx_frame.x[i][YY]),
                         formatDistance_(_gmx_frame.x[i][ZZ])};
    conf.getBead(i)->setPos(r);
    if (_gmx_frame.bF) {
      Eigen::Vector3d f = {formatForce_(_gmx_frame.f[i][XX]),
                           formatForce_(_gmx_frame.f[i][YY]),
                           formatForce_(_gmx_frame.f[i][ZZ])};
      conf.getBead(i)->setF(f);
    }
    if (_gmx_frame.bV) {
      Eigen::Vector3d v = {formatVelocity_(_gmx_frame.v[i][XX]),
                           formatVelocity_(_gmx_frame.v[i][YY]),
                           formatVelocity_(_gmx_frame.v[i][ZZ])};
      conf.getBead(i)->setVel(v);
    }
  }
  return true;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GMXTRAJECTORYREADER_PRIV_H
