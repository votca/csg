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
#ifndef VOTCA_CSG_GMXTRAJECTORYWRITER_PRIV_H
#define VOTCA_CSG_GMXTRAJECTORYWRITER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double GMXTrajectoryWriter<Topology_T>::formatTime(const double &time) {
  return converter_.convert(Topology_T::units::time_unit, time_unit) * time;
}

template <class Topology_T>
double GMXTrajectoryWriter<Topology_T>::formatDistance(const double &distance) {
  return converter_.convert(Topology_T::units::distance_unit, distance_unit) *
         distance;
}

template <class Topology_T>
double GMXTrajectoryWriter<Topology_T>::formatForce(const double &force) {
  return converter_.convert(Topology_T::units::force_unit, force_unit) * force;
}

template <class Topology_T>
double GMXTrajectoryWriter<Topology_T>::formatVelocity(const double &velocity) {
  return converter_.convert(Topology_T::units::velocity_unit, velocity_unit) *
         velocity;
}

template <class Topology_T>
void GMXTrajectoryWriter<Topology_T>::Open(std::string file, bool bAppend) {
  // char c[1] = bAppend ? "a" : "w";
  _file = open_trx((char *)file.c_str(), "w");
}

template <class Topology_T>
void GMXTrajectoryWriter<Topology_T>::Close() {
  close_trx(_file);
}

template <class Topology_T>
void GMXTrajectoryWriter<Topology_T>::Write(boost::any conf_any) {

  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gmx trajectory writer, incorrect "
        "topology type provided.");
  }
  Topology_T &conf = *boost::any_cast<Topology_T *>(conf_any);

  static int step = 0;
  int N = conf.BeadCount();
  t_trxframe frame;
  rvec *x = new rvec[N];
  rvec *v = NULL;
  rvec *f = NULL;
  Eigen::Matrix3d box = conf.getBox();

  frame.natoms = N;
  frame.bTime = true;
  frame.time = formatTime(conf.getTime());
  frame.bStep = true;
  frame.step = conf.getStep();

  frame.x = x;
  frame.bLambda = false;
  frame.bAtoms = false;
  frame.bPrec = false;
  frame.bX = true;
  frame.bF = conf.HasForce();
  frame.bBox = true;
  frame.bV = conf.HasVel();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) frame.box[j][i] = formatDistance(box(i, j));

  for (int i = 0; i < N; ++i) {
    Eigen::Vector3d pos = conf.getBead(i).getPos();
    x[i][0] = formatDistance(pos.x());
    x[i][1] = formatDistance(pos.y());
    x[i][2] = formatDistance(pos.z());
  }

  if (frame.bV) {
    v = new rvec[N];
    for (int i = 0; i < N; ++i) {
      frame.v = v;
      Eigen::Vector3d vel = conf.getBead(i).getVel();
      v[i][0] = formatVelocity(vel.x());
      v[i][1] = formatVelocity(vel.y());
      v[i][2] = formatVelocity(vel.z());
    }
  }
  if (frame.bF) {
    f = new rvec[N];
    for (int i = 0; i < N; ++i) {
      frame.f = f;
      Eigen::Vector3d force = conf.getBead(i).getF();
      f[i][0] = formatForce(force.x());
      f[i][1] = formatForce(force.y());
      f[i][2] = formatForce(force.z());
    }
  }

  write_trxframe(_file, &frame, NULL);

  step++;
  delete[] x;
  if (frame.bV) delete[] v;
  if (frame.bF) delete[] f;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GMXTRAJECTORYWRITER_PRIV_H
