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

#ifndef VOTCA_CSG_GMXTRAJECTORYWRITER_H
#define VOTCA_CSG_GMXTRAJECTORYWRITER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../trajectorywriter.h"
#include <gromacs/fileio/trxio.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <string>
// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

template <class Bead_T, class Molecule_T, class Topology_T>
class GMXTrajectoryWriter : public TrajectoryWriter {
 public:
  GMXTrajectoryWriter() {}

  void Open(std::string file, bool bAppend = false);
  void Close();
  void Write(void *conf);

 private:
  t_trxstatus *_file;
};

template <class Bead_T, class Molecule_T, class Topology_T>
void GMXTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Open(std::string file,
                                                               bool bAppend) {
  // char c[1] = bAppend ? "a" : "w";
  _file = open_trx((char *)file.c_str(), "w");
}

template <class Bead_T, class Molecule_T, class Topology_T>
void GMXTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Close() {
  close_trx(_file);
}

template <class Bead_T, class Molecule_T, class Topology_T>
void GMXTrajectoryWriter<Bead_T, Molecule_T, Topology_T>::Write(
    void *uncast_conf) {
  Topology_T *conf = static_cast<Topology_T *>(uncast_conf);
  static int step = 0;
  int N = conf->BeadCount();
  t_trxframe frame;
  rvec *x = new rvec[N];
  rvec *v = NULL;
  rvec *f = NULL;
  Eigen::Matrix3d box = conf->getBox();

  frame.natoms = N;
  frame.bTime = true;
  frame.time = conf->getTime();
  frame.bStep = true;
  frame.step = conf->getStep();
  ;
  frame.x = x;
  frame.bLambda = false;
  frame.bAtoms = false;
  frame.bPrec = false;
  frame.bX = true;
  frame.bF = conf->HasForce();
  frame.bBox = true;
  frame.bV = conf->HasVel();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) frame.box[j][i] = box(i, j);

  for (int i = 0; i < N; ++i) {
    Eigen::Vector3d pos = conf->getBead(i)->getPos();
    x[i][0] = pos.x();
    x[i][1] = pos.y();
    x[i][2] = pos.z();
  }

  if (frame.bV) {
    v = new rvec[N];
    for (int i = 0; i < N; ++i) {
      frame.v = v;
      Eigen::Vector3d vel = conf->getBead(i)->getVel();
      v[i][0] = vel.x();
      v[i][1] = vel.y();
      v[i][2] = vel.z();
    }
  }
  if (frame.bF) {
    f = new rvec[N];
    for (int i = 0; i < N; ++i) {
      frame.f = f;
      Eigen::Vector3d force = conf->getBead(i)->getF();
      f[i][0] = force.x();
      f[i][1] = force.y();
      f[i][2] = force.z();
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
#endif  // VOTCA_CSG_GMXTRAJECTORYWRITER_H
