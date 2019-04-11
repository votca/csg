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

#ifndef _gmxtrajectoryreader_H
#define _gmxtrajectoryreader_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../trajectoryreader.h"

#include <cstdlib>
#include <gromacs/fileio/oenv.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/utility/programcontext.h>
#include <stdexcept>
#include <string>
#include <votca/tools/matrix.h>
// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

/**
    \brief class for reading gromacs trajectory files

    This class provides the TrajectoryReader interface and encapsulates the
   trajectory reading function of gromacs

*/
template <class Bead_T, class Molecule_T, class Topology_T>
class GMXTrajectoryReader : public TrajectoryReader {
 public:
  GMXTrajectoryReader() {}

  /// open a trejectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(void *top);
  /// read in the next frame
  bool NextFrame(void *top);

  void Close();

 private:
  std::string _filename;

  // gmx status used in read_first_frame and _read_next_frame;
  t_trxstatus *_gmx_status;
  /// gmx frame
  t_trxframe _gmx_frame;
};

template <class Bead_T, class Molecule_T, class Topology_T>
bool GMXTrajectoryReader<Bead_T, Molecule_T, Topology_T>::Open(
    const std::string &file) {
  _filename = file;
  return true;
}

template <class Bead_T, class Molecule_T, class Topology_T>
void GMXTrajectoryReader<Bead_T, Molecule_T, Topology_T>::Close() {
  close_trx(_gmx_status);
}

template <class Bead_T, class Molecule_T, class Topology_T>
bool GMXTrajectoryReader<Bead_T, Molecule_T, Topology_T>::FirstFrame(
    void *conf) {
  Topology_T *conf_cast = static_cast<Topology_T *>(conf);
  gmx_output_env_t *oenv;
  output_env_init(&oenv, gmx::getProgramContext(), time_ps, FALSE, exvgNONE, 0);
  if (!read_first_frame(oenv, &_gmx_status, (char *)_filename.c_str(),
                        &_gmx_frame, TRX_READ_X | TRX_READ_V | TRX_READ_F))
    throw std::runtime_error(std::string("cannot open ") + _filename);
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m(i, j) = _gmx_frame.box[j][i];
  conf_cast->setBox(m);
  conf_cast->setTime(_gmx_frame.time);
  conf_cast->setStep(_gmx_frame.step);

  if (_gmx_frame.natoms != (int)conf_cast->BeadCount())
    throw std::runtime_error(
        "number of beads in trajectory do not match topology");

  for (int i = 0; i < _gmx_frame.natoms; i++) {
    Eigen::Vector3d r = {_gmx_frame.x[i][XX], _gmx_frame.x[i][YY],
                         _gmx_frame.x[i][ZZ]};
    conf_cast->getBead(i)->setPos(r);
    if (_gmx_frame.bF) {
      Eigen::Vector3d f = {_gmx_frame.f[i][XX], _gmx_frame.f[i][YY],
                           _gmx_frame.f[i][ZZ]};
      conf_cast->getBead(i)->setF(f);
    }
    if (_gmx_frame.bV) {
      Eigen::Vector3d v = {_gmx_frame.v[i][XX], _gmx_frame.v[i][YY],
                           _gmx_frame.v[i][ZZ]};
      conf_cast->getBead(i)->setVel(v);
    }
  }
  return true;
}

template <class Bead_T, class Molecule_T, class Topology_T>
bool GMXTrajectoryReader<Bead_T, Molecule_T, Topology_T>::NextFrame(
    void *conf) {
  Topology_T *conf_cast = static_cast<Topology_T *>(conf);
  gmx_output_env_t *oenv;
  output_env_init(&oenv, gmx::getProgramContext(), time_ps, FALSE, exvgNONE, 0);
  if (!read_next_frame(oenv, _gmx_status, &_gmx_frame)) return false;
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m(i, j) = _gmx_frame.box[j][i];
  conf_cast->setTime(_gmx_frame.time);
  conf_cast->setStep(_gmx_frame.step);
  conf_cast->setBox(m);

  for (int i = 0; i < _gmx_frame.natoms; i++) {
    Eigen::Vector3d r = {_gmx_frame.x[i][XX], _gmx_frame.x[i][YY],
                         _gmx_frame.x[i][ZZ]};
    conf_cast->getBead(i)->setPos(r);
    if (_gmx_frame.bF) {
      Eigen::Vector3d f = {_gmx_frame.f[i][XX], _gmx_frame.f[i][YY],
                           _gmx_frame.f[i][ZZ]};
      conf_cast->getBead(i)->setF(f);
    }
    if (_gmx_frame.bV) {
      Eigen::Vector3d v = {_gmx_frame.v[i][XX], _gmx_frame.v[i][YY],
                           _gmx_frame.v[i][ZZ]};
      conf_cast->getBead(i)->setVel(v);
    }
  }
  return true;
}

}  // namespace csg
}  // namespace votca
#endif /* _gmxtrajectoryreader_H */
