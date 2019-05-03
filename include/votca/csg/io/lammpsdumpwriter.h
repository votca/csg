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

namespace votca {
namespace csg {

template <class Topology_T>
class LAMMPSDumpWriter : public TrajectoryWriter {
 public:
  void Open(std::string file, bool bAppend = false);
  void Close();

  void RegisteredAt(
      tools::ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(boost::any conf);

 private:
  FILE *_out;
};

template <class Topology_T>
void LAMMPSDumpWriter<Topology_T>::Open(std::string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

template <class Topology_T>
void LAMMPSDumpWriter<Topology_T>::Close() {
  fclose(_out);
}

template <class Topology_T>
void LAMMPSDumpWriter<Topology_T>::Write(boost::any conf_any) {

  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps dump writer, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  Eigen::Matrix3d box = top.getBox();
  fprintf(_out, "ITEM: TIMESTEP\n%i\n", top.getStep());
  fprintf(_out, "ITEM: NUMBER OF ATOMS\n%i\n", (int)top.BeadCount());
  fprintf(_out, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(_out, "0 %f\n0 %f\n0 %f\n", box(0, 0), box(1, 1), box(2, 2));

  fprintf(_out, "ITEM: ATOMS id type x y z");
  bool v = top.HasVel();
  if (v) {
    fprintf(_out, " vx vy vz");
  }
  bool f = top.HasForce();
  if (f) {
    fprintf(_out, " fx fy fz");
  }
  fprintf(_out, "\n");

  std::vector<int> bead_ids = top.getBeadIds();
  // Sort the beads before outputing them
  std::sort(bead_ids.begin(), bead_ids.end());
  for (const int bead_id : bead_ids) {
    typename Topology_T::bead_t *bead = top.getBead(bead_id);
    int bead_type_id = top.getBeadTypeId(bead_id);

    fprintf(_out, "%i %i", bead->getId() + 1, bead_type_id);
    fprintf(_out, " %f %f %f", bead->getPos().x(), bead->getPos().y(),
            bead->getPos().z());
    if (v) {
      fprintf(_out, " %f %f %f", bead->getVel().x(), bead->getVel().y(),
              bead->getVel().z());
    }
    if (f) {
      fprintf(_out, " %f %f %f", bead->getF().x(), bead->getF().y(),
              bead->getF().z());
    }
    fprintf(_out, "\n");
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_LAMMPSDUMPWRITER_H
