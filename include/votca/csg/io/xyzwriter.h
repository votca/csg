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
#ifndef VOTCA_CSG_XYZWRITER_H
#define VOTCA_CSG_XYZWRITER_H

#include "../csgtopology.h"
#include "../molecule.h"
#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <stdio.h>
#include <votca/tools/constants.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

template <class Topology_T>
class XYZWriter : public TrajectoryWriter {
 public:
  void Write(boost::any conf);

  template <class T>
  void WriteContainer(T &container);

  void WriteHeader(std::string header, int number_atoms);

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;

 private:
  tools::UnitConverter converter_;

  void formatType(std::string &atom_type);
  void formatPosition(Eigen::Vector3d &position);

  template <class T>
  void WriteContainer_(T &container);

  template <class T>
  int getAtomCount(T &container);

  int getAtomCount(Topology_T &container);

  template <class T>
  T &getIterable(T &container);

  std::vector<typename Molecule::bead_t> getIterable(Molecule &container);

  template <typename T>
  T *ptr(T *obj);

  template <typename T>
  T *ptr(T &obj);
};

}  // namespace csg
}  // namespace votca

#include "xyzwriter_priv.h"

#endif
