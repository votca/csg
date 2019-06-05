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
#ifndef VOTCA_CSG_PDBWRITER_H
#define VOTCA_CSG_PDBWRITER_H

#include "../molecule.h"
#include "../trajectorywriter.h"
#include <Eigen/Dense>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <stdio.h>
#include <type_traits>
#include <votca/tools/constants.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/unitconverter.h>
namespace votca {
namespace csg {

/**
 * @brief
 *
 * Note that by default the coordinates stored in the pdb file are in units
 * of Angstroms while the following repositories have default units of:
 *
 * csg: nm
 * xtp: bohr
 *
 * Further note that pdb files store the ids, serical numbers and indices
 * starting at 1, votca uses ids, indices and serial numbers starting at 0.
 */
template <class Topology_T>
class PDBWriter : public TrajectoryWriter {
 public:
  void Write(boost::any conf);

  template <class T>
  void WriteContainer(T &container);

  void WriteHeader(std::string header);

  template <class T>
  void WriteBox(const T &cont);

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;

 private:
  tools::UnitConverter converter_;

  void formatType_(std::string &atom_type);
  void formatElement_(std::string &element);
  /**
   * @brief Get the residue type
   *
   * In the case that the residue type is unknown it will follow the pdb
   * convention and assign it a value of "UNK", see the url for more info:
   *
   * www.wwpdb.org/documentation/file-format-content/format33/sect4.html
   *
   * @tparam Atom
   * @param atom
   *
   * @return
   */
  void formatResidueType_(std::string &restype);
  /**
   * @brief Gets the id of the internal VOTCA class instance and adds one
   *
   * This is necessary becasue ids in VOTCA are stored starting with id 0,
   * whereas the pdb files store the ids starting at integer values of 1.
   *
   * @tparam T
   * @param item
   *
   * @return equivalent id for pdb file
   */
  void formatId_(int &id) noexcept;
  void formatResId_(int &resId) noexcept;
  void formatPos_(Eigen::Vector3d &pos);

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

#include "pdbwriter_priv.h"

#endif  // VOTCA_CSG_PDBWRITER_H
