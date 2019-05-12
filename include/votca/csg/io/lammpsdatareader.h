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
#ifndef VOTCA_CSG_LAMMPSDATAREADER_H
#define VOTCA_CSG_LAMMPSDATAREADER_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../topologyreader.h"
#include "../trajectoryreader.h"

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

#include <votca/tools/elements.h>
#include <votca/tools/getline.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

/**
    \brief class for reading lammps data files

    This class provides the TrajectoryReader + Topology reader interface
    for lammps data files

*/
template <class Topology_T>
class LAMMPSDataReader : public TrajectoryReader, public TopologyReader {
 public:
  LAMMPSDataReader() {}
  ~LAMMPSDataReader() {}

  /// open, read and close topology file
  bool ReadTopology(const std::string &file, boost::any top);

  /// open a trajectory file
  bool Open(const std::string &file);
  /// read in the first frame of trajectory file
  bool FirstFrame(boost::any top);
  /// read in the next frame of trajectory file
  bool NextFrame(boost::any top);
  /// close the topology file
  void Close();

  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;
  const tools::TimeUnit time_unit = tools::TimeUnit::femtoseconds;
  const tools::MassUnit mass_unit = tools::MassUnit::grams_per_mole;
  const tools::EnergyUnit energy_unit =
      tools::EnergyUnit::kilocalories_per_mole;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::ForceUnit force_unit =
      tools::ForceUnit::kilocalories_per_mole_ansgtrom;

 private:
  std::ifstream fl_;
  std::string fname_;
  bool topology_;

  tools::UnitConverter converter_;
  std::map<std::string, std::vector<std::vector<std::string>>> data_;

  // int - atom type starting index of 0
  // .at(0) - element symbol or bead
  // .at(1) - atom name may be the same as the element symbol or bead depending
  //          on if there is more than one atom type for a given element
  std::map<int, std::vector<std::string>> atomtypes_;

  // String is the type .e.g. "atom","bond" etc
  // int is the number of different types
  std::map<std::string, int> numberOfDifferentTypes_;

  // String is the type .e.g. "atom", "bond" etc
  // int is the number of them
  std::map<std::string, int> numberOf_;

  // First int is the molecule id
  // Second int is the molecule ptr
  std::map<int, Molecule *> molecules_;

  // First int is the atom id second the molecule id
  std::map<int, int> atomIdToMoleculeId_;

  // First int is the atom id second the atom index
  std::map<int, int> atomIdToIndex_;

  void formatDistance(double &distance);
  void formatId(int &id);
  void formatMass(double &mass);
  void formatCharge(double &mass);

  bool MatchOneFieldLabel_(std::vector<std::string> fields, Topology_T &top);
  bool MatchTwoFieldLabels_(std::vector<std::string> fields, Topology_T &top);
  bool MatchThreeFieldLabels_(std::vector<std::string> fields, Topology_T &top);
  bool MatchFourFieldLabels_(std::vector<std::string> fields, Topology_T &top);
  bool MatchFieldsTimeStepLabel_(std::vector<std::string> fields,
                                 Topology_T &top);

  void ReadBox_(std::vector<std::string> fields, Topology_T &top);
  void SortIntoDataGroup_(std::string tag);
  void ReadNumTypes_(std::vector<std::string> fields, std::string type);

  void ReadNumOfAtoms_(std::vector<std::string> fields, Topology_T &top);
  void ReadNumOfBonds_(std::vector<std::string> fields);
  void ReadNumOfAngles_(std::vector<std::string> fields);
  void ReadNumOfDihedrals_(std::vector<std::string> fields);
  void ReadNumOfImpropers_(std::vector<std::string> fields);

  void ReadAtoms_(Topology_T &top);
  void ReadBonds_(Topology_T &top);
  void ReadAngles_(Topology_T &top);
  void ReadDihedrals_(Topology_T &top);
  void ReadImpropers_(Topology_T &top);

  enum lammps_format {
    style_angle_bond_molecule = 0,
    style_atomic = 1,
    style_full = 2
  };
  lammps_format determineDataFileFormat_(std::string line);

  /**
   * \brief Determines atom and bead types based on masses in lammps files
   *
   * The purpose of this function is to take lammps output where there are more
   * than a single atom type of the same element. For instance there may be 4
   * atom types with mass of 12.01. Well this means that they are all carbon but
   * are treated differently in lammps. It makes since to keep track of this. If
   * a mass cannot be associated with an element we will assume it is pseudo
   * atom or course grained watom which we will represent with as a bead. So
   * when creating the atom names we will take into account so say we have the
   * following masses in the lammps .data file:
   *
   * Masses
   *
   * 1 1.0
   * 2 12.01
   * 3 12.01
   * 4 16.0
   * 5 12.01
   * 6 15.2
   * 7 12.8
   * 8 15.2
   *
   * Then we would translate this to the following atom names
   * 1 H
   * 2 C1
   * 3 C2
   * 4 O
   * 5 C3
   * 6 Bead1 BeadType 1
   * 7 Bead2 BeadType
   * 8 Bead1 BeadType 2
   *
   * Note that we do not append a number if it is singular, in such cases the
   * element and the atom name is the same.
   **/
  void InitializeAtomAndBeadTypes_();
  std::map<std::string, double> determineBaseNameAssociatedWithMass_();
  std::map<std::string, int> determineAtomAndBeadCountBasedOnMass_(
      std::map<std::string, double> baseNamesAndMasses);

  /// Remove comments and padded spaces
  void stripLine(std::string &line);
  //  std::vector<std::string> TrimCommentsFrom_(std::vector<std::string>
  //  fields); void ltrim_(std::string &s); void rtrim_(std::string &s); void
  //  trim_(std::string &s);
  bool withinTolerance_(double value1, double value2, double tolerance);
  std::string getStringGivenDoubleAndMap_(
      double value, std::map<std::string, double> nameValue, double tolerance);
};

}  // namespace csg
}  // namespace votca

#include "../../../../src/libcsg/modules/io/lammpsdatareader_priv.h"
#endif  // VOTCA_CSG_LAMMPSDATAREADER_H
