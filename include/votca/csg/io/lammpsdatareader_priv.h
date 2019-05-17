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
#ifndef VOTCA_CSG_LAMMPSDATAREADER_PRIV_H
#define VOTCA_CSG_LAMMPSDATAREADER_PRIV_H
namespace votca {
namespace csg {

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::formatDistance(double &distance) {
  distance *=
      converter_.convert(this->distance_unit, Topology_T::distance_unit);
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::formatId(int &id) {
  --id;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::formatMass(double &mass) {
  mass *= converter_.convert(this->mass_unit, Topology_T::mass_unit);
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::formatCharge(double &charge) {
  charge *= converter_.convert(this->charge_unit, Topology_T::charge_unit);
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::stripLine(std::string &line) {
  // Hash (#) comments
  line = line.substr(0, line.find("#", 0));
  boost::algorithm::trim(line);
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::withinTolerance_(double value1,
                                                    double value2,
                                                    double tolerance) {
  return abs(value1 - value2) / std::min(value1, value2) < tolerance;
}

template <class Topology_T>
std::string LAMMPSDataReader<Topology_T>::getStringGivenDoubleAndMap_(
    double value, std::map<std::string, double> nameValue, double tolerance) {

  for (auto string_value_pair : nameValue) {
    if (withinTolerance_(value, string_value_pair.second, tolerance)) {
      return string_value_pair.first;
    }
  }
  throw std::runtime_error(
      "getStringGivenDoubleAndMap_ function fails. This method "
      "is meant to be passed a double that is to be matched within a tolerance"
      " with a double in a map<string,double> and then return the string. It is"
      " likely that none of the doubles were a close enough match.");
}

/*****************************************************************************
 * Public Facing Methods                                                     *
 *****************************************************************************/

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::ReadTopology(const std::string &file,
                                                boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps data reader read topology, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);

  std::cout << std::endl;
  std::cout << "WARNING: The votca lammps data reader is only able to read ";
  std::cout << "lammps files formatted in the following styles:" << std::endl;
  std::cout << "angle" << std::endl;
  std::cout << "atom" << std::endl;
  std::cout << "bond" << std::endl;
  std::cout << "full" << std::endl;
  std::cout << "molecule" << std::endl;
  std::cout << std::endl;
  std::cout << "These styles use the following formats in the atom block:"
            << std::endl;
  std::cout << "atom-ID molecule-ID atom-type charge x y z" << std::endl;
  std::cout << "atom-ID molecule-ID atom-type charge x y z nx ny nz"
            << std::endl;
  std::cout << "atom-ID molecule-ID atom-type x y z" << std::endl;
  std::cout << "atom-ID molecule-ID atom-type x y z nx ny nz" << std::endl;
  std::cout << "atom-ID atom-type x y z" << std::endl;
  std::cout << "atom-ID atom-type x y z nx ny nz" << std::endl;
  std::cout << std::endl;
  std::cout << "WARNING: Currently the votca lammps data reader only suppots ";
  std::cout << "lammps units specified by the 'units real' ";
  std::cout << "command." << std::endl;
  std::cout << "mass: grams/mole" << std::endl;
  std::cout << "distance: angstroms" << std::endl;
  std::cout << "time: femtoseconds" << std::endl;
  std::cout << "energy: Kcal/mole" << std::endl;
  std::cout << "charge: e" << std::endl;
  std::cout << std::endl;

  topology_ = true;
  top.Cleanup();
  fl_.open(file.c_str());
  if (!fl_.is_open())
    throw std::ios_base::failure("Error on open topology file: " + file);

  fname_ = file;

  NextFrame(top_any);

  fl_.close();

  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::Open(const std::string &file) {
  fl_.open(file.c_str());
  if (!fl_.is_open())
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  fname_ = file;
  return true;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::Close() {
  fl_.close();
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::FirstFrame(boost::any top) {
  topology_ = false;
  NextFrame(top);
  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::NextFrame(boost::any top_any) {
  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps data reader next frame, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  std::string line;
  getline(fl_, line);
  while (!fl_.eof()) {

    // Remove all comments
    stripLine(line);
    bool labelMatched = false;
    tools::Tokenizer tok(line, " ");
    std::vector<std::string> fields;
    tok.ToVector(fields);
    // fields = TrimCommentsFrom_(fields);
    // If not check the size of the vector and parse according
    // to the number of fields
    if (fields.size() == 1) {
      labelMatched = MatchOneFieldLabel_(fields, top);
    } else if (fields.size() == 2) {
      labelMatched = MatchTwoFieldLabels_(fields, top);
    } else if (fields.size() == 3) {
      labelMatched = MatchThreeFieldLabels_(fields, top);
    } else if (fields.size() == 4) {
      labelMatched = MatchFourFieldLabels_(fields, top);
    } else if (fields.size() != 0) {

      // See if the line is the lammps .data header/info line
      labelMatched = MatchFieldsTimeStepLabel_(fields, top);

      if (!labelMatched) {
        std::string err = "Unrecognized line in lammps .data file:\n" + line;
        throw std::runtime_error(err);
      }
    }
    getline(fl_, line);
  }
  return !fl_.eof();
}

/*****************************************************************************
 * Private Facing Methods                                                    *
 *****************************************************************************/

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::MatchOneFieldLabel_(
    std::vector<std::string> fields, Topology_T &top) {

  if (fields.at(0) == "Masses") {
    SortIntoDataGroup_("Masses");
    InitializeAtomAndBeadTypes_();
  } else if (fields.at(0) == "Atoms") {
    ReadAtoms_(top);
  } else if (fields.at(0) == "Bonds") {
    ReadBonds_(top);
  } else if (fields.at(0) == "Angles") {
    ReadAngles_(top);
  } else if (fields.at(0) == "Dihedrals") {
    ReadDihedrals_(top);
  } else if (fields.at(0) == "Impropers") {
  } else {
    return false;
  }
  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::MatchTwoFieldLabels_(
    std::vector<std::string> fields, Topology_T &top) {

  std::string label = fields.at(0) + " " + fields.at(1);

  if (fields.at(1) == "atoms") {
    ReadNumOfAtoms_(fields, top);
  } else if (fields.at(1) == "bonds") {
    ReadNumOfBonds_(fields);
  } else if (fields.at(1) == "angles") {
    ReadNumOfAngles_(fields);
  } else if (fields.at(1) == "dihedrals") {
    ReadNumOfDihedrals_(fields);
  } else if (fields.at(1) == "impropers") {
    ReadNumOfImpropers_(fields);
  } else if (label == "Pair Coeffs") {
    SortIntoDataGroup_("Pair Coeffs");
  } else if (label == "Bond Coeffs") {
    SortIntoDataGroup_("Bond Coeffs");
  } else if (label == "Angle Coeffs") {
    SortIntoDataGroup_("Angle Coeffs");
  } else if (label == "Improper Coeffs") {
    SortIntoDataGroup_("Improper Coeffs");
  } else {
    return false;
  }
  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::MatchThreeFieldLabels_(
    std::vector<std::string> fields, Topology_T &top) {
  std::string label = fields.at(1) + " " + fields.at(2);
  if (label == "atom types") {
    ReadNumTypes_(fields, "atom");
  } else if (label == "bond types") {
    ReadNumTypes_(fields, "bond");
  } else if (label == "angle types") {
    ReadNumTypes_(fields, "angle");
  } else if (label == "dihedral types") {
    ReadNumTypes_(fields, "Dihedral");
  } else if (label == "improper types") {
    ReadNumTypes_(fields, "Improper");
  } else {
    return false;
  }
  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::MatchFourFieldLabels_(
    std::vector<std::string> fields, Topology_T &top) {
  std::string label = fields.at(2) + " " + fields.at(3);
  if (label == "xlo xhi") {
    ReadBox_(fields, top);
  } else {
    return false;
  }
  return true;
}

template <class Topology_T>
bool LAMMPSDataReader<Topology_T>::MatchFieldsTimeStepLabel_(
    std::vector<std::string> fields, Topology_T &top) {
  size_t index = 0;
  for (auto field : fields) {
    if (field == "timestep" && (index + 2) < fields.size()) {
      top.setStep(stoi(fields.at(index + 2)));
      return true;
    }
    ++index;
  }
  return false;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::InitializeAtomAndBeadTypes_() {
  if (!data_.count("Masses")) {
    std::string err =
        "Masses must first be parsed before the atoms can be read.";
    throw std::runtime_error(err);
  }

  auto baseNamesMasses = determineBaseNameAssociatedWithMass_();
  auto baseNamesCount = determineAtomAndBeadCountBasedOnMass_(baseNamesMasses);

  tools::Elements elements;
  // If there is more than one atom type of the same element append a number
  // to the atom type name
  std::map<std::string, int> baseNameIndices;
  int index = 0;

  for (auto mass : data_["Masses"]) {
    // Determine the mass associated with the atom
    double mass_atom_bead = stod(mass.at(1));

    auto baseName =
        getStringGivenDoubleAndMap_(mass_atom_bead, baseNamesMasses, 0.01);

    std::string label = baseName;
    if (baseNamesCount[baseName] > 1) {
      if (baseNameIndices.count(baseName) == 0) {
        label += label + " Type 1";
        baseNameIndices[baseName] = 1;
      } else {
        baseNameIndices[baseName]++;
        label += label + " Type " + std::to_string(baseNameIndices[baseName]);
      }
    }

    std::string name_all_caps = boost::to_upper_copy<std::string>(baseName);
    std::string element = tools::topology_constants::unassigned_element;
    if (elements.isEleShort(baseName)) {
      element = baseName;
    } else if (elements.isEleFull(name_all_caps)) {
      element = elements.getEleShort(name_all_caps);
    }
    atomtypes_[index].push_back(baseName);
    atomtypes_[index].push_back(label);
    atomtypes_[index].push_back(element);
    ++index;
  }
}

template <class Topology_T>
std::map<std::string, double>
    LAMMPSDataReader<Topology_T>::determineBaseNameAssociatedWithMass_() {

  tools::Elements elements;
  std::map<std::string, double> baseNamesAndMasses;
  int bead_index_type = 1;
  for (auto mass : data_["Masses"]) {
    double mass_atom_bead = stod(mass.at(1));
    std::string beadElementName;
    if (elements.isMassAssociatedWithElement(mass_atom_bead, 0.01)) {
      beadElementName = elements.getEleShortClosestInMass(mass_atom_bead, 0.01);
    } else {
      beadElementName = "Bead" + std::to_string(bead_index_type);
      ++bead_index_type;
    }
    baseNamesAndMasses[beadElementName] = mass_atom_bead;
  }
  return baseNamesAndMasses;
}

template <class Topology_T>
std::map<std::string, int>
    LAMMPSDataReader<Topology_T>::determineAtomAndBeadCountBasedOnMass_(
        std::map<std::string, double> baseNamesAndMasses) {

  std::map<std::string, int> countSameElementOrBead;
  for (auto mass : data_["Masses"]) {
    double mass_atom_bead = stod(mass.at(1));
    auto baseName =
        getStringGivenDoubleAndMap_(mass_atom_bead, baseNamesAndMasses, 0.01);

    if (countSameElementOrBead.count(baseName) == 0) {
      countSameElementOrBead[baseName] = 1;
    } else {
      countSameElementOrBead[baseName]++;
    }
  }
  return countSameElementOrBead;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadBox_(std::vector<std::string> fields,
                                            Topology_T &top) {
  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
  m(0, 0) = stod(fields.at(1)) - stod(fields.at(0));

  for (int i = 1; i < 3; ++i) {
    std::string line;
    getline(fl_, line);
    stripLine(line);
    tools::Tokenizer tok(line, " ");
    tok.ConvertToVector(fields);
    if (fields.size() != 4) {
      throw std::runtime_error("invalid box format in the lammps data file");
    }

    double distance = stod(fields.at(1)) - stod(fields.at(0));
    formatDistance(distance);
    m(i, i) = distance;
  }
  top.setBox(m);
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::SortIntoDataGroup_(std::string tag) {
  std::string line;
  getline(fl_, line);
  getline(fl_, line);
  stripLine(line);

  std::vector<std::vector<std::string>> group;
  std::string data_elem;
  while (!line.empty()) {
    std::vector<std::string> mini_group;
    tools::Tokenizer tok(line, " ");
    std::vector<std::string> fields;
    tok.ToVector(fields);
    for (auto field : fields) {
      mini_group.push_back(field);
    }
    group.push_back(mini_group);
    getline(fl_, line);
    stripLine(line);
  }

  data_[tag] = group;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumTypes_(
    std::vector<std::string> fields, std::string type) {
  numberOfDifferentTypes_[type] = stoi(fields.at(0));
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumOfAtoms_(
    std::vector<std::string> fields, Topology_T &top) {
  numberOf_["atoms"] = stoi(fields.at(0));
  if (!topology_ && static_cast<size_t>(numberOf_["atoms"]) != top.BeadCount())
    std::runtime_error("Number of beads in topology and trajectory differ");
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumOfBonds_(
    std::vector<std::string> fields) {
  numberOf_["bonds"] = stoi(fields.at(0));
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumOfAngles_(
    std::vector<std::string> fields) {
  numberOf_["angles"] = stoi(fields.at(0));
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumOfDihedrals_(
    std::vector<std::string> fields) {
  numberOf_["dihedrals"] = stoi(fields.at(0));
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadNumOfImpropers_(
    std::vector<std::string> fields) {
  numberOf_["impropers"] = stoi(fields.at(0));
}

template <class Topology_T>
typename LAMMPSDataReader<Topology_T>::lammps_format
    LAMMPSDataReader<Topology_T>::determineDataFileFormat_(std::string line) {

  tools::Tokenizer tok(line, " ");
  std::vector<std::string> fields;
  tok.ConvertToVector(fields);
  lammps_format format;
  if (fields.size() == 5 || fields.size() == 8) {
    format = style_atomic;
  } else if (fields.size() == 6 || fields.size() == 9) {
    format = style_angle_bond_molecule;
  } else if (fields.size() == 7 || fields.size() == 10) {
    format = style_full;
  } else {
    throw std::runtime_error(
        "You have submitted a lammps data file with an "
        "unsupported format.");
  }
  return format;
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadAtoms_(Topology_T &top) {

  std::string line;
  getline(fl_, line);
  getline(fl_, line);
  stripLine(line);

  lammps_format format = determineDataFileFormat_(line);
  bool chargeRead = false;
  bool moleculeRead = false;
  if (format == style_angle_bond_molecule) moleculeRead = true;
  if (format == style_full) {
    moleculeRead = true;
    chargeRead = true;
  }

  std::map<int, std::string> sorted_file;
  int startingIndex;
  int startingIndexMolecule = 0;
  std::istringstream issfirst(line);
  issfirst >> startingIndex;
  if (moleculeRead) {
    issfirst >> startingIndexMolecule;
  }
  sorted_file[startingIndex] = line;
  getline(fl_, line);
  stripLine(line);

  int atomId = 0;
  int moleculeId = tools::topology_constants::unassigned_molecule_id;
  while (!line.empty()) {
    std::istringstream iss(line);
    iss >> atomId;
    if (moleculeRead) {
      iss >> moleculeId;
    }
    sorted_file[atomId] = line;
    getline(fl_, line);
    stripLine(line);
    if (atomId < startingIndex) startingIndex = atomId;
    if (moleculeId < startingIndexMolecule) startingIndexMolecule = moleculeId;
  }

  for (int atomIndex = startingIndex;
       static_cast<size_t>(atomIndex - startingIndex) < sorted_file.size();
       ++atomIndex) {

    int atomTypeId;
    double charge = 0;
    double x, y, z;

    std::istringstream iss(sorted_file[atomIndex]);
    iss >> atomId;
    if (moleculeRead) {
      iss >> moleculeId;
    }
    iss >> atomTypeId;
    if (chargeRead) iss >> charge;
    iss >> x;
    iss >> y;
    iss >> z;

    moleculeId = startingIndexMolecule;

    // Exclusion list assumes beads start with ids of 0
    formatId(atomId);
    formatId(atomTypeId);
    formatId(moleculeId);

    // We want to start with an index of 0 not 1
    // atomId;
    typename Topology_T::bead_t *b;
    if (topology_) {

      atomIdToIndex_[atomId] = atomIndex - startingIndex;
      atomIdToMoleculeId_[atomId] = moleculeId;
      typename Topology_T::container_t *mol;
      if (!molecules_.count(moleculeId)) {
        tools::StructureParameters params;
        params.set(tools::StructureParameter::MoleculeId, moleculeId);
        params.set(tools::StructureParameter::MoleculeType,
                   tools::topology_constants::unassigned_molecule_type);
        mol = &(top.CreateMolecule(params));
        molecules_[moleculeId] = mol;
      } else {
        mol = molecules_[moleculeId];
      }
      tools::byte_t symmetry = 1;  // spherical
      double mass =
          boost::lexical_cast<double>(data_["Masses"].at(atomTypeId).at(1));

      if (atomtypes_.count(atomTypeId) == 0) {
        std::string err =
            "Unrecognized atomTypeId, the atomtypes map "
            "may be uninitialized";
        throw std::runtime_error(err);
      }

      formatMass(mass);
      formatCharge(charge);

      std::string atom_type = atomtypes_.at(atomTypeId).at(0);
      std::string element = atomtypes_[atomTypeId].at(2);
      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::CSG_Mass, mass);
      params.set(tools::StructureParameter::CSG_Charge, charge);
      params.set(tools::StructureParameter::Element, element);
      params.set(tools::StructureParameter::BeadId, atomId);
      params.set(tools::StructureParameter::BeadType, atom_type);
      params.set(tools::StructureParameter::ResidueId,
                 tools::topology_constants::unassigned_residue_id);
      params.set(tools::StructureParameter::ResidueType,
                 tools::topology_constants::unassigned_residue_type);
      params.set(tools::StructureParameter::MoleculeId, mol->getId());
      b = &(top.CreateBead(params));
      mol->AddBead(b);
      b->setMoleculeId(mol->getId());

    } else {
      b = top.getBead(atomIndex - startingIndex);
    }

    formatDistance(x);
    formatDistance(y);
    formatDistance(z);
    Eigen::Vector3d xyz_pos(x, y, z);

    b->setPos(xyz_pos);
  }

  if (top.BeadCount() != static_cast<size_t>(numberOf_["atoms"])) {
    std::string err =
        "The number of atoms read in is not equivalent to the "
        "number of atoms indicated to exist in the lammps data file. \n"
        "Number of atoms that should exist " +
        std::to_string(numberOf_["atoms"]) +
        "\nNumber of atoms that were read in " +
        std::to_string(top.BeadCount()) + "\n";
    throw std::runtime_error(err);
  }
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadBonds_(Topology_T &top) {
  std::string line;
  getline(fl_, line);
  getline(fl_, line);
  stripLine(line);

  int bondId;
  int bondTypeId;
  int atom1Id, atom2Id;

  int bond_count = 0;
  while (!line.empty()) {

    if (topology_) {
      std::istringstream iss(line);
      iss >> bondId;
      iss >> bondTypeId;
      iss >> atom1Id;
      iss >> atom2Id;

      --atom1Id;
      --atom2Id;
      --bondTypeId;
      --bondId;

      int atom1Index = atomIdToIndex_[atom1Id];
      int atom2Index = atomIdToIndex_[atom2Id];

      typename Topology_T::bead_t *b = top.getBead(atom1Index);
      typename Topology_T::container_t *mi =
          &top.getMolecule(b->getMoleculeId());

      Interaction *ic = top.CreateInteraction(
          InteractionType::bond, "BONDS", bondId, mi->getId(),
          std::vector<int>{atom1Index, atom2Index});
      mi->AddInteraction(ic);
    }

    ++bond_count;
    getline(fl_, line);
    stripLine(line);
  }

  if (bond_count != numberOf_["bonds"]) {
    std::string err =
        "The number of bonds read in is not equivalent to the "
        "number of bonds indicated to exist in the lammps data file. \n"
        "Number of bonds that should exist " +
        std::to_string(numberOf_["bonds"]) +
        "\nNumber of bonds that were read in " + std::to_string(bond_count) +
        "\n";
    throw std::runtime_error(err);
  }
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadAngles_(Topology_T &top) {
  std::string line;
  getline(fl_, line);
  getline(fl_, line);
  stripLine(line);

  int angleId;
  int angleTypeId;
  int atom1Id, atom2Id, atom3Id;

  int angle_count = 0;

  while (!line.empty()) {

    if (topology_) {
      std::istringstream iss(line);
      iss >> angleId;
      iss >> angleTypeId;
      iss >> atom1Id;
      iss >> atom2Id;
      iss >> atom3Id;

      --angleId;
      --atom1Id;
      --atom2Id;
      --atom3Id;
      --angleTypeId;

      int atom1Index = atomIdToIndex_[atom1Id];
      int atom2Index = atomIdToIndex_[atom2Id];
      int atom3Index = atomIdToIndex_[atom3Id];

      typename Topology_T::bead_t *b = top.getBead(atom1Index);
      typename Topology_T::container_t *mi =
          &top.getMolecule(b->getMoleculeId());

      Interaction *ic = top.CreateInteraction(
          InteractionType::angle, "ANGLES", angleId, mi->getId(),
          std::vector<int>{atom1Index, atom2Index, atom3Index});
      mi->AddInteraction(ic);
    }

    ++angle_count;

    getline(fl_, line);
    stripLine(line);
  }

  if (angle_count != numberOf_["angles"]) {
    std::string err =
        "The number of angles read in is not equivalent to the "
        "number of angles indicated to exist in the lammps data file. \n"
        "Number of angles that should exist " +
        std::to_string(numberOf_["angles"]) +
        "\nNumber of angles that were read in " + std::to_string(angle_count) +
        "\n";
    throw std::runtime_error(err);
  }
}

template <class Topology_T>
void LAMMPSDataReader<Topology_T>::ReadDihedrals_(Topology_T &top) {
  std::string line;
  getline(fl_, line);
  getline(fl_, line);
  stripLine(line);

  int dihedralId;
  int dihedralTypeId;
  int atom1Id, atom2Id, atom3Id, atom4Id;

  int dihedral_count = 0;
  while (!line.empty()) {

    if (topology_) {
      std::istringstream iss(line);
      iss >> dihedralId;
      iss >> dihedralTypeId;
      iss >> atom1Id;
      iss >> atom2Id;
      iss >> atom3Id;
      iss >> atom4Id;

      --dihedralId;
      --atom1Id;
      --atom2Id;
      --atom3Id;
      --atom4Id;
      --dihedralTypeId;

      int atom1Index = atomIdToIndex_[atom1Id];
      int atom2Index = atomIdToIndex_[atom2Id];
      int atom3Index = atomIdToIndex_[atom3Id];
      int atom4Index = atomIdToIndex_[atom4Id];

      typename Topology_T::bead_t *b = top.getBead(atom1Index);
      typename Topology_T::container_t *mi =
          &top.getMolecule(b->getMoleculeId());
      Interaction *ic = top.CreateInteraction(
          InteractionType::dihedral, "DIHEDRALS", dihedralId, mi->getId(),
          std::vector<int>{atom1Index, atom2Index, atom3Index, atom4Index});
      mi->AddInteraction(ic);
    }
    ++dihedral_count;
    getline(fl_, line);
    stripLine(line);
  }

  if (dihedral_count != numberOf_["dihedrals"]) {
    std::string err =
        "The number of dihedrals read in is not equivalent to the "
        "number of dihedrals indicated to exist in the lammps data file. \n"
        "Number of dihedrals that should exist " +
        std::to_string(numberOf_["dihedrals"]) +
        "\nNumber of dihedrals that were read in " +
        std::to_string(dihedral_count) + "\n";
    throw std::runtime_error(err);
  }
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_LAMMPSDATAREADER_PRIV_H
