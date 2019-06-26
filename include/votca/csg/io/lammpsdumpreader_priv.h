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
#ifndef VOTCA_CSG_LAMMPSDUMPREADER_PRIV_H
#define VOTCA_CSG_LAMMPSDUMPREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
int LAMMPSDumpReader<Topology_T>::formatId_(const int &id) {
  return id - 1;
}

template <class Topology_T>
double LAMMPSDumpReader<Topology_T>::formatDistance_(const double &distance) {
  return converter_.convert(distance_unit, Topology_T::units::distance_unit) *
         distance;
}

template <class Topology_T>
double LAMMPSDumpReader<Topology_T>::formatForce_(const double &force) {
  return converter_.convert(force_unit, Topology_T::units::force_unit) * force;
}

template <class Topology_T>
double LAMMPSDumpReader<Topology_T>::formatVelocity_(const double &velocity) {
  return converter_.convert(velocity_unit, Topology_T::units::velocity_unit) *
         velocity;
}

template <class Topology_T>
double LAMMPSDumpReader<Topology_T>::formatCharge_(const double &charge) {
  return converter_.convert(charge_unit, Topology_T::units::charge_unit) *
         charge;
}

template <class Topology_T>
double LAMMPSDumpReader<Topology_T>::formatMass_(const double &mass) {
  return converter_.convert(mass_unit, Topology_T::units::mass_unit) * mass;
}

template <class Topology_T>
bool LAMMPSDumpReader<Topology_T>::ReadTopology(const std::string &file,
                                                boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps dump reader read topology, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  read_topology_data_ = true;
  top.Cleanup();

  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open topology file: " + file);
  _fname = file;

  NextFrame(top_any);

  _fl.close();

  return true;
}

template <class Topology_T>
bool LAMMPSDumpReader<Topology_T>::Open(const std::string &file) {
  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  _fname = file;
  return true;
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::Close() {
  _fl.close();
}

template <class Topology_T>
bool LAMMPSDumpReader<Topology_T>::FirstFrame(boost::any top) {
  read_topology_data_ = false;
  NextFrame(top);
  return true;
}

template <class Topology_T>
bool LAMMPSDumpReader<Topology_T>::NextFrame(boost::any top_any) {

  if (warning_msg_printed_ == false) {
    std::cout
        << "WARNING: Currently the votca lammps dump reader only supports ";
    std::cout << "lammps units specified by the 'units real' ";
    std::cout << "command." << std::endl;
    std::cout << "mass: grams/mole" << std::endl;
    std::cout << "distance: angstroms" << std::endl;
    std::cout << "time: femtoseconds" << std::endl;
    std::cout << "energy: Kcal/mole" << std::endl;
    std::cout << "charge: e" << std::endl;
    std::cout << std::endl;
    warning_msg_printed_ = true;
  }

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps dump reader next frame, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  std::string line;
  getline(_fl, line);
  while (!_fl.eof()) {

    if (line.substr(0, 5) != "ITEM:")
      throw std::ios_base::failure("unexpected line in lammps file:\n" + line);
    if (line.substr(6, 8) == "TIMESTEP") {
      ReadTimestep(top, line);
    } else if (line.substr(6, 15) == "NUMBER OF ATOMS") {
      ReadNumAtoms(top, line);
    } else if (line.substr(6, 10) == "BOX BOUNDS") {
      ReadBox(top, line);
    } else if (line.substr(6, 5) == "ATOMS") {
      ReadAtoms(top, line);
      break;
    }

    else {
      throw std::ios_base::failure("unknown item lammps file : " +
                                   line.substr(6));
    }
    getline(_fl, line);
  }
  return !_fl.eof();
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadTimestep(Topology_T &top,
                                                const std::string &itemline) {
  std::string s;
  getline(_fl, s);
  boost::algorithm::trim(s);
  top.setStep(boost::lexical_cast<int>(s));
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadBox(Topology_T &top,
                                           const std::string &itemline) {
  std::string s;

  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();

  for (int i = 0; i < 3; ++i) {
    getline(_fl, s);
    tools::Tokenizer tok(s, " ");
    std::vector<double> v;
    tok.ConvertToVector(v);
    if (v.size() != 2) throw std::ios_base::failure("invalid box format");
    m(i, i) = formatDistance_(v[1] - v[0]);
  }
  top.setBox(m);
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadNumAtoms(Topology_T &top,
                                                const std::string &itemline) {
  std::string s;
  getline(_fl, s);
  boost::algorithm::trim(s);
  number_of_atoms_ = boost::lexical_cast<int>(s);
  if (!read_topology_data_ &&
      static_cast<size_t>(number_of_atoms_) != top.BeadCount()) {
    std::runtime_error("number of beads in topology and trajectory differ");
  }
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadAtoms(Topology_T &top,
                                             std::string itemline) {

  bool pos = false;
  bool force = false;
  bool vel = false;
  int id = -1;

  std::vector<std::string> fields;

  {
    tools::Tokenizer tok(itemline.substr(12), " ");
    tok.ToVector(fields);
    int j = 0;
    for (tools::Tokenizer::iterator i = tok.begin(); i != tok.end(); ++i, ++j) {
      if (*i == "x" || *i == "y" || *i == "z")
        pos = true;
      else if (*i == "xu" || *i == "yu" || *i == "zu")
        pos = true;
      else if (*i == "xs" || *i == "ys" || *i == "zs")
        pos = true;
      else if (*i == "vx" || *i == "vy" || *i == "vz")
        vel = true;
      else if (*i == "fx" || *i == "fy" || *i == "fz")
        force = true;
      else if (*i == "id")
        id = j;
    }
  }
  if (id < 0) {
    throw std::runtime_error(
        "error, id not found in any column of the atoms section");
  }

  for (int i = 0; i < number_of_atoms_; ++i) {
    std::string s;
    getline(_fl, s);
    if (_fl.eof())
      throw std::runtime_error(
          "Error: unexpected end of lammps file '" + _fname + "' only " +
          boost::lexical_cast<std::string>(i) + " atoms of " +
          boost::lexical_cast<std::string>(number_of_atoms_) + " read.");

    tools::Tokenizer tok(s, " ");
    tools::Tokenizer::iterator itok = tok.begin();
    std::vector<std::string> fields2;
    tok.ToVector(fields2);
    // Lammps starts with ids at 1 we handle ids internally at 0

    int atom_id = boost::lexical_cast<int>(fields2[id]);

    if (atom_id > number_of_atoms_) {
      throw std::runtime_error(
          "Error: found atom with id " +
          boost::lexical_cast<std::string>(atom_id) + " but only " +
          boost::lexical_cast<std::string>(number_of_atoms_) +
          " atoms defined in header of file '" + _fname + "'");
    }

    Eigen::Matrix3d m = top.getBox();

    std::unordered_map<std::string, double> atom_attributes_double;
    std::unordered_map<std::string, int> atom_attributes_int;
    std::unordered_map<std::string, std::string> atom_attributes_string;

    atom_attributes_int["id"] = formatId_(atom_id);

    for (size_t j = 0; itok != tok.end(); ++itok, ++j) {
      if (j == fields.size()) {
        throw std::runtime_error(
            "error, wrong number of columns in atoms section");
      } else if (fields[j] == "x") {
        double x_pos = stod(*itok);
        atom_attributes_double["x"] = formatDistance_(x_pos);
      } else if (fields[j] == "y") {
        double y_pos = stod(*itok);
        atom_attributes_double["y"] = formatDistance_(y_pos);
      } else if (fields[j] == "z") {
        double z_pos = stod(*itok);
        atom_attributes_double["z"] = formatDistance_(z_pos);
      } else if (fields[j] == "xu") {
        double x_pos = stod(*itok);
        atom_attributes_double["xu"] = formatDistance_(x_pos);
      } else if (fields[j] == "yu") {
        double y_pos = stod(*itok);
        atom_attributes_double["yu"] = formatDistance_(y_pos);
      } else if (fields[j] == "zu") {
        double z_pos = stod(*itok);
        atom_attributes_double["zu"] = formatDistance_(z_pos);
      } else if (fields[j] == "xs") {
        /// Scaled units x value beteen 0-1 must be scaled by multiplying by
        /// box size
        double x_pos = stod(*itok);
        x_pos *= m(0, 0);
        atom_attributes_double["xs"] = x_pos;
      } else if (fields[j] == "ys") {
        /// Scaled units y value beteen 0-1 must be scaled by multiplying by
        /// box size
        double y_pos = stod(*itok);
        y_pos *= m(1, 1);
        atom_attributes_double["ys"] = y_pos;
      } else if (fields[j] == "zs") {
        /// Scaled units z value beteen 0-1 must be scaled by multiplying by
        /// box size
        double z_pos = stod(*itok);
        z_pos *= m(2, 2);
        atom_attributes_double["zs"] = z_pos;
      } else if (fields[j] == "vx") {
        atom_attributes_double["vx"] = formatVelocity_(stod(*itok));
      } else if (fields[j] == "vy") {
        atom_attributes_double["vy"] = formatVelocity_(stod(*itok));
      } else if (fields[j] == "vz") {
        atom_attributes_double["vz"] = formatVelocity_(stod(*itok));
      } else if (fields[j] == "fx") {
        atom_attributes_double["fx"] = formatForce_(stod(*itok));
      } else if (fields[j] == "fy") {
        atom_attributes_double["fy"] = formatForce_(stod(*itok));
      } else if (fields[j] == "fz") {
        atom_attributes_double["fz"] = formatForce_(stod(*itok));
      } else if (read_topology_data_) {
        if (fields[j] == "q") {
          atom_attributes_double["q"] = formatCharge_(stod(*itok));
        } else if (fields[j] == "mol") {
          atom_attributes_int["mol"] = boost::lexical_cast<int>(*itok);
        } else if (fields[j] == "mass") {
          atom_attributes_double["mass"] = formatMass_(stod(*itok));
        } else if (fields[j] == "element") {
          atom_attributes_string["element"] = *itok;
        } else if (fields[j] == "type") {
          atom_attributes_string["type"] = *itok;
        }
      }

    }  // for (size_t j = 0; itok != tok.end(); ++itok, ++j)

    if (read_topology_data_) {
      tools::byte_t symmetry = 1;
      std::string residue_type =
          tools::topology_constants::unassigned_residue_type;
      int residue_id = tools::topology_constants::unassigned_residue_id;
      tools::StructureParameters params;
      params.set(tools::StructureParameter::Symmetry, symmetry);
      params.set(tools::StructureParameter::CSG_Mass,
                 atom_attributes_double["mass"]);
      params.set(tools::StructureParameter::CSG_Charge,
                 atom_attributes_double["q"]);
      params.set(tools::StructureParameter::BeadId, atom_attributes_int["id"]);
      params.set(tools::StructureParameter::BeadType,
                 atom_attributes_string["type"]);
      params.set(tools::StructureParameter::Element,
                 atom_attributes_string["element"]);
      params.set(tools::StructureParameter::ResidueId, residue_id);
      params.set(tools::StructureParameter::ResidueType, residue_type);
      params.set(tools::StructureParameter::MoleculeId,
                 atom_attributes_string["mol"]);
      top.CreateBead(params);
    }
    typename Topology_T::bead_t &b = top.getBead(atom_attributes_int["id"]);
    b.HasPos(pos);
    b.HasF(force);
    b.HasVel(vel);

    if (pos) {
      if (atom_attributes_double.count("x")) {
        b.Pos().x() = atom_attributes_double["x"];
        b.Pos().y() = atom_attributes_double["y"];
        b.Pos().z() = atom_attributes_double["z"];
      } else if (atom_attributes_double.count("xs")) {
        b.Pos().x() = atom_attributes_double["xs"];
        b.Pos().y() = atom_attributes_double["ys"];
        b.Pos().z() = atom_attributes_double["zs"];
      } else if (atom_attributes_double.count("xu")) {
        b.Pos().x() = atom_attributes_double["xu"];
        b.Pos().y() = atom_attributes_double["yu"];
        b.Pos().z() = atom_attributes_double["zu"];
      }
    }
    if (force) {
      b.F().x() = atom_attributes_double["fx"];
      b.F().y() = atom_attributes_double["fy"];
      b.F().z() = atom_attributes_double["fz"];
    }
    if (vel) {
      b.Vel().x() = atom_attributes_double["vx"];
      b.Vel().y() = atom_attributes_double["vy"];
      b.Vel().z() = atom_attributes_double["vz"];
    }

  }  // for (int i = 0; i < number_of_atoms_; ++i)
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_LAMMPSDUMPREADER_PRIV_H
