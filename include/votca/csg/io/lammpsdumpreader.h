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
#ifndef _VOTCA_CSG_LAMMPSDUMPREADER_H
#define _VOTCA_CSG_LAMMPSDUMPREADER_H

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

#include <votca/tools/constants.h>
#include <votca/tools/getline.h>
#include <votca/tools/structureparameters.h>

#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace votca {
namespace csg {

/**
    \brief class for reading lammps dump files

    This class provides the TrajectoryReader + Topology reader interface
    for lammps dump files

*/
template <class Topology_T>
class LAMMPSDumpReader : public TrajectoryReader, public TopologyReader {
 public:
  LAMMPSDumpReader() {}
  ~LAMMPSDumpReader() {}

  /// open a topology file
  bool ReadTopology(std::string file, boost::any top);

  /// open a trejectory file
  bool Open(const std::string &file);
  /// read in the first frame
  bool FirstFrame(boost::any top);
  /// read in the next frame
  bool NextFrame(boost::any top);

  void Close();

 private:
  void ReadTimestep(Topology_T &top, const std::string &itemline);
  void ReadBox(Topology_T &top, const std::string &itemline);
  void ReadNumAtoms(Topology_T &top, const std::string &itemline);
  void ReadAtoms(Topology_T &top, std::string itemline);

  std::ifstream _fl;
  std::string _fname;
  bool read_topology_data_ = false;
  int number_of_atoms_ = 0;
};

template <class Topology_T>
bool LAMMPSDumpReader<Topology_T>::ReadTopology(std::string file,
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

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using lammps dump reader next frame, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  std::string line;
  getline(_fl, line);
  std::cout << "Reading lammps dump file" << std::endl;
  while (!_fl.eof()) {
    if (line.substr(0, 5) != "ITEM:")
      throw std::ios_base::failure("unexpected line in lammps file:\n" + line);
    if (line.substr(6, 8) == "TIMESTEP") {
      ReadTimestep(top, line);
    } else if (line.substr(6, 15) == "NUMBER OF ATOMS") {
      ReadNumAtoms(top, line);
    } else if (line.substr(6, 10) == "BOX BOUNDS") {
      std::cout << "Reading box bounds " << std::endl;
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
  if (read_topology_data_) {
    std::cout << "WARNING: topology created from .dump file, masses, charges, "
                 "types, residue names are wrong!\n";
  }
  return !_fl.eof();
  ;
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadTimestep(Topology_T &top,
                                                const std::string &itemline) {
  std::string s;
  getline(_fl, s);
  top.setStep(boost::lexical_cast<int>(s));
  std::cout << "Reading frame, timestep " << top.getStep() << std::endl;
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
    m(i, i) = v[1] - v[0];
  }
  std::cout << "Reading dump file box and setting " << std::endl;
  top.setBox(m);
}

template <class Topology_T>
void LAMMPSDumpReader<Topology_T>::ReadNumAtoms(Topology_T &top,
                                                const std::string &itemline) {
  std::string s;
  getline(_fl, s);
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
    int atom_id = boost::lexical_cast<int>(fields2[id]) - 1;
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

    atom_attributes_int["id"] = atom_id;

    for (size_t j = 0; itok != tok.end(); ++itok, ++j) {
      if (j == fields.size()) {
        throw std::runtime_error(
            "error, wrong number of columns in atoms section");
      } else if (fields[j] == "x") {
        atom_attributes_double["x"] = stod(*itok);
      } else if (fields[j] == "y") {
        atom_attributes_double["y"] = stod(*itok);
      } else if (fields[j] == "z") {
        atom_attributes_double["z"] = stod(*itok);
      } else if (fields[j] == "xu") {
        atom_attributes_double["xu"] = stod(*itok);
      } else if (fields[j] == "yu") {
        atom_attributes_double["yu"] = stod(*itok);
      } else if (fields[j] == "zu") {
        atom_attributes_double["zu"] = stod(*itok);
      } else if (fields[j] == "xs") {
        atom_attributes_double["xs"] = stod(*itok) * m(0, 0);
      } else if (fields[j] == "ys") {
        atom_attributes_double["ys"] = stod(*itok) * m(1, 1);
      } else if (fields[j] == "zs") {
        atom_attributes_double["zs"] = stod(*itok) * m(2, 2);
      } else if (fields[j] == "vx") {
        atom_attributes_double["vx"] = stod(*itok);
      } else if (fields[j] == "vy") {
        atom_attributes_double["vy"] = stod(*itok);
      } else if (fields[j] == "vz") {
        atom_attributes_double["vz"] = stod(*itok);
      } else if (fields[j] == "fx") {
        atom_attributes_double["fx"] = stod(*itok);
      } else if (fields[j] == "fy") {
        atom_attributes_double["fy"] = stod(*itok);
      } else if (fields[j] == "fz") {
        atom_attributes_double["fz"] = stod(*itok);
      } else if (read_topology_data_) {
        if (fields[j] == "q") {
          atom_attributes_double["q"] = stod(*itok);
        } else if (fields[j] == "mol") {
          atom_attributes_int["mol"] = boost::lexical_cast<int>(*itok);
        } else if (fields[j] == "mass") {
          atom_attributes_double["mass"] = stod(*itok);
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
      params.set(tools::StructureParameter::Mass,
                 atom_attributes_double["mass"]);
      params.set(tools::StructureParameter::Charge,
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
    typename Topology_T::bead_t *b = top.getBead(atom_attributes_int["id"]);
    b->HasPos(pos);
    b->HasF(force);
    b->HasVel(vel);

    if (pos) {
      if (atom_attributes_double.count("x")) {
        b->Pos().x() = atom_attributes_double["x"];
        b->Pos().y() = atom_attributes_double["y"];
        b->Pos().z() = atom_attributes_double["z"];
      } else if (atom_attributes_double.count("xs")) {
        b->Pos().x() = atom_attributes_double["xs"];
        b->Pos().y() = atom_attributes_double["ys"];
        b->Pos().z() = atom_attributes_double["zs"];
      } else if (atom_attributes_double.count("xu")) {
        b->Pos().x() = atom_attributes_double["xu"];
        b->Pos().y() = atom_attributes_double["yu"];
        b->Pos().z() = atom_attributes_double["zu"];
      }
    }
    if (force) {
      b->F().x() = atom_attributes_double["fx"];
      b->F().y() = atom_attributes_double["fy"];
      b->F().z() = atom_attributes_double["fz"];
    }
    if (vel) {
      b->Vel().x() = atom_attributes_double["vx"];
      b->Vel().y() = atom_attributes_double["vy"];
      b->Vel().z() = atom_attributes_double["vz"];
    }

  }  // for (int i = 0; i < number_of_atoms_; ++i)
}

}  // namespace csg
}  // namespace votca
#endif  // _VOTCA_CSG_LAMMPSDUMPREADER_H