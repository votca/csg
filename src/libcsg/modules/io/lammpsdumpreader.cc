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

#include "lammpsdumpreader.h"
#include <boost/lexical_cast.hpp>
#include <memory>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/tools/getline.h>

namespace votca {
namespace csg {
using namespace boost;
using namespace std;
using namespace votca::tools;
bool LAMMPSDumpReader::ReadTopology(string file, CSG_Topology &top) {
  read_topology_data_ = true;
  top.Cleanup();

  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open topology file: " + file);
  _fname = file;

  NextFrame(top);

  _fl.close();

  return true;
}

bool LAMMPSDumpReader::Open(const string &file) {
  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  _fname = file;
  return true;
}

void LAMMPSDumpReader::Close() { _fl.close(); }

bool LAMMPSDumpReader::FirstFrame(CSG_Topology &top) {
  read_topology_data_ = false;
  NextFrame(top);
  return true;
}

bool LAMMPSDumpReader::NextFrame(CSG_Topology &top) {
  string line;
  getline(_fl, line);
  cout << "Reading lammps dump file" << endl;
  while (!_fl.eof()) {
    if (line.substr(0, 5) != "ITEM:")
      throw std::ios_base::failure("unexpected line in lammps file:\n" + line);
    if (line.substr(6, 8) == "TIMESTEP") {
      ReadTimestep(top, line);
    } else if (line.substr(6, 15) == "NUMBER OF ATOMS") {
      ReadNumAtoms(top, line);
    } else if (line.substr(6, 10) == "BOX BOUNDS") {
      cout << "Reading box bounds " << endl;
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
    cout << "WARNING: topology created from .dump file, masses, charges, "
            "types, residue names are wrong!\n";
  }
  return !_fl.eof();
  ;
}

void LAMMPSDumpReader::ReadTimestep(CSG_Topology &top, const string &itemline) {
  string s;
  getline(_fl, s);
  top.setStep(boost::lexical_cast<int>(s));
  cout << "Reading frame, timestep " << top.getStep() << endl;
}

void LAMMPSDumpReader::ReadBox(CSG_Topology &top, const string &itemline) {
  string s;

  matrix m;
  m.ZeroMatrix();

  for (int i = 0; i < 3; ++i) {
    getline(_fl, s);
    Tokenizer tok(s, " ");
    vector<double> v;
    tok.ConvertToVector(v);
    if (v.size() != 2) throw std::ios_base::failure("invalid box format");
    m[i][i] = v[1] - v[0];
  }
  cout << "Reading dump file box and setting " << endl;
  top.setBox(m);
}

void LAMMPSDumpReader::ReadNumAtoms(CSG_Topology &top, const string &itemline) {
  string s;
  getline(_fl, s);
  number_of_atoms_ = boost::lexical_cast<int>(s);
  if (!read_topology_data_ &&
      static_cast<size_t>(number_of_atoms_) != top.BeadCount()) {
    std::runtime_error("number of beads in topology and trajectory differ");
  }
}

void LAMMPSDumpReader::ReadAtoms(CSG_Topology &top, string itemline) {

  bool pos = false;
  bool force = false;
  bool vel = false;
  int id = -1;

  vector<string> fields;

  {
    Tokenizer tok(itemline.substr(12), " ");
    tok.ToVector(fields);
    int j = 0;
    for (Tokenizer::iterator i = tok.begin(); i != tok.end(); ++i, ++j) {
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
    string s;
    getline(_fl, s);
    if (_fl.eof())
      throw std::runtime_error(
          "Error: unexpected end of lammps file '" + _fname + "' only " +
          boost::lexical_cast<string>(i) + " atoms of " +
          boost::lexical_cast<string>(number_of_atoms_) + " read.");

    Tokenizer tok(s, " ");
    Tokenizer::iterator itok = tok.begin();
    vector<string> fields2;
    tok.ToVector(fields2);
    // Lammps starts with ids at 1 we handle ids internally at 0
    int atom_id = boost::lexical_cast<int>(fields2[id]) - 1;
    if (atom_id > number_of_atoms_) {
      throw std::runtime_error(
          "Error: found atom with id " + boost::lexical_cast<string>(atom_id) +
          " but only " + boost::lexical_cast<string>(number_of_atoms_) +
          " atoms defined in header of file '" + _fname + "'");
    }

    matrix m = top.getBox();

    unordered_map<string, double> atom_attributes_double;
    unordered_map<string, int> atom_attributes_int;
    unordered_map<string, string> atom_attributes_string;

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
        atom_attributes_double["xs"] = stod(*itok) * m[0][0];
      } else if (fields[j] == "ys") {
        atom_attributes_double["ys"] = stod(*itok) * m[1][1];
      } else if (fields[j] == "zs") {
        atom_attributes_double["zs"] = stod(*itok) * m[2][2];
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
      byte_t symmetry = 1;
      string residue_type = topology_constants::unassigned_residue_type;
      int residue_id = topology_constants::unassigned_residue_id;
      top.CreateBead(
          symmetry, atom_attributes_string["type"], atom_attributes_int["id"],
          atom_attributes_int["mol"], residue_id, residue_type,
          atom_attributes_string["element"], atom_attributes_double["mass"],
          atom_attributes_double["q"]);
    }
    cout << "Getting bead of id " << atom_attributes_int["id"] << endl;
    Bead *b = top.getBead(atom_attributes_int["id"]);
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
