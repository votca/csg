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

#ifndef VOTCA_CSG_XMLTOPOLOGYREADER_H
#define VOTCA_CSG_XMLTOPOLOGYREADER_H

#include <boost/unordered_map.hpp>
#include <stack>
#include <string>
#include <votca/csg/topologyreader.h>
#include <votca/tools/parsexml.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class BondBead {
 public:
  BondBead(std::string &line) {
    TOOLS::Tokenizer tok(line, ":");
    std::vector<std::string> tmp_vec;
    tok.ToVector(tmp_vec);
    if (tmp_vec.size() != 2)
      throw std::runtime_error("Wrong number of elements in bead: " + line);
    molecule_type_ = tmp_vec[0];
    atom_type_ = tmp_vec[1];
    molecule_type_.erase(molecule_type_.find_last_not_of(" \n\r\t") + 1);
    atom_type_.erase(atom_type_.find_last_not_of(" \n\r\t") + 1);
  }

  std::string molecule_type_;
  std::string atom_type_;
};

class XMLBead {
 public:
  XMLBead(std::string _name, std::string _type, int _residue_number,
          double _mass = 1.0, double _q = 0.0)
      : name(_name),
        type(_type),
        residue_number(_residue_number),
        mass(_mass),
        q(_q){};
  XMLBead(){};

  int pid;
  std::string name;
  std::string type;
  int residue_number;
  double mass;
  double q;
};

class XMLMolecule {
 public:
  XMLMolecule(std::string _name, int _nmols) : name(_name), nmols(_nmols) {}
  std::string name;
  int nmols;
  int pid;
  // std::vector<XMLBead> beads;
  std::map<std::string, XMLBead> name2beads;
  Molecule *mi;
};

/**
 *  Reads in an xml topology
 *
 * \todo this is a sloppy implementation using expat, is just reads attributes
 * \todo should be extended to also read beads, ...
 *
 */
class XMLTopologyReader : public TopologyReader {
 public:
  /// read a topology file
  bool ReadTopology(std::string file, CSG_Topology &top);
  ~XMLTopologyReader();

 private:
  typedef boost::unordered_multimap<std::string, XMLMolecule> MoleculesMap;

  void ReadTopolFile(std::string file);

  void ParseRoot(TOOLS::Property &el);
  void ParseMolecules(TOOLS::Property &el);
  void ParseBeadTypes(TOOLS::Property &el);
  void ParseBonded(TOOLS::Property &el);
  void ParseBox(TOOLS::Property &p);
  void ParseMolecule(TOOLS::Property &p, std::string molecule_type_, int nbeads,
                     int nmols);
  void ParseBond(TOOLS::Property &p);
  void ParseAngle(TOOLS::Property &p);
  void ParseDihedral(TOOLS::Property &p);

 private:
  TOOLS::ParseXML _parser;

  CSG_Topology *_top;
  MoleculesMap _molecules;
  int _mol_index;
  int _bead_index;

  bool _has_base_topology;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_XMLTOPOLOGYREADER_H
