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

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

#include <fstream>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include "../topologyreader.h"
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/parsexml.h>

namespace votca {
namespace csg {

class BondBead {
 public:
  BondBead(std::string &line) {
    tools::Tokenizer tok(line, ":");
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

template <class Molecule_T>
class XMLMolecule {
 public:
  XMLMolecule(std::string _name, int _nmols) : name(_name), nmols(_nmols) {}
  std::string name;
  int nmols;
  int pid;
  // std::vector<XMLBead> beads;
  std::map<std::string, XMLBead> name2beads;
  Molecule_T *mi;
};

/**
 *  Reads in an xml topology
 *
 * \todo this is a sloppy implementation using expat, is just reads attributes
 * \todo should be extended to also read beads, ...
 *
 */
template <class Bead_T, class Molecule_T, class Topology_T>
class XMLTopologyReader : public TopologyReader {
 public:
  /// read a topology file
  bool ReadTopology(std::string file, void *top);
  ~XMLTopologyReader();

 private:
  typedef boost::unordered_multimap<std::string, XMLMolecule<Molecule_T>>
      MoleculesMap;

  void ReadTopolFile(std::string file);

  void ParseRoot(tools::Property &el);
  void ParseMolecules(tools::Property &el);
  void ParseBeadTypes(tools::Property &el);
  void ParseBonded(tools::Property &el);
  void ParseBox(tools::Property &p);
  void ParseMolecule(tools::Property &p, std::string molecule_type_, int nbeads,
                     int nmols);
  void ParseBond(tools::Property &p);
  void ParseAngle(tools::Property &p);
  void ParseDihedral(tools::Property &p);

 private:
  tools::ParseXML _parser;

  Topology_T *_top;
  MoleculesMap _molecules;
  int _mol_index;
  int _bead_index;

  bool _has_base_topology;
};

template <class Bead_T, class Molecule_T, class Topology_T>
bool XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ReadTopology(
    std::string filename, void *top) {
  _top = static_cast<Topology_T *>(top);

  tools::Property options;
  load_property_from_xml(options, filename);
  ParseRoot(options.get("topology"));

  _top->RebuildExclusions();
  return true;
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ReadTopolFile(
    std::string file) {
  TopologyReader *reader;
  reader = TopReaderFactory().Create(file);
  if (!reader) throw std::runtime_error(file + ": unknown topology format");

  reader->ReadTopology(file, *_top);

  delete reader;
  // Clean XML molecules and beads.
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseRoot(
    tools::Property &property) {
  _has_base_topology = false;
  if (property.hasAttribute("base")) {
    ReadTopolFile(property.getAttribute<std::string>("base"));
    _has_base_topology = true;
  }

  // Iterate over keys at first level.
  for (tools::Property::iterator it = property.begin(); it != property.end();
       ++it) {
    if (it->name() == "h5md_particle_group") {
      _top->setParticleGroup(it->getAttribute<std::string>("name"));
    } else if (it->name() == "molecules") {
      _mol_index = 1;
      _bead_index = 0;
      ParseMolecules(*it);
    } else if (it->name() == "bonded") {
      ParseBonded(*it);
    } else if (it->name() == "box") {
      ParseBox(*it);
    } else if (it->name() == "beadtypes") {
      ParseBeadTypes(*it);
    } else {
      throw std::runtime_error("unknown tag: topology." + it->name());
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseBox(
    tools::Property &p) {
  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
  m(0, 0) = p.getAttribute<double>("xx");
  m(1, 1) = p.getAttribute<double>("yy");
  m(2, 2) = p.getAttribute<double>("zz");
  _top->setBox(m);
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseMolecules(
    tools::Property &p) {
  for (tools::Property::iterator it = p.begin(); it != p.end(); ++it) {
    if (it->name() == "clear") {
      _top->ClearMoleculeList();
    } else if (it->name() == "rename") {
      std::string molecule_type_ = it->getAttribute<std::string>("name");
      std::string range = it->getAttribute<std::string>("range");
      _top->RenameMoleculesType(range, molecule_type_);
    } else if (it->name() == "define" || it->name() == "molecule") {
      std::string molecule_type_ = it->getAttribute<std::string>("name");
      int first = 0;
      if (it->name() == "define") first = it->getAttribute<int>("first");
      int nbeads = it->getAttribute<int>("nbeads");
      int nmols = it->getAttribute<int>("nmols");
      if (it->name() == "define" && first < 1)
        throw std::runtime_error(
            "Attribute first is suppose to be > 0, but found " +
            boost::lexical_cast<std::string>(
                it->getAttribute<std::string>("first")));
      if (nbeads < 1)
        throw std::runtime_error(
            "Attribute nbeads is suppose to be > 0, but found " +
            boost::lexical_cast<std::string>(
                it->getAttribute<std::string>("nbeads")));
      if (nmols < 1)
        throw std::runtime_error(
            "Attribute nmols is suppose to be > 0, but found " +
            boost::lexical_cast<std::string>(
                it->getAttribute<std::string>("nmols")));
      if (it->name() == "define") {
      } else {
        if (_has_base_topology)
          throw std::runtime_error(
              "The defined list of beads only works for pure xml topology, "
              "without 'base' attribute.");
        ParseMolecule(*it, molecule_type_, nbeads, nmols);
      }
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseMolecule(
    tools::Property &p, std::string molecule_type_, int nbeads, int nmols) {
  std::vector<XMLBead> xmlBeads;
  std::vector<int> xmlResidues;
  for (tools::Property::iterator it = p.begin(); it != p.end(); ++it) {
    if (it->name() == "bead") {
      std::string atom_type_ = it->getAttribute<std::string>("name");
      std::string attype = it->getAttribute<std::string>("type");
      double atmass, atq;
      int resid;
      try {
        atmass = it->getAttribute<double>("mass");
      } catch (std::runtime_error &) {
        atmass = 1.0;
      }
      try {
        atq = it->getAttribute<double>("q");
      } catch (std::runtime_error &) {
        atq = 0.0;
      }
      try {
        resid = it->getAttribute<int>("resid");
        if (resid <= 0) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule has to "
              "be greater than zero");
        }
      } catch (std::runtime_error &) {
        resid = -1;
      }
      if (!xmlResidues.empty()) {
        if (xmlResidues.back() != resid && xmlResidues.back() != resid - 1) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule does not "
              "increase monotonically in steps of 1");
        }
        if (xmlResidues.back() != -1 && resid == -1) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule has to "
              "be declared for all beads or for none");
        }
      }
      xmlResidues.push_back(resid);
      XMLBead xmlBead = XMLBead(atom_type_, attype, resid, atmass, atq);
      xmlBeads.push_back(xmlBead);
    } else {
      throw std::runtime_error(
          "Wrong element under topology.molecules.molecule: " + it->name());
    }
  }
  if (xmlResidues.size() != xmlBeads.size())
    throw std::runtime_error(
        "Number of elements in bead-vector and residue-vector are not "
        "identical");

  tools::Elements elements;
  for (int mn = 0; mn < nmols; mn++) {
    Molecule_T *mi =
        _top->CreateMolecule(_top->MoleculeCount(), molecule_type_);
    XMLMolecule<Molecule_T> xmlMolecule =
        XMLMolecule<Molecule_T>(molecule_type_, nmols);
    xmlMolecule.pid = mi->getId();
    xmlMolecule.mi = mi;
    _molecules.insert(make_pair(molecule_type_, xmlMolecule));
    std::unordered_map<std::string, int> residuename_residuenumber;
    for (XMLBead &xml_bead : xmlBeads) {
      XMLBead &b = xml_bead;

      if (residuename_residuenumber.count(molecule_type_) == 0) {
        residuename_residuenumber[molecule_type_] = 1;
      } else {
        ++residuename_residuenumber[molecule_type_];
      }
      tools::byte_t symmetry = 1;

      std::string element = tools::topology_constants::unassigned_element;
      std::string name_all_caps =
          boost::algorithm::to_upper_copy<std::string>(b.name);
      if (elements.isEleShort(b.name)) {
        element = b.name;
      } else if (elements.isEleFull(name_all_caps)) {
        element = elements.getEleShort(name_all_caps);
      }
      Bead_T *bead = _top->CreateBead(
          symmetry, b.type, _top->BeadCount(), xmlMolecule.pid,
          b.residue_number, tools::topology_constants::unassigned_residue_type,
          element, b.mass, b.q);

      bead->setMoleculeId(_mol_index);
      mi->AddBead(bead);

      // Data for bonded terms.
      XMLBead b_rep = XMLBead(b);
      b_rep.pid = _bead_index;
      if (xmlMolecule.name2beads.count(b.name) != 0)
        throw std::runtime_error("Atom " + b.name + " in molecule " +
                                 molecule_type_ + " already exists.");
      xmlMolecule.name2beads.insert(make_pair(b.name, b_rep));
      _bead_index++;
    }
  }
  _mol_index++;

  // clean up
  /*  for (std::vector<XMLBead *>::iterator itb = xmlBeads.begin();
         itb != xmlBeads.end(); ++itb) {
      delete (*itb);
    }*/
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseBeadTypes(
    tools::Property &el) {
  for (tools::Property::iterator it = el.begin(); it != el.end(); ++it) {
    if (it->name() == "rename") {
      std::string bead_type = it->getAttribute<std::string>("name");
      std::string bead_new_type = it->getAttribute<std::string>("newname");
      if (bead_type == "" || bead_new_type == "") {
        throw std::runtime_error(
            "invalid rename tag, name or newname are empty.");
      }
      _top->RenameBeadsType(bead_type, bead_new_type);
    } else if (it->name() == "mass") {
      std::string bead_type = it->getAttribute<std::string>("name");
      double value = it->getAttribute<double>("value");
      _top->setBeadOfGivenTypeToNewMass(bead_type, value);
    } else {
      throw std::runtime_error("Wrong element under beadtypes: " + it->name());
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseBonded(
    tools::Property &el) {
  for (tools::Property::iterator it = el.begin(); it != el.end(); ++it) {
    if (it->name() == "bond") {
      ParseBond(*it);
    } else if (it->name() == "angle") {
      ParseAngle(*it);
    } else if (it->name() == "dihedral") {
      ParseDihedral(*it);
    } else {
      throw std::runtime_error("Wrong element under bonded: " + it->name());
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseBond(
    tools::Property &p) {
  std::string interaction_group = p.get("name").as<std::string>();
  std::string beads = p.get("beads").as<std::string>();
  tools::Tokenizer tok(beads, " \n\t");
  std::vector<std::string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 2 == 1)
    throw std::runtime_error("Wrong number of beads in bond: " +
                             interaction_group);
  Interaction *ic = NULL;
  typedef std::pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (std::vector<std::string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    if (b1.molecule_type_ == b2.molecule_type_) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molecule_type_);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule<Molecule_T> &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        ic = _top->CreateInteraction(
            InteractionType::bond, interaction_group, bond_index,
            xmlMolecule.pid, std::vector<int>{xmlBead1.pid, xmlBead2.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseAngle(
    tools::Property &p) {
  std::string interaction_group = p.get("name").as<std::string>();
  std::string beads = p.get("beads").as<std::string>();
  tools::Tokenizer tok(beads, " \n\t");
  std::vector<std::string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 3 == 1)
    throw std::runtime_error("Wrong number of beads in angle: " +
                             interaction_group);
  Interaction *ic = NULL;
  typedef std::pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (std::vector<std::string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    if ((b1.molecule_type_ == b2.molecule_type_) &&
        (b2.molecule_type_ == b3.molecule_type_)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molecule_type_);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule<Molecule_T> &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        XMLBead &xmlBead3 = xmlMolecule.name2beads[b3.atom_type_];
        ic = _top->CreateInteraction(
            InteractionType::angle, interaction_group, bond_index,
            xmlMolecule.pid,
            std::vector<int>{xmlBead1.pid, xmlBead2.pid, xmlBead3.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}
template <class Bead_T, class Molecule_T, class Topology_T>
void XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::ParseDihedral(
    tools::Property &p) {
  std::string interaction_group = p.get("name").as<std::string>();
  std::string beads = p.get("beads").as<std::string>();
  tools::Tokenizer tok(beads, " \n\t");
  std::vector<std::string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 4 == 1)
    throw std::runtime_error("Wrong number of beads in dihedral: " +
                             interaction_group);
  Interaction *ic = NULL;
  typedef std::pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (std::vector<std::string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    BondBead b4(*(it++));
    if ((b1.molecule_type_ == b2.molecule_type_) &&
        (b3.molecule_type_ == b4.molecule_type_) &&
        (b1.molecule_type_ == b4.molecule_type_)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molecule_type_);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule<Molecule_T> &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        XMLBead &xmlBead3 = xmlMolecule.name2beads[b3.atom_type_];
        XMLBead &xmlBead4 = xmlMolecule.name2beads[b4.atom_type_];
        ic = _top->CreateInteraction(
            InteractionType::dihedral, interaction_group, bond_index,
            xmlMolecule.pid,
            std::vector<int>{xmlBead1.pid, xmlBead2.pid, xmlBead3.pid,
                             xmlBead4.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
XMLTopologyReader<Bead_T, Molecule_T, Topology_T>::~XMLTopologyReader() {
  // Clean _molecules map
  /*  for (MoleculesMap::iterator it = _molecules.begin(); it !=
    _molecules.end();
         ++it) {
      XMLMolecule &xmlMolecule = it->second;
      for (std::vector<XMLBead *>::iterator itb = xmlMolecule->beads.begin();
           itb != xmlMolecule->beads.end(); ++itb) {
        delete (*itb);
      }
      xmlMolecule->beads.clear();
      delete xmlMolecule;
    }*/
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_XMLTOPOLOGYREADER_H
