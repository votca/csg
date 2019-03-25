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

#include "xmltopologyreader.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <votca/csg/bead.h>
#include <votca/csg/csgtopology.h>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

using namespace votca::tools;

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;

bool XMLTopologyReader::ReadTopology(string filename, CSG_Topology &top) {
  _top = &top;

  Property options;
  load_property_from_xml(options, filename);
  ParseRoot(options.get("topology"));

  _top->RebuildExclusions();
  return true;
}

void XMLTopologyReader::ReadTopolFile(string file) {
  TopologyReader *reader;
  reader = TopReaderFactory().Create(file);
  if (!reader) throw runtime_error(file + ": unknown topology format");

  reader->ReadTopology(file, *_top);

  delete reader;
  // Clean XML molecules and beads.
}

void XMLTopologyReader::ParseRoot(Property &property) {
  _has_base_topology = false;
  if (property.hasAttribute("base")) {
    ReadTopolFile(property.getAttribute<string>("base"));
    _has_base_topology = true;
  }

  // Iterate over keys at first level.
  for (Property::iterator it = property.begin(); it != property.end(); ++it) {
    if (it->name() == "h5md_particle_group") {
      _top->setParticleGroup(it->getAttribute<string>("name"));
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
      throw runtime_error("unknown tag: topology." + it->name());
    }
  }
}

void XMLTopologyReader::ParseBox(Property &p) {
  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
  m(0, 0) = p.getAttribute<double>("xx");
  m(1, 1) = p.getAttribute<double>("yy");
  m(2, 2) = p.getAttribute<double>("zz");
  _top->setBox(m);
}

void XMLTopologyReader::ParseMolecules(Property &p) {
  for (Property::iterator it = p.begin(); it != p.end(); ++it) {
    if (it->name() == "clear") {
      _top->ClearMoleculeList();
    } else if (it->name() == "rename") {
      string molecule_type_ = it->getAttribute<string>("name");
      string range = it->getAttribute<string>("range");
      _top->RenameMoleculesType(range, molecule_type_);
    } else if (it->name() == "define" || it->name() == "molecule") {
      string molecule_type_ = it->getAttribute<string>("name");
      int first = 0;
      if (it->name() == "define") first = it->getAttribute<int>("first");
      int nbeads = it->getAttribute<int>("nbeads");
      int nmols = it->getAttribute<int>("nmols");
      if (it->name() == "define" && first < 1)
        throw std::runtime_error(
            "Attribute first is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it->getAttribute<string>("first")));
      if (nbeads < 1)
        throw std::runtime_error(
            "Attribute nbeads is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it->getAttribute<string>("nbeads")));
      if (nmols < 1)
        throw std::runtime_error(
            "Attribute nmols is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it->getAttribute<string>("nmols")));
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

void XMLTopologyReader::ParseMolecule(Property &p, string molecule_type_,
                                      int nbeads, int nmols) {
  vector<XMLBead> xmlBeads;
  vector<int> xmlResidues;
  for (Property::iterator it = p.begin(); it != p.end(); ++it) {
    if (it->name() == "bead") {
      string atom_type_ = it->getAttribute<string>("name");
      string attype = it->getAttribute<string>("type");
      double atmass, atq;
      int resid;
      try {
        atmass = it->getAttribute<double>("mass");
      } catch (runtime_error &) {
        atmass = 1.0;
      }
      try {
        atq = it->getAttribute<double>("q");
      } catch (runtime_error &) {
        atq = 0.0;
      }
      try {
        resid = it->getAttribute<int>("resid");
        if (resid <= 0) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule has to "
              "be greater than zero");
        }
      } catch (runtime_error &) {
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

  Elements elements;
  for (int mn = 0; mn < nmols; mn++) {
    Molecule *mi = _top->CreateMolecule(_top->MoleculeCount(), molecule_type_);
    XMLMolecule xmlMolecule = XMLMolecule(molecule_type_, nmols);
    xmlMolecule.pid = mi->getId();
    xmlMolecule.mi = mi;
    _molecules.insert(make_pair(molecule_type_, xmlMolecule));
    unordered_map<string, int> residuename_residuenumber;
    for (XMLBead &xml_bead : xmlBeads) {
      XMLBead &b = xml_bead;

      if (residuename_residuenumber.count(molecule_type_) == 0) {
        residuename_residuenumber[molecule_type_] = 1;
      } else {
        ++residuename_residuenumber[molecule_type_];
      }
      byte_t symmetry = 1;

      string element = topology_constants::unassigned_element;
      string name_all_caps = boost::algorithm::to_upper_copy<string>(b.name);
      if (elements.isEleShort(b.name)) {
        element = b.name;
      } else if (elements.isEleFull(name_all_caps)) {
        element = elements.getEleShort(name_all_caps);
      }
      Bead *bead = _top->CreateBead(symmetry, b.type, _top->BeadCount(),
                                    xmlMolecule.pid, b.residue_number,
                                    topology_constants::unassigned_residue_type,
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

void XMLTopologyReader::ParseBeadTypes(Property &el) {
  for (Property::iterator it = el.begin(); it != el.end(); ++it) {
    if (it->name() == "rename") {
      string bead_type = it->getAttribute<string>("name");
      string bead_new_type = it->getAttribute<string>("newname");
      if (bead_type == "" || bead_new_type == "") {
        throw runtime_error("invalid rename tag, name or newname are empty.");
      }
      _top->RenameBeadsType(bead_type, bead_new_type);
    } else if (it->name() == "mass") {
      string bead_type = it->getAttribute<string>("name");
      double value = it->getAttribute<double>("value");
      _top->setBeadOfGivenTypeToNewMass(bead_type, value);
    } else {
      throw std::runtime_error("Wrong element under beadtypes: " + it->name());
    }
  }
}

void XMLTopologyReader::ParseBonded(Property &el) {
  for (Property::iterator it = el.begin(); it != el.end(); ++it) {
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

void XMLTopologyReader::ParseBond(Property &p) {
  string interaction_group = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  Tokenizer tok(beads, " \n\t");
  vector<string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 2 == 1)
    throw runtime_error("Wrong number of beads in bond: " + interaction_group);
  Interaction *ic = NULL;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    if (b1.molecule_type_ == b2.molecule_type_) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molecule_type_);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        ic = _top->CreateInteraction(InteractionType::bond, interaction_group,
                                     bond_index, xmlMolecule.pid,
                                     vector<int>{xmlBead1.pid, xmlBead2.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

void XMLTopologyReader::ParseAngle(Property &p) {
  string interaction_group = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  Tokenizer tok(beads, " \n\t");
  vector<string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 3 == 1)
    throw runtime_error("Wrong number of beads in angle: " + interaction_group);
  Interaction *ic = NULL;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
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
        XMLMolecule &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        XMLBead &xmlBead3 = xmlMolecule.name2beads[b3.atom_type_];
        ic = _top->CreateInteraction(
            InteractionType::angle, interaction_group, bond_index,
            xmlMolecule.pid,
            vector<int>{xmlBead1.pid, xmlBead2.pid, xmlBead3.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}
void XMLTopologyReader::ParseDihedral(Property &p) {
  string interaction_group = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  Tokenizer tok(beads, " \n\t");
  vector<string> bead_list;
  tok.ToVector(bead_list);
  if (bead_list.size() % 4 == 1)
    throw runtime_error("Wrong number of beads in dihedral: " +
                        interaction_group);
  Interaction *ic = NULL;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int bond_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
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
        XMLMolecule &xmlMolecule = itm->second;
        XMLBead &xmlBead1 = xmlMolecule.name2beads[b1.atom_type_];
        XMLBead &xmlBead2 = xmlMolecule.name2beads[b2.atom_type_];
        XMLBead &xmlBead3 = xmlMolecule.name2beads[b3.atom_type_];
        XMLBead &xmlBead4 = xmlMolecule.name2beads[b4.atom_type_];
        ic = _top->CreateInteraction(InteractionType::dihedral,
                                     interaction_group, bond_index,
                                     xmlMolecule.pid,
                                     vector<int>{xmlBead1.pid, xmlBead2.pid,
                                                 xmlBead3.pid, xmlBead4.pid});
        xmlMolecule.mi->AddInteraction(ic);
        bond_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

XMLTopologyReader::~XMLTopologyReader() {
  // Clean _molecules map
  /*  for (MoleculesMap::iterator it = _molecules.begin(); it !=
    _molecules.end();
         ++it) {
      XMLMolecule &xmlMolecule = it->second;
      for (vector<XMLBead *>::iterator itb = xmlMolecule->beads.begin();
           itb != xmlMolecule->beads.end(); ++itb) {
        delete (*itb);
      }
      xmlMolecule->beads.clear();
      delete xmlMolecule;
    }*/
}

}  // namespace csg
}  // namespace votca
