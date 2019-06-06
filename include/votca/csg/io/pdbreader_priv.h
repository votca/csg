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
#ifndef VOTCA_CSG_PDBREADER_PRIV_H
#define VOTCA_CSG_PDBREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
void PDBReader<Topology_T>::formatId_(int &id) {
  --id;
}

template <class Topology_T>
bool PDBReader<Topology_T>::Open(const std::string &file) {
  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open pdb trajectory file: " + file);
  return true;
}

template <class Topology_T>
void PDBReader<Topology_T>::Close() {
  _fl.close();
}

template <class Topology_T>
void PDBReader<Topology_T>::formatElement_(std::string &element_symbol,
                                           const std::string &atom_type) {
  if (!elements_.isEleShort(element_symbol)) {
    if (elements_.isEleShort(atom_type)) {
      element_symbol = atom_type;
    } else {
      element_symbol = tools::topology_constants::unassigned_element;
    }
  }
}

template <class Topology_T>
void PDBReader<Topology_T>::formatDistance_(double &distance) {
  distance *=
      converter_.convert(distance_unit, Topology_T::units::distance_unit);
}

template <class Topology_T>
bool PDBReader<Topology_T>::ReadTopology(const std::string &file,
                                         boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using pdb reader read topology, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);

  _topology = true;
  top.Cleanup();

  _fl.open(file.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open topology file: " + file);

  NextFrame(top_any);
  _fl.close();

  return true;
}

template <class Topology_T>
bool PDBReader<Topology_T>::FirstFrame(boost::any top) {
  _topology = false;
  return NextFrame(top);
}

template <class Topology_T>
bool PDBReader<Topology_T>::NextFrame(boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using pdb reader next frame, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
  std::string line;
  // Two column vector for storing all bonds
  // 1 - id of first atom
  // 2 - id of second atom
  std::vector<std::vector<int>> bond_pairs;
  ////////////////////////////////////////////////////////////////////////////////
  // Read in information from .pdb file
  ////////////////////////////////////////////////////////////////////////////////
  size_t bead_count = 0;
  // Can only support a pdb file with a single step per file at the moment
  bool step_set = false;
  while (std::getline(_fl, line)) {
    if (tools::wildcmp("MODEL*", line.c_str())) {
      // The model number is the same as the step number
      // Internallu step numbers start at 0 however pdb files store model
      // numbers starting at 1
      // 11 - 14    Model number/Step Number
      std::string model_num = std::string(line, (11 - 1), 4);
      boost::algorithm::trim(model_num);
      int step = stoi(model_num) - 1;
      top.setStep(step);
      if (step_set) {
        throw std::runtime_error(
            "Current pdbreader does not support reading more than a single "
            "model from a pdb file at a time.");
      }
      step_set = true;
    } else if (tools::wildcmp("CRYST1*", line.c_str())) {
      std::string a, b, c, alpha, beta, gamma = "";
      try {
        // 1 -  6       Record name    "CRYST1"
        a = std::string(line, (7 - 1), 9);
        // 7 - 15       Real(9.3)      a           (Angstroms)
        b = std::string(line, (16 - 1), 9);
        // 16 - 24       Real(9.3)     b           (Angstroms)
        c = std::string(line, (25 - 1), 9);
        // 25 - 33       Real(9.3)     c           (Angstroms)
        alpha = std::string(line, (34 - 1), 7);
        // 34 - 40       Real(7.2)     alpha (degrees)
        beta = std::string(line, (41 - 1), 7);
        // 41 - 47       Real(7.2)     beta        (degrees)
        gamma = std::string(line, (48 - 1), 7);
        // 48 - 54       Real(7.2)     gamma (degrees)
        // 56 - 66       LString       Space group
        // 67 - 70       Integer       Z value
      } catch (std::out_of_range &err) {
        throw std::runtime_error("Misformated pdb file in CRYST1 line");
      }
      boost::algorithm::trim(a);
      boost::algorithm::trim(b);
      boost::algorithm::trim(c);
      boost::algorithm::trim(alpha);
      boost::algorithm::trim(beta);
      boost::algorithm::trim(gamma);
      if ((!tools::wildcmp("90*", alpha.c_str())) ||
          (!tools::wildcmp("90*", alpha.c_str())) ||
          (!tools::wildcmp("90*", alpha.c_str()))) {
        throw std::runtime_error(
            "Non cubical box in pdb file not implemented, yet!");
      }
      double aa = stod(a);
      double bb = stod(b);
      double cc = stod(c);
      formatDistance_(aa);
      formatDistance_(bb);
      formatDistance_(cc);

      Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
      box.diagonal() = Eigen::Vector3d(aa, bb, cc);
      top.setBox(box);
    }
    // Only read the CONECT keyword if the topology is set too true
    if (_topology && tools::wildcmp("CONECT*", line.c_str())) {
      std::vector<std::string> bonded_atms;
      std::string atm1;
      // Keep track of the number of bonds
      int num_bonded_atoms = 0;
      try {
        // If the CONECT keyword is found then there must be at least
        // two atom identifiers, more than that is optional.

        std::stringstream ss;
        ss << std::string(line.substr(6));
        // 1 -  6       Record name    "CONECT"
        // 11 -  7       Real(5)        atm1           (ID)
        ss >> atm1;
        std::string temp_atm;
        num_bonded_atoms = 1;
        while (ss >> temp_atm) {
          boost::algorithm::trim(temp_atm);
          bonded_atms.push_back(temp_atm);
          ++num_bonded_atoms;
        }
      } catch (std::out_of_range &err) {
        throw std::runtime_error("Misformated pdb file in CONECT line\n" +
                                 line);
      }

      boost::algorithm::trim(atm1);
      // Atom ids are stored internally starting at 0 but are stored in pdb
      // files starting at 1
      int at1 = boost::lexical_cast<int>(atm1);
      formatId_(at1);

      size_t index = 0;
      while (index < bonded_atms.size()) {

        // Atom ids are stored internally starting at 0 but are stored in pdb
        // files starting at 1
        int at2 = boost::lexical_cast<int>(bonded_atms.at(index)) - 1;
        std::vector<int> row = {at1, at2};
        // Because every bond will be counted twice in a .pdb file
        // we will only add bonds where the id (atm1) is less than the
        // bonded_atm
        if (at1 < at2) {
          std::cout << "Adding bond pair " << row.at(0) << " " << row.at(1)
                    << std::endl;
          bond_pairs.push_back(row);
        }
        ++index;
      }
    }

    if (tools::wildcmp("ATOM*", line.c_str()) ||
        tools::wildcmp("HETATM*", line.c_str())) {

      // according to PDB format
      std::string atom_id_pdb, x, y, z, residue_id_pdb_str, residue_type,
          atom_type;
      std::string charge_str, element_symbol;
      // std::string atom_id_pdb;
      try {
        // Some pdb don't include all this, read only what we really need
        // leave this here in case we need more later

        // str       ,  "ATOM", "HETATM"
        // std::string recType    (line,( 1-1),6);
        // int       , Atom serial number
        atom_id_pdb = std::string(line, (7 - 1), 6);
        // str       , Atom name
        atom_type = std::string(line, (13 - 1), 4);
        // char      , Alternate location indicator
        // std::string atAltLoc   (line,(17-1),1);
        // str       , Residue name
        residue_type = std::string(line, (18 - 1), 3);
        // char      , Chain identifier
        // std::string chainID    (line,(22-1),1);
        // int       , Residue sequence number
        residue_id_pdb_str = std::string(line, (23 - 1), 4);
        // char      , Code for insertion of res
        // std::string atICode    (line,(27-1),1);
        // float 8.3 , x
        x = std::string(line, (31 - 1), 8);
        // float 8.3 , y
        y = std::string(line, (39 - 1), 8);
        // float 8.3 , z
        z = std::string(line, (47 - 1), 8);
        // float 6.2 , Occupancy
        // std::string atOccup    (line,(55-1),6);
        // float 6.2 , Temperature factor
        // std::string atTFactor  (line,(61-1),6);
        // str       , Segment identifier
        // std::string segID      (line,(73-1),4);
        // str       , Element symbol
        element_symbol = std::string(line, (77 - 1), 2);
        // str       , Charge on the atom
        charge_str = std::string(line, (79 - 1), 2);
      } catch (std::out_of_range &err) {
        std::string err_msg =
            "Misformated pdb file in atom line # " +
            boost::lexical_cast<std::string>(bead_count) +
            "\n the correct pdb file format requires 80 "
            "characters in width. Furthermore, " +
            "\n to read the topology in from a .pdb file the "
            "following attributes must be " +
            "\n specified:                                        "
            "                        " +
            "\n Atom Name, Residue Name, Residue Number, x, y, z, "
            "charge (optional)     \n";
        throw std::runtime_error(err_msg);
      }
      boost::algorithm::trim(atom_id_pdb);
      boost::algorithm::trim(atom_type);
      boost::algorithm::trim(residue_type);
      boost::algorithm::trim(residue_id_pdb_str);
      boost::algorithm::trim(x);
      boost::algorithm::trim(y);
      boost::algorithm::trim(z);
      boost::algorithm::trim(element_symbol);
      boost::algorithm::trim(charge_str);

      // Atom number is stored internally starting at 0 but .pdb files start at
      // ids of 1
      int atom_number = boost::lexical_cast<int>(atom_id_pdb);
      ++bead_count;

      formatElement_(element_symbol, atom_type);
      formatId_(atom_number);

      // Only read the CONECT keyword if the topology is set too true
      if (_topology) {
        int residue_id;
        try {
          residue_id = boost::lexical_cast<int>(residue_id_pdb_str);
        } catch (boost::bad_lexical_cast &) {
          throw std::runtime_error(
              "Cannot convert residue_id_pdb='" + residue_id_pdb_str +
              "' to int, that usallly means: misformated pdb file");
        }
        if (residue_id < 1)
          throw std::runtime_error(
              "Misformated pdb file, residue_id has to be > 0");

        // Determine if the charge_str has been provided in the .pdb file or if
        // we will be assuming it is 0
        double charge = 0;
        if (charge_str != "") {
          charge = boost::lexical_cast<double>(charge_str);
        }

        formatId_(residue_id);
        // CreateBead takes 6 parameters in the following order
        // 1 - symmetry of the bead (1-indicates sphere, 3-indicates
        // ellipsoidal)
        // 2 - name of the bead     (string)
        // 3 - bead type            (BeadType *)
        // 4 - residue number       (int)
        // 5 - mass                 (double)
        // 6 - charge               (double)
        //
        // res -1 as internal number starts with 0
        tools::byte_t symmetry = 1;

        tools::StructureParameters params;
        params.set(tools::StructureParameter::Symmetry, symmetry);
        params.set(tools::StructureParameter::CSG_Mass,
                   elements_.getMass(element_symbol));
        params.set(tools::StructureParameter::CSG_Charge, charge);
        params.set(tools::StructureParameter::Element, element_symbol);
        params.set(tools::StructureParameter::BeadId, atom_number);
        params.set(tools::StructureParameter::BeadType, atom_type);
        params.set(tools::StructureParameter::MoleculeId,
                   tools::topology_constants::unassigned_molecule_id);
        params.set(tools::StructureParameter::ResidueId, residue_id);
        params.set(tools::StructureParameter::ResidueType, residue_type);
        std::cout << "Creating Bead " << atom_number << std::endl;
        top.CreateBead(params);
      }
      typename Topology_T::bead_t &b = top.getBead(atom_number);
      // convert to nm from A
      double x_pos = stod(x);
      double y_pos = stod(y);
      double z_pos = stod(z);

      formatDistance_(x_pos);
      formatDistance_(y_pos);
      formatDistance_(z_pos);
      b.setPos(Eigen::Vector3d(x_pos, y_pos, z_pos));
    }

    if ((line == "ENDMDL") || (line == "END") || (_fl.eof())) {
      break;
    }
  }  // while std::getline
  if (!_topology && (bead_count > 0) && bead_count != top.BeadCount())
    throw std::runtime_error(
        "number of beads in topology and trajectory differ");

  ////////////////////////////////////////////////////////////////////////////////
  // Sort data and determine atom structure, connect with top.molecules,
  // bonds)
  ////////////////////////////////////////////////////////////////////////////////
  // Extra processing is done if the file is a topology file, in which case the
  // atoms must be sorted into molecules and the bond interactions recorded
  if (_topology) {
    // Now we need to add the bond pairs
    // WARNING We are assuming the atom ids are contiguous with no gaps

    // First int  - is the index of the atom
    // Second int - is the index of the molecule
    std::map<int, int> atm_molecule;

    // First int  - is the index of the molecule
    // list<int>  - is a list of the atoms in the molecule
    std::unordered_map<int, std::list<int>> molecule_atms;

    // Keep track of the number of molecules we have created through an index
    int mol_index = 0;

    // Cycle through all bonds
    for (auto row = bond_pairs.begin(); row != bond_pairs.end(); row++) {

      int atm_id1 = row->at(0);
      int atm_id2 = row->at(1);
      // Check to see if either atm referred to in the bond is already
      // attached to a molecule
      auto mol_iter1 = atm_molecule.find(atm_id1);
      auto mol_iter2 = atm_molecule.find(atm_id2);

      // This means neither atom is attached to a molecule
      if (mol_iter1 == atm_molecule.end() && mol_iter2 == atm_molecule.end()) {
        // We are going to create a new row for a new molecule
        std::list<int> atms_in_mol;
        atms_in_mol.push_back(atm_id1);
        atms_in_mol.push_back(atm_id2);
        //        std::cout << "Molecule " << mol_index << "Adding atom " <<
        //        atms_in_mol << std::endl;
        molecule_atms[mol_index] = atms_in_mol;
        // Associate atm1 and atm2 with the molecule index
        atm_molecule[atm_id1] = mol_index;
        atm_molecule[atm_id2] = mol_index;
        // Increment the molecule index
        mol_index++;

        // This means only atm2 is attached to a molecule
      } else if (mol_iter1 == atm_molecule.end()) {
        // Add atm1 to the molecule that contains atm2
        molecule_atms[mol_iter2->second].push_back(atm_id1);
        // Associate atm1 with the molecule it is now part of
        atm_molecule[atm_id1] = mol_iter2->second;

        // This means only atm1 is attached to a molecule
      } else if (mol_iter2 == atm_molecule.end()) {
        // Add atm2 to the molecule that contains atm1
        molecule_atms[mol_iter1->second].push_back(atm_id2);
        // Associate atm1 with the molecule it is now part of
        atm_molecule[atm_id2] = mol_iter1->second;

      } else if (mol_iter1 != mol_iter2) {
        // This means both atm1 and atm2 are attached to a molecule
        // But if they are already attached to the same molecule there is
        // nothing else to be done.
        int chosen_mol;
        int obsolete_mol;
        // We will merge the atms to the molecule with the smallest index
        if (mol_iter1->second < mol_iter2->second) {
          chosen_mol = mol_iter1->second;
          obsolete_mol = mol_iter2->second;
        } else {
          chosen_mol = mol_iter2->second;
          obsolete_mol = mol_iter1->second;
        }

        // Now we will proceed to cycle through the atms that were in the now
        // obsolete molecule and make sure they are pointing to the new molecule
        for (auto atm_temp = molecule_atms[obsolete_mol].begin();
             atm_temp != molecule_atms[obsolete_mol].end(); atm_temp++) {

          atm_molecule[*atm_temp] = chosen_mol;
        }

        // Splicing will remove atoms from the now obsolete molecule and place
        // them on the chosen molecule.
        molecule_atms[chosen_mol].splice(molecule_atms[chosen_mol].end(),
                                         molecule_atms[obsolete_mol]);

        // Finally we will clear out the record of the obsolete molecule
        molecule_atms.erase(obsolete_mol);
      }
    }
#ifndef NDEBUG
    std::cerr << "Consistency check for pdbreader" << std::endl;
    int i = 0;
    for (auto lis = molecule_atms.begin(); lis != molecule_atms.end(); lis++) {
      std::cerr << "Molecule " << i << std::endl;
      std::cerr << "Atoms: ";
      for (auto atm_ind = lis->second.begin(); atm_ind != lis->second.end();
           atm_ind++) {
        std::cerr << *atm_ind << " ";
      }
      std::cerr << std::endl;
      i++;
    }
    std::cerr << std::endl;
#endif
    // Now that we know which interactions belong to which molecules we can:
    // 1 Add the molecules
    // 2 Add the bond interactions

    // Molecule map
    // First int - is the index of the molecule
    // Molecule* - is a pointer to the Molecule object
    std::map<int, typename Topology_T::container_t *> mol_map;
    std::unordered_map<int, std::vector<typename Topology_T::bead_t *>>
        molecule_beads;
    for (const std::pair<const int, std::list<int>> &mol_and_atom_ids :
         molecule_atms) {

      int molecule_id = mol_and_atom_ids.first;
      tools::StructureParameters params_mol;
      params_mol.set(tools::StructureParameter::MoleculeId, molecule_id);
      params_mol.set(tools::StructureParameter::MoleculeType,
                     tools::topology_constants::unassigned_molecule_type);
      mol_map[molecule_id] = &(top.CreateMolecule(params_mol));

      // Add all the atoms to the appropriate molecule object
      for (int atm_temp : mol_and_atom_ids.second) {
        molecule_beads[molecule_id].push_back(&top.getBead(atm_temp));
      }
    }

    for (std::pair<int, std::vector<typename Topology_T::bead_t *>> mol_b :
         molecule_beads) {
      for (typename Topology_T::bead_t *bead_f : mol_b.second) {
        mol_map[mol_b.first]->AddBead(*bead_f);
      }
    }
    // Cyle through the bonds and add them to the appropriate molecule
    for (std::vector<int> &bond_pair : bond_pairs) {

      int atm_id1 = bond_pair.at(0);
      int bond_index = 1;
      while (bond_index < static_cast<int>(bond_pair.size())) {
        int atm_id2 = bond_pair.at(bond_index);

        // Should be able to just look at one of the atoms the bond is attached
        // too because the other will also be attached to the same molecule.
        int mol_ind = atm_molecule[atm_id1];

        typename Topology_T::container_t *mi = mol_map[mol_ind];

        // Grab the id of the bead associated with the atom
        // It may be the case that the atom id's and bead id's are different
        int bead_id1 = atm_id1;
        int bead_id2 = atm_id2;
        mi->ConnectBeads(bead_id1, bead_id2);

        Interaction *ic = top.CreateInteraction(
            InteractionType::bond, "BONDS", bond_index, mi->getId(),
            std::vector<int>{bead_id1, bead_id2});
        mi->AddInteraction(ic);
        ++bond_index;
      }
    }

    // Finally we want to build an exclusion matrix
    top.RebuildExclusions();
  }
  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_PDBREADER_PRIV_H
