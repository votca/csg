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
#ifndef VOTCA_CSG_ATOMTOCGCONVERTER_H
#define VOTCA_CSG_ATOMTOCGCONVERTER_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "beadmap.h"
#include "cgmoleculestencil.h"
#include "exclusionlist.h"
#include "molecule.h"
#include <boost/bimap.hpp>
#include <votca/tools/property.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {

/**
 * @brief Converter class provides methods for converting a atomistic molecule
 * to a coarse grained molecule
 *
 * Note the converter should not have ownership of any of the topology objects
 * because it should be able to convert between any topology object as long as
 * it has the stencil, this means that for each operation that it performs
 * an atomistic (input) object and a coarse grained (output) object must be
 * provided.
 */
class AtomCGConverter {
 public:
  /// Constructor
  AtomCGConverter(std::vector<std::string> ignore_atom_mol_types);

  /// Define parameterless constructor
  AtomCGConverter(){};

  /**
   * @brief Loads the .xml file describing the cg molecule
   *
   * Will load the information from the cg.xml file and create a stencil from
   * the input. For each cg molecule defined.
   *
   * @param[in] filename should be an .xml file
   */
  void LoadMoleculeStencil(const std::string &filename);

  /**
   * @brief Provided the atomic type of the molecule returns the cg type of the
   * molecule
   *
   * @param[in] atomistic_mol_type
   *
   * @return cg molecule type as a string
   */
  const std::string &getCGMoleculeType(const std::string &atom_mol_type) const;

  /**
   * @brief Provided the cg type of the molecule returns the atomic type of the
   * molecule
   *
   * @param cg_mol_type
   *
   * @return atomistic molecule type as a string
   */
  const std::string &getAtomisticMoleculeType(
      const std::string &cg_mol_type) const;

  /**
   * @brief Converts a atomic topology to a coarsegrained topology
   *
   * @param atomic_top_in
   * @param cg_top_out
   */
  Topology Convert(const Topology &atom_top_in);

  /**
   * @brief Updates the cg topology object using the atomistic topology
   *
   * The properties updated include the time step, bead positions, forces and
   * velocities etc...
   *
   * @param atomic_top_in
   * @param cg_top_out
   */
  void Update(const Topology &atom_top_in, Topology &cg_top_out) const;

  std::vector<std::string> getAtomicBeadNamesOfCGBead(
      const std::string &cg_mol_type, const std::string &cg_bead_type) const;

  bool AtomisticMoleculeTypeExist(const std::string &atom_mol_type) const;

 private:
  std::map<int, std::map<int, std::vector<std::pair<std::string, int>>>>
      cgmolid_cgbeadid_atombeadnames_and_ids_;

  std::unordered_set<std::string> atom_mol_types_to_ignore_;
  /**
   * @brief Atomic and cg bimap stores the relationship between the two types
   * of molecules
   */
  boost::bimap<std::string, std::string> atom_and_cg_mol_types_;

  /**
   * @brief Stores the stencil of the cg molecule
   *
   * string is the cg_mol_type
   */
  std::unordered_map<std::string, CGMoleculeStencil> cg_mol_and_stencil_;

  // First string atomic_molecule_type
  // Second string cg_mol_type
  // The molecule map
  std::unordered_map<std::string,
                     std::unordered_map<std::string, AtomToCGMoleculeMapper>>
      mol_names_and_maps_;

  /**
   * @brief Maps a atomistic molecule to a cg molecule and adds it to the cg
   * topology
   *
   * Note that the cg and atomistic descriptions should have the same molecule
   * ids, this means that the molecule ids will remain consistent between the
   * coarse grained and atomistic descriptions.
   *
   * @param[in] atomistic_mol
   * @param[in,out] cg_top_out
   */
  //  void ConvertAtomisticMoleculeToCGAndAddToCGTopology_(
  Molecule &CoarseGrainMolecule_(const Topology &atom_top, Topology &cg_top_out,
                                 const Molecule &atom_mol);

  /**
   * @brief
   *
   * Must pass all the atomic molecules bead ids in, the bead ids must increase
   * in magnitude in correlation to how the .xml files are defined
   *
   * @param cg_or_atomic_mol_type
   * @param atomic_bead_ids
   *
   * @return
   */
  std::unordered_map<int, std::string> MapAtomicBeadIdsToAtomicBeadNames_(
      const std::string &cg_or_atom_mol_type,
      const std::vector<int> &atom_bead_ids) const;

  bool CheckThatBeadCountAndInteractionTypeAreConsistent_(
      const std::string &interaction_type, const size_t &bead_count) const
      noexcept;

  std::vector<CGBeadStencil> ParseBeads_(tools::Property &options);

  std::vector<CGInteractionStencil> ParseBonded_(tools::Property &options);

  void ParseBeadMaps_(
      tools::Property &options_in,
      std::unordered_map<std::string, CGBeadStencil> &bead_maps_info);

  std::map<int, std::vector<std::pair<std::string, int>>>
      MapAtomicBeadsToCGBeads_(const Topology &atom_top, Topology &cg_top_out,
                               Molecule &cg_mol,
                               const CGMoleculeStencil &stencil) const;

  void CreateInteractions_(Topology &cg_top_out, Molecule &cg_mol,
                           const CGMoleculeStencil &stencil) const;

  std::map<int, std::vector<std::pair<std::string, int>>>
      MapAtomicMoleculeToCGMolecule_(const Topology &atom_top,
                                     Topology &cg_top_out, const int &mol_id);
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_ATOMTOCGCONVERTER_H
