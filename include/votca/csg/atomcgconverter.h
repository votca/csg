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

#ifndef VOTCA_CSG_ATOMTOCGCONVERTER_H
#define VOTCA_CSG_ATOMTOCGCONVERTER_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "exclusionlist.h"
#include "map.h"
#include "molecule.h"
#include <boost/bimap.hpp>
#include <votca/tools/property.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BoundaryCondition;
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
  AtomCGConverter(){};

  /// Destructor
  ~AtomCGConverter(){};

  /**
   * @brief Loads the .xml file describing the cg molecule
   *
   * Will load the information from the cg.xml file and create a stencil from
   * the input. For each cg molecule defined.
   *
   * @param[in] filename should be an .xml file
   */
  void LoadConversionStencil(std::string filename);

  /**
   * @brief Provided the atomic name of the molecule returns the cg name of the
   * molecule
   *
   * @param[in] atomistic_molecule_type
   *
   * @return cg molecule type as a string
   */
  const std::string &getCGMoleculeType(string atomistic_molecule_type) const;

  /**
   * @brief Provided the cg name of the molecule returns the atomic name of the
   * molecule
   *
   * @param cg_molecule_type
   *
   * @return atomistic molecule type as a string
   */
  const std::string &getAtomisticMoleculeType(string cg_molecule_type) const;

  /**
   * @brief Maps a atomistic molecule to a cg molecule and adds it to the cg
   * topology
   *
   * Note that the cg and atomistic descriptions should have the same molecule
   * ids, this means that the molecule ids will remain consistent between the
   * coarse grained and atomistic descriptions.
   *
   * @param[in] atomistic_molecule
   * @param[in,out] cg_top_out
   */
  void ConvertAtomisticMoleculeToCGAndAddToCGTopology(
      Molecule &atomistic_molecule, CSG_Topology &cg_top_out);

  /**
   * @brief
   *
   * Must pass all the atomic molecules bead ids in, the bead ids must increase
   * in magnitude in correlation to how the .xml files are defined
   *
   * @param cg_or_atomic_molecule_type
   * @param atomic_bead_ids
   *
   * @return
   */
  std::unordered_map<int, std::string> MapAtomicBeadIdsToAtomicBeadNames(
      string cg_or_atomic_molecule_type, vector<int> atomic_bead_ids);

  std::vector<std::string> getAtomicBeadNamesOfCGBead(string cg_molecule_type,
                                                      string cg_bead_type);

 private:
  /**
   * @brief Atomic and cg bimap stores the relationship between the two types
   * of molecules
   */
  boost::bimap<std::string, std::string> atomic_and_cg_molecule_types_;

  /**
   * @brief Stores the stencil of the cg molecule
   */
  std::unordered_map<std::string, CGStencil> cg_molecule_and_stencil;

  void CheckThatBeadCountAndInteractionTypeAreConsistent_(
      string interaction_type, size_t bead_count) const;

  void ParseBeads_(Property &options);

  void ParseBonded_(Property &options);

  std::unordered_map<std::string, int> CreateBeads_(Molecule *cg_mol,
                                                    CGStencil stencil,
                                                    CSG_Topology &cg_top_out);

  void CreateInteractions_(
      Molecule *cg_mol, CGStencil stencil, CSG_Topology &cg_top_out,
      std::unordered_map<std::string, int> bead_name_to_id);

  Molecule *CreateMolecule_(std::string cg_molecule_type, int molecule_id,
                            CSG_Topology &cg_top_out);
  // void ParseMapping(Property &options);

  /*beaddef_t *getBeadByCGName(const std::string &cg_name);
  Property *getMapByName(const std::string &map_name);

  TopologyMap topology_map_;*/
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_ATOMTOCGCONVERTER_H
