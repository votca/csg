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

#ifndef VOTCA_CSG_CGMOLECULESTENCIL_H
#define VOTCA_CSG_CGMOLECULESTENCIL_H

#include "cgbeadstencil.h"
#include "cginteractionstencil.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <votca/tools/property.h>
#include <votca/tools/types.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

namespace votca {
namespace csg {

typedef boost::bimap<boost::bimaps::multiset_of<std::string>,
                     boost::bimaps::set_of<std::string>>
    multi_bimap;

/**
 * @brief Definition of coarse grained molecule
 *
 * This class is to define a stencil which can be used to create a coarse
 * grained molecule, by itself the stencil can furnish a Molecule class
 * instance with enough information to create a shell of the coarse grained
 * molecule. The stencil however, does not contain positions, velocities, force
 * and orientation data.
 *
 * The reason for this design choice is that a stencil can be reused to create
 * multiple coarse grained molecules of a given type. To create a stencil the
 * atomistic representation of a molecule and the coarse grained representation
 * of the same molecule type must be known. This is why the constructor
 * requires that both the atomistic and coarse grained types be provided on
 * instantiation.
 *
 * This class also provides functionality for understanding how the beads in
 * a coarse grained representation is related to it atomic version.
 *
 */
class CGMoleculeStencil {
 public:
  CGMoleculeStencil(std::string cg_molecule_type,
                    std::string atomistic_molecule_type)
      : cg_molecule_type_(cg_molecule_type),
        atomistic_molecule_type_(atomistic_molecule_type){};

  ~CGMoleculeStencil(){};

  // This assumes that the vector is ordered according to the creation of the
  // beads in the molecule Such the the first cg_bead points to the atoms in the
  // atomistic molecule with the smallest ids And the second cg_bead in the
  // vector points to the beads in the atomistic molecule with the next largest
  // ids etc...
  void AddBeadStencil(const std::vector<CGBeadStencil> &bead_info);

  /**
   * @brief Maps a vector of atomic bead id to their names
   *
   * Note this method makes several imporant assumptions.
   * 1. It assumes when the bead ids passed in are sorted in order that they
   * will line up with the atomic bead names that were passed in when calling
   * the AddBeadStencil method.
   *
   * E.g. If I had a propane molecule
   *
   * cg_bead_name "A1"
   * atomic_beads "C1" "H4" "H5" "H6"
   *
   * cg_bead_name "B1"
   * atomic_beads "C2" "H7" "H8"
   *
   * cg_bead_name "A2"
   * atomic_beads "C3" "H9" "H10" "H11"
   *
   * And I passed them to the CGMoleculeStencil in a vector in the depicted
   * order, and then were to pass in some atomic ids of a propane molecule, it
   * would not matter what the ids were as long as they belonged to the same
   * propane molecule. As an example we will consider two separate examples
   *
   * propane molecule 1 has the following bead ids in a vector
   *
   * 11 1 2 10 3 4 5 6 8 7 9
   *
   * Notice that they are not in order. The method will return them as follows:
   *
   * 1 C1
   * 2 H4
   * 3 H5
   * 4 H6
   * 5 C2
   * 6 H7
   * 7 H8
   * 8 C3
   * 9 H9
   * 10 H10
   * 11 H11
   *
   * Notice this matches the order in which the CGBeadStencils were passed in
   *
   * Now if I were to pass in a second propane molecule with the following
   * bead ids:
   *
   * 23 24 33 25 32 30 31 26 27 28 29
   *
   * The method would return them as follows:
   *
   * 23 C1
   * 24 H4
   * 25 H5
   * 26 H6
   * 27 C2
   * 28 H7
   * 29 H8
   * 30 C3
   * 31 H9
   * 32 H10
   * 33 H11
   *
   * @param[in] bead_ids - the number of bead ids passed in must be equal to the
   * number of atomic beads in the molecule.
   *
   * @return an unordered map containing the bead ids with their matched atomic
   * name
   */
  std::unordered_map<int, std::string> MapAtomicBeadIdsToAtomicBeadNames(
      std::vector<int> bead_ids) const;

  /**
   * @brief This function attempts to map the bead ids of a specific coarse
   * grained instantation to its bead names.
   *
   * To be clear bead names are not the same a bead types below. For instance
   * if I have 10 propane molecules in my system and I want to know how the
   * coarse grained bead names of propane number 3 than I could provide the
   * bead ids associated with propane 3, lets assume the ids are as follows,
   * where -- indicate where the coarse grained beads have been split up.
   *
   *                               H2         H5      H8
   *                                |         |        |
   *                           H1 - C0   --   C4   --  C7 - H10
   *                                |         |        |
   *                               H3         H6      H9
   *
   * Coarse Grained Bead Id:       9         10        11
   *
   * Coarse Grained Type:          A          B        A
   *
   * Coarse Grained Bead Name:     A1         B2       A3
   *
   * The input would be the bead ids { 9, 10 , 13 }, and the output would be a
   * map of pairs of the form:
   *
   * {{ 9, A1 }
   *  {10, B2 }
   *  {11, A2 }}
   *
   * @param[in] bead_ids to coarse grained beads, to work correctly all the
   * beads ids for a particular molecule instance must be provided.
   *
   * @return map to the names of the coarse grained beads
   */
  std::unordered_map<int, std::string> MapCGBeadIdsToCGBeadNames(
      std::vector<int> bead_ids) const;

  std::vector<std::string> getAtomicBeadNames(std::string cg_bead_name) const;

  std::vector<std::string> getCGBeadNames() const;

  std::string getCGBeadName(std::string atom_bead_name) const;

  void AddInteractionStencil(
      const std::vector<CGInteractionStencil> &interaction_info);

  const std::vector<CGBeadStencil> &getBeadStencil() const;

  const std::vector<CGInteractionStencil> &getInteractionStencil() const;

  const std::string &getCGMoleculeType() const;

  const std::string &getAtomisticMoleculeType() const;

 private:
  std::string cg_molecule_type_;
  std::string atomistic_molecule_type_;
  multi_bimap cg_and_atom_names_;
  // Must appear in the order it is read in from the .xml file
  std::vector<CGBeadStencil> bead_info_;
  std::vector<CGInteractionStencil> interaction_info_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGMOLECULESTENCIL_H
