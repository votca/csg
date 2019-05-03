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
#ifndef VOTCA_CSG_BEADMAP_H
#define VOTCA_CSG_BEADMAP_H

#include "cgbeadstencil.h"
#include "csgtopology.h"
#include "molecule.h"
#include <string>
#include <vector>
#include <votca/tools/property.h>

namespace votca {
namespace csg {

class BoundaryCondition;

/**
 * @brief BeadMap for coarse graining an atomistic representation to a coarser
 * representation
 *
 * This abstract interface class that sets the requirements for creating other
 * classes that define how a set of atomistic beads could be used to update a
 * coarse grained representation. For instance if one was to use H2O as an
 * example we could create a single bead to determine how a water molecule
 * behaves as opposed to using 3 separate atomic particles. Though this class
 * does not create the coarse grained representation of water it does store the
 * relationship between a coarse grained and atomistic representation.
 *
 * The relationship defined by the map allows one to update the positions,
 * velocities, forces, masses and orientations of a coarse grained bead with
 * its atomistic pieces.
 *
 * Atomistic Representation          Coarse Grained Representation
 *
 *  H - O                                   H2O
 *      |
 *      H
 *
 * Note that the bead map does not have ownership over either coarse grained
 * representation or atomisitic representation it can only manipulate them
 * when they are passed in, it does not and should not have responsibility for
 * their memory managment.
 */
class BeadMap {
 public:
  virtual ~BeadMap(){};

  /**
   * @brief Applies mapping from atomistic positions to coarse grained bead
   *
   * Note: the cg bead cannot simply be returned via the return value. The Bead
   * map does not own the CGBead and it is only updateing some of the
   * information in the cg bead, thus it makes sense to pass it in as a
   * pointer.
   *
   * @param boundaries - required to be a pointer to take advantage of
   * polymorphic behavior
   * @param atomistic_beads - pointers to the atomistic beads
   * @param cg_bead - pointer to the cg bead
   */
  virtual void UpdateCGBead(const BoundaryCondition *boundaries,
                            std::map<std::string, const Bead *> atomistic_beads,
                            Bead *cg_bead) const = 0;

  /**
   * @brief Method initializes the paramters important for describing how the
   * atomistic beads will be mapped to the cg bead.
   *
   * More specifically determines how the forces acting on the individual
   * atomistic beads will affect the coarse grained bead, by calculating force
   * weights.
   *
   * @param subbeads
   * @param weights - passed by copy, because they are temporary renormalized
   * within the method
   * @param ds
   */
  virtual void InitializeBeadMap(const std::vector<std::string> &subbeads,
                                 std::vector<double> weights,
                                 std::vector<double> ds) = 0;

  void InitializeBeadMap(const std::vector<std::string> &subbeads,
                         const std::vector<double> &weights) {
    std::vector<double> empty;
    InitializeBeadMap(subbeads, weights, empty);
  }

  std::vector<std::string> getAtomicBeadNames() const;

  virtual std::unique_ptr<BeadMap> Clone() const = 0;

 protected:
  /**
   * @brief BeadMap constructor is hidden to help with encapsulation, only the
   * AtomToCGMoleculeMapper class can now manipulate the BeadMap.
   */
  BeadMap(){};
  struct element_t {
    double weight_;
    double force_weight_;
  };
  std::unordered_map<std::string, element_t> matrix_;
  friend class AtomToCGMoleculeMapper;
};

/*******************************************************
    Linear map for spherical beads
*******************************************************/

class Map_Sphere : public BeadMap {
 public:
  void UpdateCGBead(const BoundaryCondition *boundaries,
                    std::map<std::string, const Bead *> atomistic_beads,
                    Bead *cg_bead) const override;

  void InitializeBeadMap(const std::vector<std::string> &subbeads,
                         std::vector<double> weights,
                         std::vector<double> ds) override;

  virtual std::unique_ptr<BeadMap> Clone() const override {
    return std::unique_ptr<BeadMap>(new Map_Sphere(*this));
  }

 protected:
  Map_Sphere() {}

  /**
   * @brief Adds an atomistic bead with its weight to the coarse grained bead
   *
   * @param atomic_bead_name
   * @param weight
   * @param force_weight
   */
  void AddAtomisticBead(const std::string &atomic_bead_name,
                        const double &weight, const double &force_weight);

  friend class AtomToCGMoleculeMapper;
};

/*******************************************************
    Linear map for ellipsoidal bead
*******************************************************/
class Map_Ellipsoid : public Map_Sphere {
 public:
  void UpdateCGBead(const BoundaryCondition *boundaries,
                    std::map<std::string, const Bead *> atomistic_beads,
                    Bead *cg_bead) const override;

  virtual std::unique_ptr<BeadMap> Clone() const override {
    return std::unique_ptr<BeadMap>(new Map_Ellipsoid(*this));
  }

 protected:
  Map_Ellipsoid() {}
  friend class AtomToCGMoleculeMapper;
};

/*******************************************************
    Mapper class, collection of maps
*******************************************************/

// Typedef is defined to reduce the verbosity of the code
// First element of the pair
// cg_molecule_id
// Second element of the pair is the map with key
// cg_bead_id
// Finally the map stores
// vector<int> { atom_id1,atom_id2,etc... }
typedef std::pair<int, std::map<int, std::vector<std::pair<std::string, int>>>>
    CGMolToAtom;

/**
 * @brief Bead Mapper for a whole molecule
 *
 * Unlike the BeadMap class which is responsible for the relationship of a few
 * atomistic bead to a single coarse grained bead, the AtomToMoleculeMapper is
 * a container for the beadmaps used to describe a full molecule. Using propane
 * as an illustration, where -- indicate where one coarse grained bead starts
 * and another begins:
 *
 *          H6       H7      H8
 *           |       |       |
 *      H1 - C2  --  C3  --  C4 - H5
 *           |       |       |
 *          H9       H10     H11
 *
 * Bead Map 1   Bead Map 2   Bead Map 3
 *
 * Here all three bead maps would be stored in this class.
 */
class AtomToCGMoleculeMapper {
 public:
  AtomToCGMoleculeMapper(std::string atom_molecule_type,
                         std::string cg_molecule_type)
      : atom_molecule_type_(atom_molecule_type),
        cg_molecule_type_(cg_molecule_type){};

  void InitializeMoleculeMap(
      const std::unordered_map<std::string, CGBeadStencil> &bead_maps_info,
      const std::vector<std::string> &bead_order);

  // Pass in a map containing the names of all the atomistic beads in the
  // molecule and pointers to them
  void UpdateCGMolecule(CSG_Topology &atom_top, CSG_Topology &cg_top,
                        CGMolToAtom cgmolid_cgbeadid_atomicbeadnames_ids);

  /***************************
   * Gang of 5
   **************************/
  /**
   * @brief Note that the gang of 5 must be explicity defined because the bead
   * maps are stored as unique_ptrs, unique pointers are used so that the
   * BeadMaps can be treated polymorphically.
   */

  /// Copy Constructor
  AtomToCGMoleculeMapper(const AtomToCGMoleculeMapper &other);
  /// Move Constructor
  AtomToCGMoleculeMapper(AtomToCGMoleculeMapper &&other) noexcept;
  /// Move Assignment
  AtomToCGMoleculeMapper &operator=(AtomToCGMoleculeMapper &&other) noexcept;
  /// Copy Assignment
  AtomToCGMoleculeMapper &operator=(const AtomToCGMoleculeMapper &other);
  ~AtomToCGMoleculeMapper();

 protected:
  // Molecule atomistic_molecule_, cg_molecule_;
  std::string atom_molecule_type_;
  std::string cg_molecule_type_;

  // contains each bead types and a vector of all the bead names of that type
  std::unordered_map<std::string, std::vector<std::string>>
      bead_type_and_names_;
  // Needs to be a unique_ptr to take advantage of polymorphism
  std::unordered_map<std::string, std::unique_ptr<BeadMap>>
      cg_bead_name_and_maps_;

  /**
   * @brief Returns pointers of the atomic beads that are within the cg bead
   *
   * @param[in] atom_top
   * @param[in] cgmolid_cgbeadid_atomicbeadids
   * @param[in] cg_bead_id
   *
   * @return map with the atomic bead names and the pointers to them
   */
  std::map<std::string, const Bead *> getAtomicNamesAndBeads_(
      const CSG_Topology &atom_top,
      const std::vector<std::pair<std::string, int>> &atomic_names_and_ids)
      const;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMAP_H
