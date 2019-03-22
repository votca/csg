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

#ifndef VOTCA_CSG_BEADMAP_H
#define VOTCA_CSG_BEADMAP_H

#include "cgbeadstencil.h"
#include "csgtopology.h"
#include "molecule.h"
#include <string>
#include <vector>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

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
  virtual void Apply(const BoundaryCondition *boundaries,
                     std::map<std::string, const Bead *> atomistic_beads,
                     Bead *cg_bead) = 0;
  virtual void Initialize(const std::vector<std::string> subbeads,
                          std::vector<double> weights,
                          std::vector<double> ds) = 0;

  void Initialize(std::vector<std::string> subbeads,
                  std::vector<double> weights) {
    std::vector<double> empty;
    Initialize(subbeads, weights, empty);
  }

  std::vector<std::string> getAtomicBeadNames() const;

  virtual std::unique_ptr<BeadMap> Clone() const = 0;

 protected:
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
  void Apply(const BoundaryCondition *boundaries,
             std::map<std::string, const Bead *> atomistic_beads,
             Bead *cg_bead) override;

  void Initialize(const std::vector<std::string> subbeads,
                  std::vector<double> weights, std::vector<double> ds) override;

  virtual std::unique_ptr<BeadMap> Clone() const override {
    return std::unique_ptr<BeadMap>(new Map_Sphere(*this));
  }

 protected:
  Map_Sphere() {}
  void AddElem(std::string atomic_bead_name, double weight,
               double force_weight);

  friend class AtomToCGMoleculeMapper;
};

inline void Map_Sphere::AddElem(std::string atomic_bead_name, double weight,
                                double force_weight) {
  element_t el;
  el.weight_ = weight;
  el.force_weight_ = force_weight;
  matrix_[atomic_bead_name] = el;
}

/*******************************************************
    Linear map for ellipsoidal bead
*******************************************************/
class Map_Ellipsoid : public Map_Sphere {
 public:
  void Apply(const BoundaryCondition *boundaries,
             std::map<std::string, const Bead *> atomistic_beads,
             Bead *cg_bead) override;

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
  AtomToCGMoleculeMapper(){};
  AtomToCGMoleculeMapper(std::string atom_molecule_type,
                         std::string cg_molecule_type)
      : atom_molecule_type_(atom_molecule_type),
        cg_molecule_type_(cg_molecule_type){};
  ~AtomToCGMoleculeMapper();

  void Initialize(std::unordered_map<std::string, CGBeadStencil> bead_maps_info,
                  std::vector<std::string> bead_order);

  // Pass in a map containing the names of all the atomistic beads in the
  // molecule and pointers to them
  void Apply(
      CSG_Topology &atom_top, CSG_Topology &cg_top,
      std::pair<int, std::map<int, std::vector<std::pair<std::string, int>>>>
          cgmolid_cgbeadid_atomicbeadnames_ids);

  /***************************
   * Gang of 3
   **************************/
  /**
   * @brief Note that the gang of 3 must be explicity defined because the bead
   * maps are stored as unique_ptrs, unique pointers are used so that the
   * BeadMaps can be treated polymorphically.
   */
  AtomToCGMoleculeMapper(const AtomToCGMoleculeMapper &other);
  AtomToCGMoleculeMapper &operator=(AtomToCGMoleculeMapper &&other);
  AtomToCGMoleculeMapper &operator=(const AtomToCGMoleculeMapper other);

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
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMAP_H
