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
#ifndef VOTCA_CSG_TOPOLOGY_H
#define VOTCA_CSG_TOPOLOGY_H

#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "basemolecule.h"
#include "boundarycondition.h"
#include "exclusionlist.h"
#include "interaction.h"
#include "openbox.h"
#include "orthorhombicbox.h"
#include "topologytypecontainer.h"
#include "triclinicbox.h"

#include <boost/container/stable_vector.hpp>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

class BaseBead;
/**
 * @brief Creates a template topology class
 *
 * This class essentially is in control of all topology objects that can be
 * created. This means it controls the memory managment of any of all the
 * topology objects.
 *
 * @tparam Bead_T a bead/atom class
 * @tparam Molecule_T molecule class
 */
template <class Bead_T, class Molecule_T>
class TemplateTopology {
 public:
  /// constructor
  TemplateTopology() : time_(0.0), has_vel_(false), has_force_(false) {
    bc_ = OpenBox().Clone();
  }
  ~TemplateTopology();

  typedef Bead_T bead_t;
  typedef Molecule_T container_t;
  typename boost::container::stable_vector<Molecule_T>::iterator begin() {
    return molecules_.begin();
  }
  typename boost::container::stable_vector<Molecule_T>::iterator end() {
    return molecules_.end();
  }
  typename boost::container::stable_vector<Molecule_T>::const_iterator begin()
      const {
    return molecules_.begin();
  }
  typename boost::container::stable_vector<Molecule_T>::const_iterator end()
      const {
    return molecules_.end();
  }
  /**
   * \brief Cleans up all the stored data
   */
  void Cleanup();

  /**
   *  \brief checks weather molecules with the same name really contain the same
   * number of beads
   */
  void CheckMoleculeNaming(void) const;

  /**
   * \brief number of molecules in the system
   * @return number of molecule in topology
   */
  size_t MoleculeCount() const { return molecules_.size(); }

  /**
   * number of beads in the system
   * @return number of beads in the system
   */
  size_t BeadCount() const { return beads_.size(); }

  size_t InteractionCount() const { return interactions_.size(); }
  /**
   * access containter with all beads
   * @return bead container
   */
  // BeadContainer &Beads() { return beads_; }

  /**
   * access  containter with all molecules
   * @return molecule container
   */
  const std::vector<Molecule_T> &Molecules() { return molecules_; }

  /**
   * access containter with all bonded interactions
   * @return bonded interaction container
   */
  const std::vector<std::unique_ptr<Interaction>> &BondedInteractions() const {
    return interactions_;
  }

  /**
   * \brief Returns a pointer to the bead with index i
   *
   * @param[in] int index is the index of the bead in the internal bead
   * container
   *
   * @return Bead * is a pointer to the bead
   **/
  Bead_T &getBead(const int id) {
    assert(beads_.count(id) &&
           "Cannot access bead with provided id because it is not stored in "
           "the topology obeject");
    return beads_.at(id);
  }

  std::unordered_map<int, const Eigen::Vector3d *> getBeadPositions(
      std::vector<int> bead_ids) const {
    std::unordered_map<int, const Eigen::Vector3d *> bead_id_and_positions;
    for (int &bead_id : bead_ids) {
      bead_id_and_positions[bead_id] = &(beads_.at(bead_id).getPos());
    }
    return bead_id_and_positions;
  }

  /**
   * @brief Determine if a molecule with the specified id exists or not
   *
   * @param[in] molecule_id
   *
   * @return true if it does exist, else false
   */
  bool MoleculeExist(const int molecule_id) const noexcept {
    return molecules_map_.count(molecule_id);
  }

  bool BoundaryExist() const noexcept { return bc_ != nullptr; }
  bool BeadExist(const int bead_id) const { return beads_.count(bead_id); }

  bool BeadTypeExist(const std::string bead_type) const {
    return type_container_.BeadTypeExist(bead_type);
  }

  bool MoleculeTypeExist(const std::string molecule_type) const {
    return type_container_.MoleculeTypeExist(molecule_type);
  }
  /**
   * \brief Returns a constant pointer to the bead with index i
   *
   * @param[in] int index is the index of the bead in the internal bead
   * container
   *
   * @return const Bead * is a pointer to the bead
   **/
  const Bead_T &getBeadConst(const int id) const {
    assert(beads_.count(id) &&
           "Cannot access bead with provided id because it is not stored in "
           "the topology obeject");
    return beads_.at(id);
  }

  /**
   * @brief Grabs a molecule with the specified id
   *
   * @param[in] id
   *
   * @return raw pointer to the molecule
   */
  Molecule_T &getMolecule(const int &id) {
    assert(molecules_map_.count(id) &&
           "Cannot access molecule with provided id because it is not stored "
           "in the topology object.");
    return *(molecules_map_[id]);
  }

  /**
   * @brief Graps a const pointer to a molecule
   *
   * @param[in] id of the molecule
   *
   * @return const pointer to the molecule
   */
  const Molecule_T &getMolecule(const int &id) const {
    assert(molecules_map_.count(id) &&
           "Cannot access const molecule with provided id because it is not "
           "stored in the topology object.");
    return *(molecules_map_.at(id));
  }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() {
    molecules_.clear();
    molecules_map_.clear();
  }

  /**
   * \brief copy topology data of different topology
   *
   * Copies everything, such that everything in the topology object is a clope
   * of **top**, it will have ownership of all its attributes.
   * \param top topology to copy from
   */
  void Copy(const TemplateTopology<Bead_T, Molecule_T> &top);

  void CopyBoundaryConditions(const TemplateTopology<Bead_T, Molecule_T> &top);
  /**
   *  \brief rename all the molecules in range
   * \param range range string of type 1:2:10 = 1, 3, 5, 7, ...
   * \param name new name of molecule
   * range is a string which is parsed by RangeParser,
   */
  void RenameMoleculesType(std::string range, const std::string name);

  /**
   *  \brief rename all the bead types
   * \param name current rame of the bead type
   * \param newname new name of bead type
   */
  void RenameBeadsType(const std::string old_type, const std::string new_type);

  /**
   *  \brief set the mass of all the beads of a certain type
   * \param name the bead type
   * \param value mass value
   */
  void setBeadOfGivenTypeToNewMass(const std::string type, const double mass);

  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const Eigen::Matrix3d &box, BoundaryCondition::eBoxtype boxtype =
                                              BoundaryCondition::typeAuto);

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const Eigen::Matrix3d &getBox() const { return bc_->getBox(); }

  /**
   * set the time of current frame
   * \param t simulation time in ns
   */
  void setTime(const double &time) { time_ = time; }

  /**
   * get the time of current frame
   * \return simulation time in ns
   */
  double getTime() const {
    assert(time_ >= 0 && "Cannot get time has not been set.");
    return time_;
  }

  /**
   * set the step number of current frame
   * \param s step number
   */
  void setStep(const int &step) { step_ = step; }

  /**
   * get the step number of current frame
   * \return step number
   */
  int getStep() const {
    assert(step_ > -1 && "Cannot get step it has not been set.");
    return step_;
  };

  /**
   * Sets the particle group. (For the H5MD file format)
   * \param particle_group The name of a particle group.
   */
  void setParticleGroup(std::string particle_group) {
    particle_group_ = particle_group;
  }

  /**
   * Gets the particle group.
   * \return The name of a particle group.
   */
  std::string getParticleGroup() const { return particle_group_; }

  /**
   * \brief pbc correct distance of two beads
   * \param bead1 index of first bead
   * \param bead2 index of second bead
   * \return distance vector
   *
   * calculates the smallest distance between two beads with correct treatment
   * of pbc
   */
  Eigen::Vector3d getDist(const int bead1, const int bead2) const;

  /**
   * \brief calculate shortest vector connecting two points
   * \param r1 first point
   * \param r2 second point
   * \return distance vector
   *
   * calculates the smallest distance between two points with correct treatment
   * of pbc
   */
  Eigen::Vector3d BCShortestConnection(const Eigen::Vector3d &r_i,
                                       const Eigen::Vector3d &r_j) const;

  /**
   * \brief return the shortest box size
   * \return shortest size
   *
   * Calculates the shortest length to connect two sides of the box
   */
  double getShortestBoxDimension() const;

  /**
   *  calculates the box volume
   *  \return box volume
   */
  double BoxVolume() const { return bc_->BoxVolume(); }

  /**
   *  rebuild exclusion list
   */
  void RebuildExclusions();

  /**
   * @brief Grabs all the interactions in a specified group
   *
   * @param[in] group
   *
   * @return returns a list of Interactions
   */
  std::list<Interaction *> InteractionsInGroup(const std::string &group) const;

  /**
   * access exclusion list
   * \return exclusion list
   */
  ExclusionList &getExclusions() { return exclusions_; }

  const ExclusionList &getExclusionsConst() const { return exclusions_; }

  /**
   * @brief Grabs all the ids of all the beads stored in the topology
   *
   * @return vector of the bead ids
   */
  std::vector<int> getBeadIds() const;

  /**
   * @brief Grabs all the ids of the molecules stored in the topology
   *
   * @return vector of all the molecule ids
   */
  std::vector<int> getMoleculeIds() const;

  /**
   * @brief Returns a map of all the residues and their respective types in the
   * molecule specified by `molecule_id`
   *
   * @param[in] molecule_id
   *
   * @return residue ids "each one is unique" and the type of residue
   */
  std::map<int, std::string> getResidueIdsAndTypesInMolecule(
      int molecule_id) const;

  /**
   * @brief The boundary condition type of the topology
   *
   * @return box type
   */
  BoundaryCondition::eBoxtype getBoxType() const { return bc_->getBoxType(); }

  /**
   * @brief Returns a pointer to the boundary condition object
   *
   * @return const pointer
   */
  const BoundaryCondition *getBoundaryCondition() const { return bc_.get(); }

  template <typename iteratable>
  void InsertExclusion(Bead_T *bead1, iteratable &l);

  bool HasVel() { return has_vel_; }
  void SetHasVel(const bool v) { has_vel_ = v; }

  bool HasForce() { return has_force_; }
  void SetHasForce(const bool v) { has_force_ = v; }

  /**
   * @brief Returns the type id of the bead with with id `bead_id`
   *
   * Each bead has an id, and each bead has a type, the types can also be given
   * an id it is this id that is returned.
   *
   * @param[in] bead_id
   *
   * @return bead_type_id
   */
  int getBeadTypeId(int bead_id) const {
    return type_container_.getBeadTypeId(beads_.at(bead_id).getType());
  }

  std::vector<std::string> getBeadTypes() const {
    return type_container_.getBeadTypes();
  }

  std::vector<std::string> getResidueTypes() const {
    return type_container_.getResidueTypes();
  }
  /**
   * @brief Returns the available molecule types and their type ids
   *
   * @return
   */
  const std::vector<std::string> getMoleculeTypes() const {
    return type_container_.getMoleculeTypes();
  }

 protected:
  std::unique_ptr<BoundaryCondition> bc_ = nullptr;

  BoundaryCondition::eBoxtype autoDetectBoxType_(
      const Eigen::Matrix3d &box) const;

  /// beads in the topology, must be a stable map, (unordered map does not work
  // because the pointers get invalidated any time the unordered_map is resized)
  std::map<int, Bead_T> beads_;

  /// molecules in the topology
  // std::unordered_map<int, Molecule_T *> molecules_map_;
  boost::container::stable_vector<Molecule_T> molecules_;
  std::unordered_map<int, Molecule_T *> molecules_map_;

  /// bonded interactions in the topology
  std::vector<std::unique_ptr<Interaction>> interactions_;

  ExclusionList exclusions_;

  std::map<std::string, int> interaction_groups_;

  std::map<std::string, std::list<Interaction *>> interactions_by_group_;

  double time_ = -1.0;
  int step_ = -1;
  bool has_vel_ = false;
  bool has_force_ = false;

  /// The particle group (For H5MD file format)
  std::string particle_group_ = "unassigned";

  /**
   * @brief Keeps track of the various different topology objects stored in
   * the topology object
   */
  TopologyTypeContainer type_container_;
};

template <class Bead_T, class Molecule_T>
template <typename iteratable>
void TemplateTopology<Bead_T, Molecule_T>::InsertExclusion(Bead_T *bead1,
                                                           iteratable &l) {
  exclusions_.InsertExclusion(bead1, l);
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::Cleanup() {

  beads_.clear();
  molecules_.clear();
  molecules_map_.clear();
  type_container_.Clear();
  interactions_.clear();
  bc_ = OpenBox().Clone();
}

template <class Bead_T, class Molecule_T>
TemplateTopology<Bead_T, Molecule_T>::~TemplateTopology() {
  Cleanup();
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::CopyBoundaryConditions(
    const TemplateTopology<Bead_T, Molecule_T> &top) {
  assert(bc_ != nullptr &&
         "Cannot copy boundary condition it has not been assigned");
  bc_ = top.bc_->Clone();
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::Copy(
    const TemplateTopology<Bead_T, Molecule_T> &top) {
  step_ = top.step_;
  time_ = top.time_;
  has_vel_ = top.has_vel_;
  has_force_ = top.has_force_;
  beads_ = top.beads_;

  molecules_.resize(top.molecules_.size());
  size_t index = 0;
  for (auto &molecule : top.molecules_) {
    molecules_.at(index) = molecule;
    molecules_map_[molecule.getId()] = &(molecules_.at(index));
    ++index;
  }

  bc_ = top.bc_->Clone();
  type_container_ = top.type_container_;
  particle_group_ = top.particle_group_;

  for (const std::unique_ptr<Interaction> &interaction : top.interactions_) {
    interactions_.push_back(interaction->Clone());
  }

  interaction_groups_ = top.interaction_groups_;

  for (const std::pair<std::string, std::list<Interaction *>>
           &interaction_name_and_list : top.interactions_by_group_) {
    std::string interaction_name = interaction_name_and_list.first;
    for (const Interaction *interaction : interaction_name_and_list.second) {
      int index_of_interaction = interaction->getIndex();
      interactions_by_group_[interaction_name].push_back(
          interactions_[index_of_interaction].get());
    }
  }

  RebuildExclusions();
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::setBox(
    const Eigen::Matrix3d &box, BoundaryCondition::eBoxtype boxtype) {
  // determine box type automatically in case boxtype==typeAuto
  if (boxtype == BoundaryCondition::typeAuto) {
    boxtype = autoDetectBoxType_(box);
  }

  switch (boxtype) {
    case BoundaryCondition::typeTriclinic:
      bc_ = TriclinicBox().Clone();
      break;
    case BoundaryCondition::typeOrthorhombic:
      bc_ = OrthorhombicBox().Clone();
      break;
    default:
      bc_ = OpenBox().Clone();
      break;
  }

  bc_->setBox(box);
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::RenameMoleculesType(
    std::string range_molecule_ids, const std::string type) {
  tools::RangeParser rp;
  tools::RangeParser::iterator molecule_id_ptr;

  rp.Parse(range_molecule_ids);
  for (molecule_id_ptr = rp.begin(); molecule_id_ptr != rp.end();
       ++molecule_id_ptr) {
    if ((unsigned int)*molecule_id_ptr > molecules_.size()) {
      throw std::runtime_error(
          std::string("RenameMoleculesType: num molecules smaller than"));
    }
    getMolecule(*molecule_id_ptr).setType(type);
  }
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::RenameBeadsType(
    const std::string old_type, std::string new_type) {
  for (std::pair<const int, Bead_T> &id_and_bead : beads_) {
    std::string bead_type = id_and_bead.second.getType();
    if (tools::wildcmp(bead_type.c_str(), old_type.c_str())) {
      id_and_bead.second.setType(new_type);
    }
  }
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::setBeadOfGivenTypeToNewMass(
    std::string type, double mass) {

  for (std::pair<const int, Bead_T> &id_and_bead : beads_) {
    std::string bead_type = id_and_bead.second.getType();
    if (tools::wildcmp(bead_type.c_str(), type.c_str())) {
      id_and_bead.second.setMass(mass);
    }
  }
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::CheckMoleculeNaming(void) const {
  std::unordered_map<std::string, size_t>
      number_of_beads_in_each_molecular_type;

  // for (const std::pair<const int, Molecule_T> &id_and_molecule : molecules_)
  // {
  //  for (const std::pair<int, Molecule_T> &molecule : molecules_) {
  for (const Molecule_T &molecule : molecules_) {
    std::string type = molecule.getType();
    size_t bead_count = molecule.BeadCount();
    if (number_of_beads_in_each_molecular_type.count(type)) {
      if (number_of_beads_in_each_molecular_type.at(type) != bead_count) {
        throw std::runtime_error(
            "There are molecules which have the same name but different number "
            "of bead please check the section manual topology handling in the "
            "votca manual");
      }
    } else {
      number_of_beads_in_each_molecular_type[type] = bead_count;
    }
  }
}

template <class Bead_T, class Molecule_T>
std::list<Interaction *>
    TemplateTopology<Bead_T, Molecule_T>::InteractionsInGroup(
        const std::string &group) const {
  auto iter = interactions_by_group_.find(group);
  if (iter == interactions_by_group_.end()) return std::list<Interaction *>();
  return iter->second;
}

template <class Bead_T, class Molecule_T>
Eigen::Vector3d TemplateTopology<Bead_T, Molecule_T>::BCShortestConnection(
    const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const {
  return bc_->BCShortestConnection(r_i, r_j);
}

template <class Bead_T, class Molecule_T>
Eigen::Vector3d TemplateTopology<Bead_T, Molecule_T>::getDist(
    const int bead1, const int bead2) const {
  return BCShortestConnection(getBead(bead1).getPos(), getBead(bead2).getPos());
}

template <class Bead_T, class Molecule_T>
BoundaryCondition::eBoxtype
    TemplateTopology<Bead_T, Molecule_T>::autoDetectBoxType_(
        const Eigen::Matrix3d &box) const {
  // set the box type to OpenBox in case "box" is the zero matrix,
  // to OrthorhombicBox in case "box" is a diagonal matrix,
  // or to TriclinicBox otherwise
  if (box.isApproxToConstant(0)) {
    return BoundaryCondition::typeOpen;
  } else if ((box - Eigen::Matrix3d(box.diagonal().asDiagonal()))
                 .isApproxToConstant(0)) {
    return BoundaryCondition::typeOrthorhombic;
  } else {
    return BoundaryCondition::typeTriclinic;
  }
  return BoundaryCondition::typeOpen;
}

template <class Bead_T, class Molecule_T>
double TemplateTopology<Bead_T, Molecule_T>::getShortestBoxDimension() const {
  return bc_->getShortestBoxDimension();
}

template <class Bead_T, class Molecule_T>
std::vector<int> TemplateTopology<Bead_T, Molecule_T>::getBeadIds() const {
  std::vector<int> bead_ids;
  for (const std::pair<const int, Bead_T> id_and_bead : beads_) {
    bead_ids.push_back(id_and_bead.first);
  }
  return bead_ids;
}

template <class Bead_T, class Molecule_T>
std::vector<int> TemplateTopology<Bead_T, Molecule_T>::getMoleculeIds() const {
  std::vector<int> molecule_ids;
  // for (const std::pair<const int, Molecule_T> id_and_molecule : molecules_) {
  for (const std::pair<const int, Molecule_T *> id_and_molecule :
       molecules_map_) {
    molecule_ids.push_back(id_and_molecule.first);
  }
  return molecule_ids;
}

template <class Bead_T, class Molecule_T>
std::map<int, std::string>
    TemplateTopology<Bead_T, Molecule_T>::getResidueIdsAndTypesInMolecule(
        int molecule_id) const {

  assert(molecules_map_.count(molecule_id) &&
         "Molecule id does not exist in topology object");
  std::map<int, std::string> id_and_residue_type;
  std::vector<int> bead_ids = molecules_map_.at(molecule_id)->getBeadIds();
  for (const int &bead_id : bead_ids) {
    id_and_residue_type[beads_.at(bead_id).getResidueId()] =
        beads_.at(bead_id).getResidueType();
  }
  return id_and_residue_type;
}

template <class Bead_T, class Molecule_T>
void TemplateTopology<Bead_T, Molecule_T>::RebuildExclusions() {

  std::vector<Bead *> beads;
  for (std::unique_ptr<Interaction> &interaction : interactions_) {
    std::vector<int> bead_ids = interaction->getBeadIds();
    for (const int &bead_id : bead_ids) {
      beads.push_back(&beads_[bead_id]);
    }
  }
  exclusions_.ExcludeList(beads);
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TOPOLOGY_H
