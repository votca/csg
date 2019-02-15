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

#ifndef _VOTCA_CSG_TOPOLOGY_H
#define _VOTCA_CSG_TOPOLOGY_H

#include <cassert>
#include <list>
#include <map>
#include <vector>

#include "bead.h"
#include "boundarycondition.h"
#include "exclusionlist.h"
#include "molecule.h"
#include "openbox.h"
#include "orthorhombicbox.h"
#include "triclinicbox.h"

#include <votca/tools/matrix.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class Interaction;
class ExclusionList;

typedef std::vector<Molecule *> MoleculeContainer;
typedef std::vector<Bead *> BeadContainer;
typedef std::vector<Interaction *> InteractionContainer;

/**
 * \brief topology of the whole system
 *
 * The Topology class stores the topology of the system like the beads, bonds,
 * molecules and residues.
 *
 **/
class Topology {
 public:
  /// constructor
  Topology()
      : _time(0.0),
        _has_vel(false),
        _has_force(false),
        max_residue_id_(bead_constants::residue_number_unassigned) {
    _bc = new OpenBox();
  }
  virtual ~Topology();

  /**
   * \brief Cleans up all the stored data
   */
  virtual void Cleanup();

  /**
   * \brief Creates a new Bead
   *
   * \param[in] symmetry symmetry of the bead, 1: spherical 3: ellipsoidal
   * \param[in] name name of the bead
   * \param[in] type bead type
   * \param[in] resnr residue number
   * \param[in] m mass
   * \param[in] q charge
   * \return pointer to created bead
   *
   * The function creates a new bead and adds it to the list of beads.
   */
  template <class T>
  T *CreateBead(TOOLS::byte_t symmetry, std::string name, std::string type,
                int residue_number, std::string residue_name,
                std::string molecule_name, double m, double q);

  /**
   * \brief creates a new molecule
   * \param name name of the molecule
   * \return pointer to created molecule
   */
  virtual Molecule *CreateMolecule(std::string name);

  /**
   *  \brief checks weather molecules with the same name really contain the same
   * number of beads
   */
  void CheckMoleculeNaming(void);

  /**
   * \brief create molecules based on blocks of atoms
   * \param[in] name molecule name
   * \param[in] first first bead
   * \param[in] nbeads number of beads per molecule
   * \param[in] nmolecules number of molecules
   */
  void CreateMoleculesByRange(std::string name, int first, int nbeads,
                              int nmolecules);

  /**
   * \brief number of molecules in the system
   * @return number of molecule in topology
   */
  int MoleculeCount() { return _molecules.size(); }

  /**
   * number of beads in the system
   * @return number of beads in the system
   */
  int BeadCount() { return _beads.size(); }

  /**
   * get molecule by index
   * @param index molecule number
   * @return pointer to molecule
   */
  Molecule *MoleculeByIndex(int index);

  /**
   * access containter with all beads
   * @return bead container
   */
  BeadContainer &Beads() { return _beads; }

  /**
   * access  containter with all molecules
   * @return molecule container
   */
  MoleculeContainer &Molecules() { return _molecules; }

  /**
   * access containter with all bonded interactions
   * @return bonded interaction container
   */
  InteractionContainer &BondedInteractions() { return _interactions; }

  void AddBondedInteraction(Interaction *ic);
  std::list<Interaction *> InteractionsInGroup(const std::string &group);

  /**
   * \brief Determine if a bead type exists.
   *
   * @return bool true if it has been registered
   **/
  bool BeadTypeExist(std::string type) const;

  /**
   * \brief Register the bead type with the topology object.
   *
   * Records are kept of the different bead types in the topology object. This
   * method stores the bead type.
   **/
  void RegisterBeadType(std::string name);

  /**
   * \brief Given a bead type this method returns the id associated with the
   * type
   *
   * @param[in] string name of the type
   * @return int the id of the type
   **/
  int getBeadTypeId(std::string type) const;

  /**
   * \brief Returns a pointer to the bead with index i
   *
   * @param[in] int index is the index of the bead in the internal bead
   * container
   * @return Bead * is a pointer to the bead
   **/
  Bead *getBead(const int index) const { return _beads[index]; }
  Molecule *getMolecule(const int i) const { return _molecules[i]; }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() { _molecules.clear(); }

  /**
   * \brief copy topology data of different topology
   * \param top topology to copy from
   */
  void CopyTopologyData(Topology *top);

  /**
   *  \brief rename all the molecules in range
   * \param range range string of type 1:2:10 = 1, 3, 5, 7, ...
   * \param name new name of molecule
   * range is a string which is parsed by RangeParser,
   */
  void RenameMolecules(std::string range, std::string name);

  /**
   *  \brief rename all the bead types
   * \param name current rame of the bead type
   * \param newname new name of bead type
   */
  void RenameBeadType(std::string name, std::string newname);

  /**
   *  \brief set the mass of all the beads of a certain type
   * \param name the bead type
   * \param value mass value
   */
  void SetBeadTypeMass(std::string name, double value);

  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const matrix &box, BoundaryCondition::eBoxtype boxtype =
                                     BoundaryCondition::typeAuto) {
    // determine box type automatically in case boxtype==typeAuto
    if (boxtype == BoundaryCondition::typeAuto) {
      boxtype = autoDetectBoxType(box);
    }

    if (_bc) {
      delete (_bc);
    }

    switch (boxtype) {
      case BoundaryCondition::typeTriclinic:
        _bc = new TriclinicBox();
        break;
      case BoundaryCondition::typeOrthorhombic:
        _bc = new OrthorhombicBox();
        break;
      default:
        _bc = new OpenBox();
        break;
    }

    _bc->setBox(box);
  };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const matrix &getBox() { return _bc->getBox(); };

  /**
   * set the time of current frame
   * \param t simulation time in ns
   */
  void setTime(double t) { _time = t; };

  /**
   * get the time of current frame
   * \return simulation time in ns
   */
  double getTime() { return _time; };

  /**
   * set the step number of current frame
   * \param s step number
   */
  void setStep(int s) { _step = s; };

  /**
   * get the step number of current frame
   * \return step number
   */
  int getStep() { return _step; };

  /**
   * Sets the particle group. (For the H5MD file format)
   * \param particle_group The name of a particle group.
   */
  void setParticleGroup(std::string particle_group) {
    _particle_group = particle_group;
  };

  /**
   * Gets the particle group.
   * \return The name of a particle group.
   */
  std::string getParticleGroup() { return _particle_group; };

  /**
   * \brief pbc correct distance of two beads
   * \param bead1 index of first bead
   * \param bead2 index of second bead
   * \return distance vector
   *
   * calculates the smallest distance between two beads with correct treatment
   * of pbc
   */
  vec getDist(int bead1, int bead2) const;

  /**
   * \brief calculate shortest vector connecting two points
   * \param r1 first point
   * \param r2 second point
   * \return distance vector
   *
   * calculates the smallest distance between two points with correct treatment
   * of pbc
   */
  vec BCShortestConnection(const vec &r1, const vec &r2) const;

  /**
   * \brief return the shortest box size
   * \return shortest size
   *
   * Calculates the shortest length to connect two sides of the box
   */
  double ShortestBoxSize();

  /**
   *  calculates the box volume
   *  \return box volume
   */
  double BoxVolume();

  /**
   *  rebuild exclusion list
   */
  void RebuildExclusions();

  /**
   * access exclusion list
   * \return exclusion list
   */
  ExclusionList &getExclusions() { return _exclusions; }

  BoundaryCondition::eBoxtype getBoxType() { return _bc->getBoxType(); }

  template <typename iteratable>
  void InsertExclusion(Bead *bead1, iteratable &l);

  bool HasVel() { return _has_vel; }
  void SetHasVel(const bool v) { _has_vel = v; }

  bool HasForce() { return _has_force; }
  void SetHasForce(const bool v) { _has_force = v; }

  void setMoleculeNamesAndIds(
      std::map<std::string, int> molecule_name_and_ids) {
    assert(
        molecule_name_and_type_id_.size() == 0 &&
        "Cannot set the molecule "
        "names and ids in the topology object because the molecule names and "
        "ids are not empty");
    molecule_name_and_type_id_ = molecule_name_and_ids;
  }

  void setResidueIdsAndNames(
      std::map<int, std::set<std::pair<int, std::string>>>
          molecule_id_residue_name_and_ids) {
    assert(molecule_id_and_residue_id_and_name_.size() == 0 &&
           "Cannot set the "
           "molecules residue names and ids as it is not empty.");
    molecule_id_and_residue_id_and_name_ = molecule_id_residue_name_and_ids;
  }

  std::map<std::string, int> getMoleculeNamesAndIds() const {
    return molecule_name_and_type_id_;
  }

  std::map<int, std::set<std::pair<int, std::string>>> getResidueIdsAndNames()
      const {
    return molecule_id_and_residue_id_and_name_;
  }

 protected:
  BoundaryCondition *_bc;

  BoundaryCondition::eBoxtype autoDetectBoxType(const matrix &box);

  /// bead types in the topology
  std::map<std::string, int> beadtypes_;

  /// beads in the topology
  BeadContainer _beads;

  /// molecules in the topology
  MoleculeContainer _molecules;

  /// bonded interactions in the topology
  InteractionContainer _interactions;

  ExclusionList _exclusions;

  std::map<std::string, int> _interaction_groups;

  std::map<std::string, std::list<Interaction *>> _interactions_by_group;

  // Need some way to keep track of the unique residue ids , id of the molecule
  // type
  std::map<std::string, int> molecule_name_and_type_id_;
  std::map<int, std::set<std::pair<int, std::string>>>
      molecule_id_and_residue_id_and_name_;

  double _time;
  int _step;
  bool _has_vel;
  bool _has_force;
  int max_residue_id_;

  /// The particle group (For H5MD file format)
  std::string _particle_group;
};

template <class T>
inline T *Topology::CreateBead(byte_t symmetry, std::string name,
                               std::string type, int residue_number,
                               std::string residue_name,
                               std::string molecule_name, double m, double q) {

  int molecule_type_id;

  if (molecule_name_and_type_id_.count(molecule_name) == 0) {
    molecule_type_id = static_cast<int>(molecule_name_and_type_id_.size()) + 1;
    molecule_name_and_type_id_[molecule_name] = molecule_type_id;
  } else {
    molecule_type_id = molecule_name_and_type_id_[molecule_name];
  }

  std::pair<int, std::string> element{residue_number, residue_name};
  molecule_id_and_residue_id_and_name_[molecule_type_id].insert(element);

  T *bead = new T(this, _beads.size(), type, symmetry, name, residue_number,
                  residue_name, molecule_type_id, m, q);

  _beads.push_back(bead);
  return bead;
}

inline Molecule *Topology::CreateMolecule(std::string name) {
  Molecule *mol = new Molecule(this, _molecules.size(), name);
  _molecules.push_back(mol);
  return mol;
}

inline Molecule *Topology::MoleculeByIndex(int index) {
  return _molecules[index];
}

template <typename iteratable>
inline void Topology::InsertExclusion(Bead *bead1, iteratable &l) {
  _exclusions.InsertExclusion(bead1, l);
}

}  // namespace csg
}  // namespace votca

#include "interaction.h"

#endif /* _VOTCA_CSG_TOPOLOGY_H */
