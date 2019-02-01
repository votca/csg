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
#include <assert.h>
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
    \brief topology of the whole system

    The Topology class stores the topology of the system like the beads, bonds,
 molecules and residues.

    \todo internal management for ids and indices
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
   * \brief cleans up all the stored data
   */
  virtual void Cleanup();

  /**
   * \brief creates a new Bead
   *
   * \param symmetry symmetry of the bead, 1: spherical 3: ellipsoidal
   * \param name name of the bead
   * \param type bead type
   * \param residue_number residue number
   * \param residue_name residue name
   * \param m mass
   * \param q charge
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
   * \brief put the whole topology in one molecule
   * \param name name of the new molecule
   *
   *  This function creates one big molecule for all beads in the topology.
   */
  void CreateOneBigMolecule(std::string name);

  /**
   * \brief create molecules based on blocks of atoms
   * \param name molecule name
   * \param first first bead
   * \param nbeads number of beads per molecule
   * \param nmolecules number of molecules
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

  bool BeadTypeExist(std::string type) const;
  void RegisterBeadType(std::string name);

  int getBeadTypeId(std::string type) const;
  Bead *getBead(const int i) const { return _beads[i]; }
  Molecule *getMolecule(const int i) const { return _molecules[i]; }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() { _molecules.clear(); }

  /**
   * \brief adds all the beads+molecules+residues from other topology
   * \param top topology to add
   */
  void Add(Topology *top);

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

  std::map<std::string, std::pair<int, std::map<int, std::string>>>
      getResidueIdsAndNames() const {
    return moleculename_residue_ids_and_names_;
  }

  //  int getMaxResidueId() const { return max_residue_id_; }

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
  std::map<std::string, std::pair<int, std::map<int, std::string>>>
      moleculename_residue_ids_and_names_;
  int max_residue_id_;
  double _time;
  int _step;
  bool _has_vel;
  bool _has_force;

  /// The particle group (For H5MD file format)
  std::string _particle_group;
};

template <class T>
inline T *Topology::CreateBead(byte_t symmetry, std::string name,
                               std::string type, int residue_number,
                               std::string residue_name,
                               std::string molecule_name, double m, double q) {

  int molecule_type_id;
  if (moleculename_residue_ids_and_names_.count(molecule_name) == 0) {
    molecule_type_id =
        static_cast<int>(moleculename_residue_ids_and_names_.size()) + 1;
  } else {
    molecule_type_id = moleculename_residue_ids_and_names_[molecule_name].first;
  }

  std::map<int, std::string> res_number_and_name;
  res_number_and_name[residue_number] = residue_name;
  moleculename_residue_ids_and_names_[molecule_name] =
      std::pair<int, std::map<int, std::string>>(molecule_type_id,
                                                 res_number_and_name);

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
