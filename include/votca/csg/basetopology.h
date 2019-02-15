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
#include "topologytypecontainer.h"

#include <votca/tools/matrix.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

template<class Bead_T, class Molecule_T>
class Topology {
 public:
  /// constructor
  Topology()
      : time_(0.0),
        has_vel_(false),
        has_force_(false),
         {
    bc_ = OpenBox();
  }
  ~Topology();

  /**
   * \brief Cleans up all the stored data
   */
  void Cleanup();

  /**
   *  \brief checks weather molecules with the same name really contain the same
   * number of beads
   */
  void CheckMoleculeNaming(void);

	/**
   * \brief number of molecules in the system
   * @return number of molecule in topology
   */
  size_t MoleculeCount() { return molecules_.size(); }

  /**
   * number of beads in the system
   * @return number of beads in the system
   */
  size_t BeadCount() { return beads_.size(); }


  /**
   * access containter with all beads
   * @return bead container
   */
  //BeadContainer &Beads() { return beads_; }

  /**
   * access  containter with all molecules
   * @return molecule container
   */
  //MoleculeContainer &Molecules() { return molecules_; }

  /**
   * access containter with all bonded interactions
   * @return bonded interaction container
   */
  InteractionContainer &BondedInteractions() { return interactions_; }

  void AddBondedInteraction(Interaction *ic);
  std::list<Interaction *> InteractionsInGroup(const std::string &group);

  /**
   * \brief Determine if a bead type exists.
   *
   * @return bool true if it has been registered
   **/
  //bool BeadTypeExist(std::string type) const;

  /**
   * \brief Register the bead type with the topology object.
   *
   * Records are kept of the different bead types in the topology object. This
   * method stores the bead type.
   **/
  //void RegisterBeadType(std::string name);

  /**
   * \brief Given a bead type this method returns the id associated with the
   * type
   *
   * @param[in] string name of the type
   * @return int the id of the type
   **/
  //int getBeadTypeId(std::string type) const;

  /**
   * \brief Returns a pointer to the bead with index i
   *
   * @param[in] int index is the index of the bead in the internal bead
   * container
   * @return Bead * is a pointer to the bead
   **/
  Bead *getBead(const int id) const { return beads_.at(id); }
  Molecule *getMolecule(const int id) const { return molecules_.at(id); }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() { molecules_.clear(); }

  /**
   * \brief copy topology data of different topology
   * 
	 * Copies everything but the interactions
	 * \param top topology to copy from
   */
  void CopyTopologyData(Topology *top);

  /**
   *  \brief rename all the molecules in range
   * \param range range string of type 1:2:10 = 1, 3, 5, 7, ...
   * \param name new name of molecule
   * range is a string which is parsed by RangeParser,
   */
  void RenameMoleculesType(std::string range, std::string name);

  /**
   *  \brief rename all the bead types
   * \param name current rame of the bead type
   * \param newname new name of bead type
   */
  void RenameBeadsType(std::string name, std::string newname);

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
  void setBox(const matrix &box, BoundaryCondition::eBoxtype boxtype);

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const matrix &getBox() { return bc_->getBox(); }

  /**
   * set the time of current frame
   * \param t simulation time in ns
   */
  void setTime(double t) { time_ = t; }

  /**
   * get the time of current frame
   * \return simulation time in ns
   */
  double getTime() { return time_; }

  /**
   * set the step number of current frame
   * \param s step number
   */
  void setStep(int s) { step_ = s; }

  /**
   * get the step number of current frame
   * \return step number
   */
  int getStep() { return step_; };

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
  std::string getParticleGroup() { return particle_group_; }

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

  BoundaryCondition::eBoxtype getBoxType() { return bc_.getBoxType(); }

  template <typename iteratable>
  void InsertExclusion(Bead *bead1, iteratable &l);

  bool HasVel() { return has_vel_; }
  void SetHasVel(const bool v) { has_vel_ = v; }

  bool HasForce() { return has_force_; }
  void SetHasForce(const bool v) { has_force_ = v; }

 	protected:
  BoundaryCondition bc_;

  BoundaryCondition::eBoxtype autoDetectBoxType(const matrix &box);

  /// bead types in the topology
  //std::map<std::string, int> beadtypes_;

  /// beads in the topology
	std::unordered_map<int,Bead_T> beads_;

  /// molecules in the topology
	std::unordered_map<int,Molecule_T> molecules_;

  /// bonded interactions in the topology
  InteractionContainer interactions_;

  ExclusionList _exclusions;

  std::map<std::string, int> _interaction_groups;

  std::map<std::string, std::list<Interaction *>> interactions__by_group;

  double time_;
  int step_;
  bool has_vel_;
  bool has_force_;

  /// The particle group (For H5MD file format)
  std::string particle_group_;

	TopologyTypeContainer type_container_;
};

template <class Bead_T, class Molecule_T>
template <typename iteratable>
void Topology<Bead_T,Molecule_T>::InsertExclusion(Bead_T *bead1, iteratable &l) {
  _exclusions.InsertExclusion(bead1, l);
}

template <class Bead_T, class Molecule_T>
void Topology<Bead_T,Molecule_T>::Cleanup() {
		
		bead_.clear();
		molecule_.clear();	
   
	 	type_container.Clear();	
       
		// cleanup interactions
    {
      InteractionContainer::iterator i;
      for (i = interactions_.begin(); i < interactions_.end(); ++i) delete (*i);
      interactions_.clear();
    }
    // cleanup bc_ object
    bc_ = OpenBox();
  }

template<class Bead_T, class Molecule_T>
void Topology<Bead_T,Molecule_T>::CopyTopologyData(const Topology & top) {
	step_ = top.step_;
	time_ = top.time_;
 	has_vel_ = top.has_vel_;
	has_force_ = top.has_force_;	
	beads_ = top.beads_;
	molecules_ = top.molecules_;
	bc_ = top.bc_;
	type_container_ = top.type_container_;

}

template<class Bead_T, class Molecule_T>
void Topology<Bead_T,Molecule_T>::setBox(const matrix &box, BoundaryCondition::eBoxtype boxtype =
                                     BoundaryCondition::typeAuto) {
    // determine box type automatically in case boxtype==typeAuto
    if (boxtype == BoundaryCondition::typeAuto) {
      boxtype = autoDetectBoxType(box);
    }

    switch (boxtype) {
      case BoundaryCondition::typeTriclinic:
        bc_ = TriclinicBox();
        break;
      case BoundaryCondition::typeOrthorhombic:
        bc_ = OrthorhombicBox();
        break;
      default:
        bc_ = OpenBox();
        break;
    }

    bc_->setBox(box);
  }

template<class Bead_T, class Molecule_T>
  void Topology<Bead_T,Molecule_T>::RenameMoleculesType(string range_molecule_ids, string type) {
    RangeParser rp;
    RangeParser::iterator molecule_id_ptr;
  
    rp.Parse(range_molecule_ids);
    for (molecule_id_ptr = rp.begin(); molecule_id_ptr != rp.end(); ++molecule_id_ptr) {
      if ((unsigned int)*molecule_id_ptr > _molecules.size()) {
        throw runtime_error(
            string("RenameMoleculesType: num molecules smaller than"));
      }
      getMolecule(*molecule_id_ptr)->setType(type);
    }
  }

template<class Bead_T, class Molecule_T>
  void Topology<Bead_T,Molecule_T>::RenameBeadsType(string old_type, string new_type) {
 		for( pair<const int, Bead_T> & id_and_bead : beads_){
			string bead_type = id_and_bead.second.getType();
			if (wildcmp(bead_type.c_str(), old_type.c_str())) {
				id_and_bead.second.setType(new_type);
			}
		}	
  }

template<class Bead_T, class Molecule_T>
  void Topology<Bead_T,Molecule_T>::setBeadOfGivenTypeToNewMass(string type, double mass) {

		for( pair<const int, Bead_T> & id_and_bead : beads_ ){
			string bead_type = id_and_bead.second.getType();
			if( wildcmp(bead_type.c_str(), type.c_str())){
				id_and_bead.second.setMass(mass);
			}
		}	
  }

template<class Bead_T,class Molecule_T>
void Topology<Bead_T,Molecule_T>::CheckMoleculeNaming(void) {
	std::unordered_map<std::string, size_t> number_of_beads_in_each_molecular_type;

	for ( const std::pair<int,Molecule_T> & id_and_molecule : molecules_){
		std::string type = id_and_molecule.second.getType();
	
	}
	for (MoleculeContainer::iterator iter = _molecules.begin();
			iter != _molecules.end(); ++iter) {
		map<string, int>::iterator entry = nbeads.find((*iter)->getName());
		if (entry != nbeads.end()) {
			if (entry->second != static_cast<int>((*iter)->BeadCount()))
				throw runtime_error(
						"There are molecules which have the same name but different number "
						"of bead please check the section manual topology handling in the "
						"votca manual");
			continue;
		}
		nbeads[(*iter)->getName()] = (*iter)->BeadCount();
	}
}



} // namespace csg
}  // namespace votca

A
#endif /* _VOTCA_CSG_TOPOLOGY_H */
