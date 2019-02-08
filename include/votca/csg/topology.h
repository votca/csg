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
#include "topologytypecontainer.h"
#include "triclinicbox.h"

#include <votca/tools/matrix.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

  namespace topology_constants {
    const std::string unassigned_element = "unassigned";
  }

  namespace TOOLS = votca::tools;

  class Interaction;
  class ExclusionList;

  // typedef std::vector<Molecule *> MoleculeContainer;
  // typedef std::vector<Bead *> BeadContainer;
  typedef std::vector<Interaction *> InteractionContainer;

  /**
   * \brief topology of the whole system
   *
   * The Topology class stores the topology of the system like the beads, bonds,
   * molecules and residues.
   *
   **/
  template <class Bead_T, class Molecule_T>
  class Topology {
    public:
      /// constructor
      Topology() : _time(0.0), _has_vel(false), _has_force(false) {
        _bc = new OpenBox();
      }
      ~Topology() {
        Cleanup();
        if (_bc) delete (_bc);
        _bc = nullptr;
      }

      /**
       * \brief Cleans up all the stored data
       */
      void Cleanup() {
        // cleanup beads
        bead_.clear();
        // cleanup molecules
        molecules_.clear();
        // cleanup interactions
        {
          InteractionContainer::iterator i;
          for (i = _interactions.begin(); i < _interactions.end(); ++i) delete (*i);
          _interactions.clear();
        }
        // cleanup _bc object
        if (_bc) delete (_bc);
        _bc = new OpenBox();
      }

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
      Bead_T *CreateBead(byte_t symmetry, std::string bead_type, int bead_id,
          int molecule_id,  
          std::string residue_type, int residue_id, std::string element_name, 
          double mass,double charge) {

        assert(beads_.count(bead_id) == 0 &&
            "Cannot create bead with the provided id "
            "the id must be unique, and the one provided is already in use for "
            "this "
            "topology instance.");

        assert(
            type_container_.MoleculeTypeExist(molecule_type) == 0 &&
            "You must "
            "create the molecules before you can create a bead that is attached to "
            "a molecule");

        if (type_container_.ResidueTypeExist(residue_type) == 0) {
          type_container_.AddResidueType(residue_type);
        }
        if (type_container_.BeadTypeExist(bead_type) == 0) {
          type_container_.AddBeadType(bead_type);
        }
        beads_[bead_id] = Bead_T(bead_id, bead_type, symmetry, element_name,
            residue_id, molecule_id, mass, charge);

        //  beads_.push_back(bead);
        return &beads_[bead_id];
      }

      /**
       * \brief creates a new molecule
       * \param name name of the molecule
       * \return pointer to created molecule
       */
      Molecule_T *CreateMolecule(std::string name, int id) {
        assert(molecules_.count(id) == 0 &&
            "Trying to create a molecule that already "
            "exists the molecules id must be unique");
        //molecules_[id] = Molecule(this, id, name);
        molecules_[id] = Molecule(id, name);

        if (type_container_.MoleculeTypeExist(name) == 0) {
          type_container_.AddMoleculeType(name);
        }
        return &molecules_[id];
      }

      /**
       *  \brief checks weather molecules with the same name really contain the same
       * number of beads
       */
      void CheckMoleculeNaming(void) {
        std::unordered_map<std::string, size_t> molecule_type_and_bead_count;
        const std::unordered_set<std::string> molecule_types =
          type_container_.getMoleculeTypes();

        for (const std::pair<int, Molecule_T> &id_and_molecule : molecules_) {
          molecule_type = id_and_molecule.second.getType();
          if (molecule_type_and_bead_count.count(molecule_type)) {
            if (molecule_type_and_bead_count[molecule_type] !=
                id_and_molecule_.second.BeadCount()) {
              throw runtime_error(
                  "There are molecules which have the same name but different "
                  "number "
                  "of bead please check the section manual topology handling in "
                  "the "
                  "votca manual");
              continue;
            }
          }
          molecule_type_and_bead_count[molecule_type] =
            id_and_molecule_.second.BeadCount();
        }
      }

      /**
       * \brief create molecules based on blocks of atoms
       * \param[in] name molecule name
       * \param[in] first first bead
       * \param[in] nbeads number of beads per molecule
       * \param[in] nmolecules number of molecules
       */
      /*  void CreateMoleculesByRange(std::string name, int first, int nbeads,
          int nmolecules){

          Molecule_T *mol = CreateMolecule(name);
          int beadcount = 0;

          BeadContainer::iterator bead;
          for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
      // xml numbering starts with 1
      if (--first > 0) continue;
      // This is not 100% correct, but let's assume for now that the resnr do
      // increase
      mol->AddBead((*bead));

      if (++beadcount == nbeads) {
      if (--nmolecules <= 0) break;
      mol = CreateMolecule(name);
      beadcount = 0;
      }
      }
      }*/

      /**
       * \brief number of molecules in the system
       * @return number of molecule in topology
       */
      int MoleculeCount() { return molecules_.size(); }

      /**
       * number of beads in the system
       * @return number of beads in the system
       */
      int BeadCount() { return beads_.size(); }

      /**
       * get molecule by index
       * @param index molecule number
       * @return pointer to molecule
       */
      //Molecule_T *MoleculeByIndex(int index);
      Molecule_T * getMolecule(int molecule_id){
        assert(molecules_.count(molecule_id) && "Molecule with specified id"
            " does not exist within the topology instance.");
        return &(molecules_.at(molecule_id));
      }

      /**
       * access containter with all beads
       * @return bead container
       */
      BeadContainer &Beads() { return beads_; }

      /**
       * access  containter with all molecules
       * @return molecule container
       */
      MoleculeContainer &Molecules() { return molecules_; }

      /**
       * access containter with all bonded interactions
       * @return bonded interaction container
       */
      InteractionContainer &BondedInteractions() { return _interactions; }

      void AddBondedInteraction(Interaction *ic) {
        std::map<std::string, int>::iterator iter;
        iter = _interaction_groups.find(ic->getGroup());
        if (iter != _interaction_groups.end()) {
          ic->setGroupId((*iter).second);
        } else {
          int i = _interaction_groups.size();
          _interaction_groups[ic->getGroup()] = i;
          ic->setGroupId(i);
        }
        _interactions.push_back(ic);
        _interactions_by_group[ic->getGroup()].push_back(ic);
      }

      std::list<Interaction *> InteractionsInGroup(const std::string &group) {
        std::map<std::string, list<Interaction *>>::iterator iter;
        iter = _interactions_by_group.find(group);
        if (iter == _interactions_by_group.end()) return list<Interaction *>();
        return iter->second;
      }

      /**
       * \brief Determine if a bead type exists.
       *
       * @return bool true if it has been registered
       **/
      // bool BeadTypeExist(std::string type) const;

      /**
       * \brief Register the bead type with the topology object.
       *
       * Records are kept of the different bead types in the topology object. This
       * method stores the bead type.
       **/
      // void RegisterBeadType(std::string name);

      /**
       * \brief Given a bead type this method returns the id associated with the
       * type
       *
       * @param[in] std::string name of the type
       * @return int the id of the type
       **/
      std::string getBeadType(int bead_id) const {
        assert(beads_.count(bead_id) &&
            "Cannot get bead type as bead is not "
            "stored in the topology instance");
        return beads_.at(bead_id).getType();
      }

      /**
       * \brief Returns a pointer to the bead with index i
       *
       * @param[in] int index is the index of the bead in the internal bead
       * container
       * @return Bead * is a pointer to the bead
       **/
      Bead *getBead(const int bead_id) const { return beads_[bead_id]; }

      std::vector<int> getBeadIds() const {
        std::vector<int> bead_ids;
        for( const std::pair<const int,Bead_T> & id_and_bead : beads_ ){
          bead_ids.push_back(id_and_bead.first);
        }
        return bead_ids;
      }

      Molecule_T *getMolecule(const int i) const { return molecules_[i]; }

      /**
       * delete all molecule information
       */
      void ClearMoleculeList() { molecules_.clear(); }

      /**
       * \brief copy topology data of different topology
       *
       * WARNING does not copy interactions
       *
       * \param top topology to copy from
       */
      void CopyTopologyData(Topology<Bead_T,Molecule_T> *top) {
        BeadContainer::iterator it_bead;
        MoleculeContainer::iterator it_mol;

        _bc->setBox(top->getBox());
        _time = top->_time;
        _step = top->_step;

        // cleanup old data
        Cleanup();

        type_container_ = top->type_container_;
        // Copy residue info
        // setResidueIdsAndNames(top->getResidueIdsAndNames());
        // setMoleculeNamesAndIds(top->getMoleculeNamesAndIds());

        // create all beads
        /*for (it_bead = top->_beads.begin(); it_bead != top->_beads.end();
          ++it_bead) { Bead *bi = *it_bead; std::string type = bi->getType();
          std::string molecule_name =
          top->getMolecule(bi->getMoleculeId())->getName();
          CreateBead<Bead>(bi->getSymmetry(), bi->getName(), type,
          bi->getResidueNumber(), bi->getResidueName(),
          molecule_name, bi->getMass(), bi->getQ());
          }*/
        beads_ = top->beads_;

        // copy all molecules
        molecules_ = top->molecules_;
        /*    for (it_mol = top->_molecules.begin(); it_mol !=
              top->_molecules.end();
              ++it_mol) {
              Molecule *mi = CreateMolecule((*it_mol)->getName());
              vector<int> bead_ids = (*it_mol)->getBeadIds();
              for (const int &bead_id : bead_ids) {
              mi->AddBead(_beads[bead_id]);
              }
              }*/
      }

      /**
       *  \brief rename all the molecules in range
       * \param range range std::string of type 1:2:10 = 1, 3, 5, 7, ...
       * \param name new name of molecule
       * range is a std::string which is parsed by RangeParser,
       */
      void RenameMolecules(std::string range, std::string molecule_type) {
        rp.Parse(range);
        if (type_container_.MoleculeTypeExist(molecule_type) == false) {
          type_container_.AddMoleculeType(molecule_type);
        }
        for (i = rp.begin(); i != rp.end(); ++i) {
          if ((unsigned int)*i > _molecules.size()) {
            throw runtime_error(
                std::string("RenameMolecules: num molecules smaller than"));
          }
          //      getMolecule(*i - 1)->setName(name);
          molecules_[*i - 1].setType(molecule_type);
        }
      }

      /**
       *  \brief rename all the bead types
       * \param name current rame of the bead type
       * \param newname new name of bead type
       */
      void RenameBeadType(std::string old_bead_type, std::string new_bead_type) {
        for (std::pair<int, Bead_T> &id_and_bead : beads_) {
          std::string bead_type = id_and_bead.second.getType();
          if (wildcmp(old_bead_type.c_str(), bead_type.c_str())) {
            id_and_bead.second.setType(new_bead_type);
          }
        }
      }

      /**
       *  \brief set the mass of all the beads of a certain type
       * \param name the bead type
       * \param value mass value
       */
      void SetBeadTypeMass(std::string bead_type, double mass) {
        //  BeadContainer::iterator bead;
        //  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
        for (std::pair<int, Bead_T> &id_and_bead : beads_) {
          std::string type = id_and_bead.second.getType();
          if (wildcmp(name.c_str(), type.c_str())) {
            id_and_bead.second.setMass(mass);
          }
        }
      }

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
      }

      /**
       * get the simulation box
       * \return triclinic box matrix
       */
      const matrix &getBox() { return _bc->getBox(); }

      /**
       * set the time of current frame
       * \param t simulation time in ns
       */
      void setTime(double t) { _time = t; }

      /**
       * get the time of current frame
       * \return simulation time in ns
       */
      double getTime() { return _time; }

      /**
       * set the step number of current frame
       * \param s step number
       */
      void setStep(int s) { _step = s; }

      /**
       * get the step number of current frame
       * \return step number
       */
      int getStep() { return _step; }

      /**
       * Sets the particle group. (For the H5MD file format)
       * \param particle_group The name of a particle group.
       */
      void setParticleGroup(std::string particle_group) {
        _particle_group = particle_group;
      }

      /**
       * Gets the particle group.
       * \return The name of a particle group.
       */
      std::string getParticleGroup() { return _particle_group; }

      /**
       * \brief pbc correct distance of two beads
       * \param bead1 index of first bead
       * \param bead2 index of second bead
       * \return distance vector
       *
       * calculates the smallest distance between two beads with correct treatment
       * of pbc
       */
      TOOLS::vec getDist(int bead1, int bead2) const {
        return BCShortestConnection(getBead(bead1)->getPos(),
            getBead(bead2)->getPos());
      }

      /**
       * \brief calculate shortest vector connecting two points
       * \param r1 first point
       * \param r2 second point
       * \return distance vector
       *
       * calculates the smallest distance between two points with correct treatment
       * of pbc
       */
      TOOLS::vec BCShortestConnection(const vec &r1, const vec &r2) const {
        return _bc->BCShortestConnection(r_i, r_j);
      }

      /**
       * \brief return the shortest box size
       * \return shortest size
       *
       * Calculates the shortest length to connect two sides of the box
       */
      double ShortestBoxSize() {
        TOOLS::vec _box_a = getBox().getCol(0);
        TOOLS::vec _box_b = getBox().getCol(1);
        TOOLS::vec _box_c = getBox().getCol(2);

        // create plane normals
        TOOLS::vec _norm_a = _box_b ^ _box_c;
        TOOLS::vec _norm_b = _box_c ^ _box_a;
        TOOLS::vec _norm_c = _box_a ^ _box_b;

        _norm_a.normalize();
        _norm_b.normalize();
        _norm_c.normalize();

        double la = _box_a * _norm_a;
        double lb = _box_b * _norm_b;
        double lc = _box_c * _norm_c;

        return std::min(la, std::min(lb, lc));
      }

      /**
       *  calculates the box volume
       *  \return box volume
       */
      double BoxVolume() { return _bc->BoxVolume(); }

      /**
       *  rebuild exclusion list
       */
      void RebuildExclusions() { _exclusions.CreateExclusions(this); }

      /**
       * access exclusion list
       * \return exclusion list
       */
      ExclusionList &getExclusions() { return _exclusions; }

      BoundaryCondition::eBoxtype getBoxType() { return _bc->getBoxType(); }

      template <typename iteratable>
        void InsertExclusion(Bead *bead1, iteratable &l) {
          _exclusions.InsertExclusion(bead1, l);
        }

      bool HasVel() { return _has_vel; }
      void SetHasVel(const bool v) { _has_vel = v; }

      bool HasForce() { return _has_force; }
      void SetHasForce(const bool v) { _has_force = v; }
      /*
         void setMoleculeNamesAndIds(
         std::map<std::string, int> molecule_name_and_ids) {
         assert(
         molecule_name_and_type_id_.size() == 0 &&
         "Cannot set the molecule "
         "names and ids in the topology object because the molecule names and "
         "ids are not empty");
         molecule_name_and_type_id_ = molecule_name_and_ids;
         }*/
      /*
         void setResidueIdsAndNames(
         std::map<int, std::set<std::pair<int, std::string>>>
         molecule_id_residue_name_and_ids) {
         assert(molecule_id_and_residue_id_and_name_.size() == 0 &&
         "Cannot set the "
         "molecules residue names and ids as it is not empty.");
         molecule_id_and_residue_id_and_name_ = molecule_id_residue_name_and_ids;
         }
         */
      /*  std::map<std::string, int> getMoleculeNamesAndIds() const {
          return molecule_name_and_type_id_;
          }

          std::map<int, std::set<std::pair<int, std::string>>> getResidueIdsAndNames()
          const {
          return molecule_id_and_residue_id_and_name_;
          }*/

    protected:
      BoundaryCondition *_bc;

      BoundaryCondition::eBoxtype autoDetectBoxType(const matrix &box) {
        // set the box type to OpenBox in case "box" is the zero matrix,
        // to OrthorhombicBox in case "box" is a diagonal matrix,
        // or to TriclinicBox otherwise
        if (box.get(0, 0) == 0 && box.get(0, 1) == 0 && box.get(0, 2) == 0 &&
            box.get(1, 0) == 0 && box.get(1, 1) == 0 && box.get(1, 2) == 0 &&
            box.get(2, 0) == 0 && box.get(2, 1) == 0 && box.get(2, 2) == 0) {
          // cout << "box open\n";
          return BoundaryCondition::typeOpen;
        } else if (box.get(0, 1) == 0 && box.get(0, 2) == 0 && box.get(1, 0) == 0 &&
            box.get(1, 2) == 0 && box.get(2, 0) == 0 && box.get(2, 1) == 0) {
          // cout << "box orth\n";
          return BoundaryCondition::typeOrthorhombic;
        } else {
          // cout << "box tric\n";
          return BoundaryCondition::typeTriclinic;
        }
        return BoundaryCondition::typeOpen;
      }

      /// bead types in the topology
      // std::map<std::string, int> beadtypes_;

      /// beads in the topology
      std::map<int, Bead_T> beads_;

      /// molecules in the topology
      std::map<int, Molecule_T> molecules_;

      /// bonded interactions in the topology
      InteractionContainer _interactions;

      ExclusionList _exclusions;

      std::map<std::string, int> _interaction_groups;

      std::map<std::string, std::list<Interaction *>> _interactions_by_group;

      // Need some way to keep track of the unique residue ids , id of the molecule
      // type
      // std::map<std::string, int> molecule_name_and_type_id_;
      // std::map<int, std::set<std::pair<int, std::string>>>
      //    molecule_id_and_residue_id_and_name_;

      TopologyTypeContainer type_container_;

      double _time;
      int _step;
      bool _has_vel;
      bool _has_force;
      // int max_residue_id_;

      /// The particle group (For H5MD file format)
      std::string _particle_group;
      };

      /**
       * @brief Create Beads
       *
       * The type is the same as the name
       * Each bead should have
       * Molecule name
       * Molecule id - locally unique to topology object
       * Residue name
       * Residue number/id - locally unique to topology object
       * Bead name/Atom name
       * Bead/Atom id - locally unique to topology object
       * Element name
       *
       * @tparam T
       * @param symmetry
       * @param name
       * @param type
       * @param residue_number
       * @param residue_name
       * @param molecule_name
       * @param m
       * @param q
       *
       * @return
       */
      /*template <class Bead_T, class Molecule_T>
        Bead_T * Topology::CreateBead(TOOLS::byte_t symmetry, std::string bead_name,int
        bead_id, int molecule_id, // Don't need the molecule name that is stored when
        creating molecules std::string residue_type, int residue_id, double mass, double
        charge) {

        assert(beads_.count(bead_id)==0 && "Cannot create bead with the provided id "
        "the id must be unique, and the one provided is already in use for this "
        "topology instance.");

        assert(type_container_.MoleculeTypeExist(molecule_type)==0 && "You must "
        "create the molecules before you can create a bead that is attached to "
        "a molecule");

        if(type_container_.ResidueTypeExist(residue_type)==0){
        type_container_.AddResidueType(residue_type);
        }

        beads_[bead_id] = Bead_T(bead_id,
        bead_name,
        symmetry,
        element_name,
        residue_id,
        molecule_id,
        mass,
        charge);

      //  beads_.push_back(bead);
      return &beads_[bead_id];
      }*/
      /*
         inline Molecule *Topology::CreateMolecule(std::string name,int id) {
         assert(molecules_.count(id)==0 && "Trying to create a molecule that already "
         "exists the molecules id must be unique");
         molecules_[id] = Molecule(this,id, name);

         if(type_container_.MoleculeTypeExist(name)==0){
         type_container_.AddMoleculeType(name);
         }
         return &molecules_[id];
         }
         */
      /*inline Molecule *Topology::MoleculeByIndex(int index) {
        return molecules_[index];
        }*/

      /*template <typename iteratable>
        inline void Topology::InsertExclusion(Bead *bead1, iteratable &l) {
        _exclusions.InsertExclusion(bead1, l);
        }*/

    }  // namespace csg
}  // namespace votca

  //#include "interaction.h"
#endif /* _VOTCA_CSG_TOPOLOGY_H */
