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

#ifndef VOTCA_CSG_MAP_H
#define VOTCA_CSG_MAP_H

#include "cgbeadinfo.h"
#include "csgtopology.h"
#include "molecule.h"
#include <vector>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BoundaryCondition;
class BeadMap;

/*******************************************************
    Mapper class, collection of maps
*******************************************************/
class AtomisticToCGMoleculeMapper {
 public:
  AtomisticToCGMoleculeMapper(std::string atom_molecule_type,
                              std::string cg_molecule_type)
      : atom_molecule_type_(atom_molecule_type),
        cg_molecule_type_(cg_molecule_type){};
  ~AtomisticToCGMoleculeMapper();

  void Initialize(std::unordered_map<std::string, CGBeadInfo> bead_maps_info);

  // Pass in a map containing the names of all the atomistic beads in the molecule and pointers to them
  void Apply(CSG_Topology &atom_top,                                                     
     CSG_Topology& cg_top,                                                       
     pair<int,map<int,vector<pair<string,int>>>> cg_mol_id_cg_bead_id_atomic_bead_names_ids);


 protected:
  // Molecule atomistic_molecule_, cg_molecule_;
  std::string atom_molecule_type_;
  std::string cg_molecule_type_;
  std::unordered_map<string, std::unique_ptr<BeadMap>> bead_type_and_maps_;
};

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap {
 public:
  virtual ~BeadMap(){};
  virtual void Apply(const BoundaryCondition *boundaries,
                     std::map<std::string, Bead *> atomistic_beads, Bead *cg_bead) = 0;
  virtual void Initialize(std::vector<string> subbeads,
                          std::vector<double> weights,
                          std::vector<double> ds) = 0;
  virtual void Initialize(std::vector<string> subbeads,
                          std::vector<double> weights);

  std::vector<string> getAtomicBeadNames();

 protected:
  BeadMap(){};
  struct element_t {
    double weight_;
    double force_weight_;
  };
  std::unordered_map<std::string, element_t> matrix_;
  friend class AtomisticToCGMoleculeMapper;
};

/*******************************************************
    Linear map for spherical beads
*******************************************************/
class Map_Sphere : public BeadMap {
 public:
  void Apply(const BoundaryCondition *boundaries,
             std::map<std::string, Bead *> atomistic_beads, Bead *cg_bead);

  void Initialize(std::vector<std::string> subbeads, std::vector<double> weights, std::vector<double> ds); 

 protected:
  Map_Sphere() {}
  void AddElem(std::string atomic_bead_name, double weight, double force_weight);

  friend class AtomisticToCGMoleculeMapper;
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
             std::map<std::string, Bead *> atomistic_beads, Bead *cg_bead);

 protected:
  Map_Ellipsoid() {}
  friend class AtomisticToCGMoleculeMapper;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_MAP_H
