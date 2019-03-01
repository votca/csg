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
  AtomisticToCGMoleculeMapper(const Molecule &atomistic_molecule,
                              Molecule &cg_molecule)
      : atomistic_molecule_(atomistic_molecule), cg_molecule_(cg_molecule) {}
  ~AtomisticToCGMoleculeMapper();

  void Apply();

  BeadMap *CreateBeadMap(const byte_t symmetry,
                         const BoundaryCondition *boundaries,
                         const Molecule *atomistic_molecule, Bead *bead_out,
                         Property *opts_map, Property *opts_bead);

 protected:
  Molecule atomistic_molecule_, cg_molecule_;
  std::vector<std::unique_ptr<BeadMap>> bead_maps_;
};

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap {
 public:
  virtual ~BeadMap(){};
  virtual void Apply() = 0;
  virtual void Initialize(const BoundaryCondition *boundaries,
                          const Molecule *mol_in, Bead *bead_out,
                          Property *opts_map, Property *opts_bead);

 protected:
  BeadMap(){};
  const BoundaryCondition *boundaries_;
  const Molecule *mol_in_;
  Bead *bead_out_;
  Property *opts_map_;
  Property *opts_bead_;
  friend class AtomisticToCGMoleculeMapper;
};

inline void BeadMap::Initialize(const BoundaryCondition *boundaries,
                                const Molecule *mol_in, Bead *bead_out,
                                Property *opts_bead, Property *opts_map) {
  boundaries_ = boundaries;
  mol_in_ = mol_in;
  bead_out_ = bead_out;
  opts_map_ = opts_map;
  opts_bead_ = opts_bead;
}

/*******************************************************
    Linear map for spherical beads
*******************************************************/
class Map_Sphere : public BeadMap {
 public:
  void Apply();

  void Initialize(const BoundaryCondition *boundaries, const Molecule *mol_in,
                  Bead *bead_out, Property *opts_bead, Property *opts_map);

 protected:
  Map_Sphere() {}
  void AddElem(const Bead *bead_in, double weight, double force_weight);

  struct element_t {
    const Bead *bead_in_;
    double weight_;
    double force_weight_;
  };
  std::vector<element_t> matrix_;
  friend class AtomisticToCGMoleculeMapper;
};

inline void Map_Sphere::AddElem(const Bead *bead_in, double weight,
                                double force_weight) {
  element_t el;
  el.bead_in_ = bead_in;
  el.weight_ = weight;
  el.force_weight_ = force_weight;
  matrix_.push_back(el);
}

/*******************************************************
    Linear map for ellipsoidal bead
*******************************************************/
class Map_Ellipsoid : public Map_Sphere {
 public:
  void Apply();

 protected:
  Map_Ellipsoid() {}
  friend class AtomisticToCGMoleculeMapper;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_MAP_H
