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

#ifndef VOTCA_CSG_TOPOLOGYMAP_H
#define VOTCA_CSG_TOPOLOGYMAP_H
#include "atomcgconverter.h"
#include "csgtopology.h"
#include "map.h"
#include <vector>

namespace votca {
namespace csg {

/**
 * @brief Contains all mapping information for translating atomistic topology
 * data to cg topology data
 *
 * This consists of how to calculate the coarse grained positions,
 * velocities, and forces from the atomistic positions, velocities and
 * sources.
 */
class TopologyMap {
 public:
  ~TopologyMap();

  TopologyMap(CSG_Topology *atomistic_top, CSG_Topology *cg_top);

  //  void AddMoleculeMap(Map *map);
  //
  
  /**
   * @brief Will load a molecular coarse graining stencil from an .xml file 
   *
   * @param[in] filename
   */
  void LoadMap(const std::string & filename);

  /**
   * @brief Will map the positions, weights and d-values from an atomistic description to a coarse grained description.
   */
  void Apply(AtomCGConverter converter);

 private:
  CSG_Topology *atomistic_top_;
  CSG_Topology *cg_top_;

  std::unique_ptr<BoundaryCondition> boundaries_;
  // First - atomic molecule type
  // Second - cg molecule type
  typedef std::pair<std::string, std::string> AtomAndCGMoleculeTypes;
   

  /**
   * @brief Contains the atomic and coarse grained atom types, as well as the the mapping object that will map the atomic molecular type to a coarse graine type.
   */
  std::unordered_map<AtomAndCGMoleculeTypes, AtomisticToCGMoleculeMapper>
      molecule_names_and_maps_;
  // typedef std::vector<Map *> MapContainer;
  // MapContainer _maps;
  std::unordered_map<std::string,BeadMapInfo> ParseBeads_(TOOLS::Property & options);
  void ParseMaps_(TOOLS::Property & options, std::unordered_map<std::string, BeadMapInfo> & bead_maps_info);
};

inline TopologyMap::TopologyMap(CSG_Topology *atomistic_top,
                                CSG_Topology *cg_top)
    : atomistic_top_(atomistic_top), cg_top_(cg_top) {

  // Ensure that the boundary conditions are the same between the two topologies
  assert(atomistic_top->getBoxType() == cg_top->getBoxType() &&
         "Boundary conditions differ between atomistic and cg topologies, "
         "cannot create topology map.");

  boundaries_ = cg_top_->getBoundaryCondition()->Clone();
}

// inline void TopologyMap::AddMoleculeMap(Map *map) { _maps.push_back(map); }

}  // namespace csg
}  // namespace votca

#endif // VOTCA_CSG_TOPOLOGYMAP_H 
