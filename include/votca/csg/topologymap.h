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

#ifndef _VOTCA_CSG_TOPOLOGYMAP_H
#define _VOTCA_CSG_TOPOLOGYMAP_H

#include "csgtopology.h"
#include "map.h"
#include <vector>

namespace votca {
namespace csg {

class TopologyMap {
 public:
  ~TopologyMap();

  TopologyMap(CSG_Topology *in, CSG_Topology *out);

  //  void AddMoleculeMap(Map *map);
  void LoadMap(std::string filename);

  void Apply();

 private:
  CSG_Topology *_in;
  CSG_Topology *_out;

  typedef std::vector<Map *> MapContainer;
  MapContainer _maps;
};

inline TopologyMap::TopologyMap(CSG_Topology *in, CSG_Topology *out)
    : _in(in), _out(out) {}

// inline void TopologyMap::AddMoleculeMap(Map *map) { _maps.push_back(map); }

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_TOPOLOGYMAP_H */
