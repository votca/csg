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

#include <votca/csg/topologymap.h>

namespace votca {
namespace csg {

TopologyMap::~TopologyMap() {
  MapContainer::iterator i;

  for (i = _maps.begin(); i != _maps.end(); ++i) delete *i;
  _maps.clear();
}

void TopologyMap::Apply() {
  MapContainer::iterator iter;

  _out->setStep(_in->getStep());
  _out->setTime(_in->getTime());
  _out->setBox(_in->getBox());

  for (iter = _maps.begin(); iter != _maps.end(); ++iter) (*iter)->Apply();
}

}  // namespace csg
}  // namespace votca
