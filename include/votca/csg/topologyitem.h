/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _TOPOLOGYITEM_H
#define	_TOPOLOGYITEM_H

#include <memory>

namespace votca { namespace csg {

class Topology;

class TopologyItem
{
public:    
    virtual ~TopologyItem() {}
    std::shared_ptr<Topology> getParent() { return _parent; }
protected:
    TopologyItem(std::shared_ptr<Topology> parent)
        : _parent(parent) {}
    
    std::shared_ptr<Topology> _parent;
    
    friend class Topology;
};

}}

#endif	/* _TOPOLOGYITEM_H */

