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

#pragma once
#ifndef VOTCA_CSG_TOPOLOGYREADER_H
#define VOTCA_CSG_TOPOLOGYREADER_H

#include "basebead.h"
#include "basemolecule.h"
#include "fileformatfactory.h"
#include "templatetopology.h"
#include <boost/any.hpp>
#include <string>

namespace votca {
namespace csg {

/**
 * @brief Topology Reader reads the topology
 *
 * NOTE The topology class cannot be made into a pure abstract class because
 * it is instantiated by a factory method with a call to new and an abstract
 * class is not allowed to be instantiated. The factory method makes use of
 * the polymorphic behavior which is why it is needed.
 */

class TopologyReader {
 public:
  virtual ~TopologyReader() {}

  virtual bool ReadTopology(std::string file, boost::any top) { return false; }

  static void RegisterPlugins(void);
};

inline FileFormatFactory<TopologyReader>& TopReaderFactory() {
  static FileFormatFactory<TopologyReader> _TopReaderFactory;
  return _TopReaderFactory;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TOPOLOGYREADER_H
