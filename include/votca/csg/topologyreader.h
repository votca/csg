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

#ifndef _VOTCA_CSG_TOPOLOGYREADER_H
#define _VOTCA_CSG_TOPOLOGYREADER_H

#include "basebead.h"
#include "basemolecule.h"
#include "fileformatfactory.h"
#include "templatetopology.h"
#include <string>

namespace votca {
namespace csg {

class TopologyReader {
 public:
  virtual ~TopologyReader() {}
  /// open, read and close topology file
  // virtual bool ReadTopology(std::string file,
  // TemplateTopology<BaseBead,BaseMolecule<BaseBead>> &top) = 0;

  template <class Bead_T, class Molecule_T>
  bool ReadTopology(std::string file,
                    TemplateTopology<Bead_T, Molecule_T> &top) {
    throw std::runtime_error(
        "ReadTopology method must be derived by child class.");
  }

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TopologyReader> &TopReaderFactory() {
  static FileFormatFactory<TopologyReader> _TopReaderFactory;
  return _TopReaderFactory;
}
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_TOPOLOGYREADER_H
