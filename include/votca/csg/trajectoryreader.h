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

#ifndef VOTCA_CSG_TRAJECTORYREADER_H
#define VOTCA_CSG_TRAJECTORYREADER_H

#include "basebead.h"
#include "basemolecule.h"
#include "fileformatfactory.h"
#include "templatetopology.h"
#include <string>

namespace votca {
namespace csg {

/**
    \brief trajectoryreader interface

    This class defines the interface a trajectory reader has to implement
 */
class TrajectoryReader {
 public:
  virtual ~TrajectoryReader() {}
  /// open a trejectory file
  virtual bool Open(const std::string &file) = 0;

  virtual void Close(){};

  //  template<class Bead_T>;

  // template< template<class> class Bead_T>
  //  class Molecule_T;

  template <class Bead_T, class Molecule_T>
  bool FirstFrame(TemplateTopology<Bead_T, Molecule_T> &top) {
    throw std::runtime_error(
        "The FirstFrame method must be overwritten by a child class");
  }
  /// read in the next frame
  template <class Bead_T, class Molecule_T>
  bool NextFrame(TemplateTopology<Bead_T, Molecule_T> &top) {
    throw std::runtime_error(
        "The NextFrame method must be overwritten by a child class");
  }
  // bool NextFrame(TemplateTopology<BaseBead,BaseMolecule<BaseBead>> &top);

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryReader> &TrjReaderFactory() {
  static FileFormatFactory<TrajectoryReader> _TrjReaderFactory;
  return _TrjReaderFactory;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TRAJECTORYREADER_H
