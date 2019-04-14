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

#include "fileformatfactory.h"
#include <boost/any.hpp>
#include <string>
namespace votca {
namespace csg {

/**
    \brief trajectoryreader interface

    This typename defines the interface a trajectory reader has to implement
 */
class TrajectoryReader {
 private:
 public:
  virtual ~TrajectoryReader() {}
  /// open a trajectory file
  virtual bool Open(const std::string &file){};

  virtual void Close(){};

  virtual bool FirstFrame(boost::any top) {
    throw std::runtime_error(
        "Cannot call FirstFrame using Trajectory reader it must be derived by "
        "child class.");
  }

  virtual bool NextFrame(boost::any top) {
    throw std::runtime_error(
        "Cannot call NextFrame using Trajectory reader it must be derived by "
        "child class.");
  }

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
