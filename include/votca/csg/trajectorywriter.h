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
#ifndef VOTCA_CSG_TRAJECTORYWRITER_H
#define VOTCA_CSG_TRAJECTORYWRITER_H

#include "fileformatfactory.h"
#include <boost/any.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
namespace votca {
namespace csg {

class TrajectoryWriter {
 protected:
  std::ofstream out_;

 public:
  virtual ~TrajectoryWriter() {}

  virtual void Open(std::string file, bool bAppend = false) {
    std::cout << "Calling open" << std::endl;
    if (bAppend) {
      out_.open(file, std::ios_base::app);
    } else {
      std::cout << "And opening " << std::endl;
      out_.open(file);
    }
  }
  virtual void Close() { out_.close(); };

  virtual void Write(boost::any top) = 0;

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryWriter> &TrjWriterFactory() {
  static FileFormatFactory<TrajectoryWriter> _TrjWriterFactory;
  return _TrjWriterFactory;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TRAJECTORYWRITER_H
