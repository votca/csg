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

#ifndef VOTCA_CSG_FILEFORMATFACTORY_H
#define VOTCA_CSG_FILEFORMATFACTORY_H

#include <string>
#include <votca/tools/filesystem.h>
//#include <votca/tools/objectfactory.h>
#include "topologyobjectfactory.h"

namespace votca {
namespace csg {

template <class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
class FileFormatFactory : public csg::ObjectFactory<std::string, Bead_T,
                                                    Molecule_T, Topology_T, T> {
 public:
  FileFormatFactory() {}

  T<Bead_T, Molecule_T, Topology_T> *Create(const std::string &file);
};

template <class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
T<Bead_T, Molecule_T, Topology_T>
    *FileFormatFactory<Bead_T, Molecule_T, Topology_T, T>::Create(
        const std::string &file) {
  std::string filetype = tools::filesystem::GetFileExtension(file);
  try {
    return csg::ObjectFactory<std::string, Bead_T, Molecule_T, Topology_T,
                              T>::Create(filetype);
  } catch (std::exception &error) {
    throw std::runtime_error("Error '" + filetype +
                             "' file format of file "
                             "'" +
                             file + "' cannot be read or written");
  }
  return nullptr;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_FILEFORMATFACTORY_H
