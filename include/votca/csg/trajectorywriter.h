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

#ifndef VOTCA_CSG_TRAJECTORYWRITER_H
#define VOTCA_CSG_TRAJECTORYWRITER_H

//#include "basemolecule.h"
#include "fileformatfactory.h"
//#include "templatetopology.h"
#include <boost/any.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
namespace votca {
namespace csg {

class TrajectoryWriter {
 private:
  /*    template <class Bead_T, class Molecule_T>
        void Write_(TemplateTopology<Bead_T, Molecule_T> *top) {
          throw std::runtime_error(
              "Trajectory Write method must be derived by child class");
        }*/
 public:
  TrajectoryWriter() {}
  virtual ~TrajectoryWriter() {}

  virtual void Open(std::string file, bool bAppend = false) {}
  virtual void Close(){};

  // void Write(CSG_Topology *top) {};
  // template <class Bead_T, class Molecule_T>
  virtual void Write(boost::any top) = 0;

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryWriter> &TrjWriterFactory() {
  static FileFormatFactory<TrajectoryWriter> _TrjWriterFactory;
  return _TrjWriterFactory;
}

// template <class Bead_T,  class Molecule_T,class Topology_T>
// void TrajectoryWriter<Bead_T, Molecule_T, Topology_T>::RegisterPlugins() {

// j  TrjWriterFactory<Bead_T,Molecule_T,Topology_T>()
//    .template Register<PDBWriter<Bead_T,Molecule_T,Topology_T>>("pdb");
/*  TrjWriterFactory().Register<XYZWriter>("xyz");
  TrjWriterFactory().Register<LAMMPSDumpWriter>("dump");
  TrjWriterFactory().Register<DLPOLYTrajectoryWriter>("dlph");
  TrjWriterFactory().Register<DLPOLYTrajectoryWriter>("dlpc");
#ifdef GMX_DOUBLE
  TrjWriterFactory().Register<GMXTrajectoryWriter>("trr");
  TrjWriterFactory().Register<GMXTrajectoryWriter>("xtc");
#endif
  TrjWriterFactory().Register<GROWriter>("gro"); */
//}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TRAJECTORYWRITER_H
