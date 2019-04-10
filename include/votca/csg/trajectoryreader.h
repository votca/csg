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

//#include "basebead.h"
//#include "basemolecule.h"
//#include "pdbreader.h"
#include "fileformatfactory.h"
//#include "templatetopology.h"
#include <string>
namespace votca {
namespace csg {

/**
    \brief trajectoryreader interface

    This typename defines the interface a trajectory reader has to implement
 */
template <typename Bead_T, template <typename> class Molecule_T,
          template <typename, template <typename> class> class Topology_T>
class TrajectoryReader {
 private:
 public:
  virtual ~TrajectoryReader() {}
  /// open a trejectory file
  virtual bool Open(const std::string &file) = 0;

  virtual void Close() = 0;

  //  template<class Bead_T>;

  // template< template<class> class Bead_T>
  //  class Molecule_T;

  // virtual bool FirstFrame(TemplateTopology<BaseBead,BaseMolecule<BaseBead>>
  // &top)=0;
  virtual bool FirstFrame(void *top) = 0;
  /// read in the next frame
  // virtual bool NextFrame(TemplateTopology<BaseBead,BaseMolecule<BaseBead>>
  // *top) = 0;
  virtual bool NextFrame(void *top) = 0;
  // bool NextFrame(TemplateTopology<BaseBead,BaseMolecule<BaseBead>> &top);

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed

template <typename Bead_T, template <typename> class Molecule_T,
          template <typename, template <typename> class> class Topology_T,
          template <typename, template <typename> class,
                    template <typename, template <typename> class> class>
          class T>
inline FileFormatFactory<Bead_T, Molecule_T, Topology_T, T>
    &TrjReaderFactory() {
  static FileFormatFactory<Bead_T, Molecule_T, Topology_T, T> _TrjReaderFactory;
  return _TrjReaderFactory;
}

template <typename Bead_T, template <typename> class Molecule_T,
          template <typename, template <typename> class> class Topology_T>
class PDBReader;

template <typename Bead_T, template <typename> class Molecule_T,
          template <typename, template <typename> class> class Topology_T>
inline void TrajectoryReader<Bead_T, Molecule_T, Topology_T>::RegisterPlugins(
    void) {
  /*    TrjReaderFactory().Register<LAMMPSDumpReader>("dump");
      TrjReaderFactory().Register<LAMMPSDataReader>("data");
      TrjReaderFactory().Register<XYZReader>("xyz");
    #ifdef GMX_DOUBLE
      TrjReaderFactory().Register<GMXTrajectoryReader>("trr");
      TrjReaderFactory().Register<GMXTrajectoryReader>("xtc");
    #endif
      TrjReaderFactory().Register<GROReader>("gro");*/
  TrjReaderFactory<Bead_T, Molecule_T, Topology_T, TrajectoryReader>()
      .template Register<PDBReader<Bead_T, Molecule_T, Topology_T>>("pdb");
  /*
  TrjReaderFactory().Register<DLPOLYTrajectoryReader>("dlph");
  TrjReaderFactory().Register<DLPOLYTrajectoryReader>("dlpc");
#ifdef H5MD
  TrjReaderFactory().Register<H5MDTrajectoryReader>("h5");
#endif*/
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TRAJECTORYREADER_H
