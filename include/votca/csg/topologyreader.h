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
#include <boost/any.hpp>
#include <string>

namespace votca {
namespace csg {

class TopologyReader {
 private:
  /*  template <class Bead_T, class Molecule_T>
      bool ReadTopology_(std::string file,
          TemplateTopology<Bead_T, Molecule_T> &top) {
        throw std::runtime_error(
            "ReadTopology_ method must be defined by child class.");
      }*/
 public:
  virtual ~TopologyReader() {}
  /// open, read and close topology file
  // virtual bool ReadTopology(std::string file,
  // TemplateTopology<BaseBead,BaseMolecule<BaseBead>> &top) = 0;

  //  template <class Bead_T, class Molecule_T>
  //  virtual bool ReadTopology(std::string file,
  //                   TemplateTopology<BaseBead, BaseMolecule<BaseBead>> * top)
  //                   = 0;
  virtual bool ReadTopology(std::string file, boost::any top) { return false; }

  //    return ReadTopology_(file,top);
  // }
  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
/*template<class Bead_T,
         template <class > class Molecule_T,
         template < class, template< class > class > class Topology_T> */
/*csg::FileFormatFactory<Bead_T,Molecule_T,Topology_T,TopologyReader> &
TopReaderFactory() { static
FileFormatFactory<Bead_T,Molecule_T,Topology_T,TopologyReader>
_TopReaderFactory; return _TopReaderFactory;
}*/

inline FileFormatFactory<TopologyReader>& TopReaderFactory() {
  static FileFormatFactory<TopologyReader> _TopReaderFactory;
  return _TopReaderFactory;
}

// template <typename Bead_T,  class Molecule_T, class Topology_T>
// class PDBReader;

// template <class Bead_T,  class Molecule_T, class Topology_T>
// void TopologyReader<Bead_T, Molecule_T, Topology_T>::RegisterPlugins(void) {
//  TopReaderFactory().Register<XMLTopologyReader>("xml");
//  TopReaderFactory().Register<LAMMPSDumpReader>("dump");
//  TopReaderFactory().Register<LAMMPSDataReader>("data");
//  TopReaderFactory().Register<XYZReader>("xyz");
//  TopReaderFactory().Register<GROReader>("gro");
//#ifdef GMX_DOUBLE
//  TopReaderFactory().Register<GMXTopologyReader>("tpr");
//#endif
// TopReaderFactory<Bead_T, Molecule_T, Topology_T>()
//    .template Register<PDBReader<Bead_T, Molecule_T, Topology_T>>("pdb");
//  TopReaderFactory().Register<DLPOLYTopologyReader>("dlpf");
//}

}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_TOPOLOGYREADER_H
