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

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../../include/votca/csg/csgtopology.h"
#include "../../include/votca/csg/io/groreader.h"
#include "../../include/votca/csg/io/lammpsdumpreader.h"
#include "../../include/votca/csg/io/pdbreader.h"
#include "../../include/votca/csg/io/xmltopologyreader.h"
#include "../../include/votca/csg/io/xyzreader.h"
#include "../../include/votca/csg/topologyreader.h"
#ifdef GMX_DOUBLE
#include "../../include/votca/csg/io/gmxtopologyreader.h"
#endif
#include "../../include/votca/csg/io/dlpolytopologyreader.h"
#include "../../include/votca/csg/io/lammpsdatareader.h"
namespace votca {
namespace csg {

void TopologyReader::RegisterPlugins(void) {
  TopReaderFactory().Register<PDBReader<CSG_Topology>>("pdb");
  TopReaderFactory().Register<XMLTopologyReader<CSG_Topology>>("xml");
  TopReaderFactory().Register<XYZReader<CSG_Topology>>("xyz");
  TopReaderFactory().Register<LAMMPSDumpReader<CSG_Topology>>("dump");
  TopReaderFactory().Register<GROReader<CSG_Topology>>("gro");
#ifdef GMX_DOUBLE
  TopReaderFactory().Register<GMXTopologyReader<CSG_Topology>>("tpr");
#endif
  TopReaderFactory().Register<LAMMPSDataReader<CSG_Topology>>("data");
  TopReaderFactory().Register<DLPOLYTopologyReader<CSG_Topology>>("dlpf");
}

}  // namespace csg
}  // namespace votca
