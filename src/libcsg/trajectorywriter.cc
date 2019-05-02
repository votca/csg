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
#include "../../include/votca/csg/io/growriter.h"
#include "../../include/votca/csg/io/pdbwriter.h"
#include "../../include/votca/csg/io/xyzwriter.h"
#include "../../include/votca/csg/trajectorywriter.h"
#ifdef GMX_DOUBLE
#include "../../include/votca/csg/io/gmxtrajectorywriter.h"
#endif
#include "../../include/votca/csg/io/dlpolytrajectorywriter.h"
#include "../../include/votca/csg/io/lammpsdumpwriter.h"

namespace votca {
namespace csg {

using namespace std;
void TrajectoryWriter::RegisterPlugins() {
  TrjWriterFactory().Register<PDBWriter<CSG_Topology>>("pdb");
  TrjWriterFactory().Register<XYZWriter<CSG_Topology>>("xyz");
  TrjWriterFactory().Register<GROWriter<CSG_Topology>>("gro");
#ifdef GMX_DOUBLE
  TrjWriterFactory().Register<GMXTrajectoryWriter<CSG_Topology>>("trr");
  TrjWriterFactory().Register<GMXTrajectoryWriter<CSG_Topology>>("xtc");
#endif
  TrjWriterFactory().Register<DLPOLYTrajectoryWriter<CSG_Topology>>("dlph");
  TrjWriterFactory().Register<DLPOLYTrajectoryWriter<CSG_Topology>>("dlpc");
  TrjWriterFactory().Register<LAMMPSDumpWriter<CSG_Topology>>("dump");
}
}  // namespace csg
}  // namespace votca
