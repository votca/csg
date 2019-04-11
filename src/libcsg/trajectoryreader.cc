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
#include "../../include/votca/csg/io/dlpolytrajectoryreader.h"
#include "../../include/votca/csg/io/gmxtrajectoryreader.h"
#include "../../include/votca/csg/io/pdbreader.h"
#include "../../include/votca/csg/io/xyzreader.h"
#include "../../include/votca/csg/trajectoryreader.h"
#ifdef H5MD
#include "../../include/votca/csg/h5mdtrajectoryreader.h"
#endif
#include "../../include/votca/csg/io/groreader.h"
#include "../../include/votca/csg/io/lammpsdatareader.h"
#include "../../include/votca/csg/io/lammpsdumpreader.h"
#ifdef GMX_DOUBLE
#include "../../include/votca/csg/io/gmxtrajectoryreader.h"
#endif
/*
 */
namespace votca {
namespace csg {

void TrajectoryReader::RegisterPlugins(void) {
  TrjReaderFactory().Register<PDBReader<Bead, Molecule, CSG_Topology>>("pdb");
  TrjReaderFactory()
      .Register<GMXTrajectoryReader<Bead, Molecule, CSG_Topology>>("trr");
  TrjReaderFactory()
      .Register<GMXTrajectoryReader<Bead, Molecule, CSG_Topology>>("xtc");
  TrjReaderFactory()
      .Register<DLPOLYTrajectoryReader<Bead, Molecule, CSG_Topology>>("dlph");
  TrjReaderFactory()
      .Register<DLPOLYTrajectoryReader<Bead, Molecule, CSG_Topology>>("dlpc");
  TrjReaderFactory().Register<XYZReader<Bead, Molecule, CSG_Topology>>("xyz");
#ifdef H5MD
  TrjReaderFactory()
      .Register<H5MDTrajectoryReader<Bead, Molecule, CSG_Topology>>("h5");
#endif
  TrjReaderFactory().Register<LAMMPSDumpReader<Bead, Molecule, CSG_Topology>>(
      "dump");
  TrjReaderFactory().Register<LAMMPSDataReader<Bead, Molecule, CSG_Topology>>(
      "data");
  TrjReaderFactory().Register<GROReader<Bead, Molecule, CSG_Topology>>("gro");
#ifdef GMX_DOUBLE
  TrjReaderFactory()
      .Register<GMXTrajectoryReader<Bead, Molecule, CSG_Topology>>("trr");
  TrjReaderFactory()
      .Register<GMXTrajectoryReader<Bead, Molecule, CSG_Topology>>("xtc");
#endif
}

}  // namespace csg
}  // namespace votca
