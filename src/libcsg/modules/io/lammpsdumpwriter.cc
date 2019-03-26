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

#include "lammpsdumpwriter.h"
#include "../../../../include/votca/csg/bead.h"
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;
void LAMMPSDumpWriter::Open(std::string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void LAMMPSDumpWriter::Close() { fclose(_out); }

void LAMMPSDumpWriter::Write(CSG_Topology *conf) {
  CSG_Topology *top = conf;
  Eigen::Matrix3d box = conf->getBox();
  fprintf(_out, "ITEM: TIMESTEP\n%i\n", top->getStep());
  fprintf(_out, "ITEM: NUMBER OF ATOMS\n%i\n", (int)top->BeadCount());
  fprintf(_out, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(_out, "0 %f\n0 %f\n0 %f\n", box(0, 0), box(1, 1), box(2, 2));

  fprintf(_out, "ITEM: ATOMS id type x y z");
  bool v = top->HasVel();
  if (v) {
    fprintf(_out, " vx vy vz");
  }
  bool f = top->HasForce();
  if (f) {
    fprintf(_out, " fx fy fz");
  }
  fprintf(_out, "\n");

  vector<int> bead_ids = conf->getBeadIds();
  // Sort the beads before outputing them
  sort(bead_ids.begin(), bead_ids.end());
  for (const int bead_id : bead_ids) {
    Bead *bead = conf->getBead(bead_id);
    int bead_type_id = conf->getBeadTypeId(bead_id);

    fprintf(_out, "%i %i", bead->getId() + 1, bead_type_id);
    fprintf(_out, " %f %f %f", bead->getPos().x(), bead->getPos().y(),
            bead->getPos().z());
    if (v) {
      fprintf(_out, " %f %f %f", bead->getVel().x(), bead->getVel().y(),
              bead->getVel().z());
    }
    if (f) {
      fprintf(_out, " %f %f %f", bead->getF().x(), bead->getF().y(),
              bead->getF().z());
    }
    fprintf(_out, "\n");
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
