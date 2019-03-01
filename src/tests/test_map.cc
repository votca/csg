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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE map_test
#include "../../include/votca/csg/csgtopology.h"
#include "../../include/votca/csg/map.h"
#include "../../include/votca/csg/molecule.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/tools/property.h>
using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(map_test)

BOOST_AUTO_TEST_CASE(test_map_constructor) {

  Molecule mol_in;
  Molecule mol_out;

  Map maptest(mol_in, mol_out);
}

BOOST_AUTO_TEST_CASE(test_load) {

  string file_cg = "cg_molecule.xml";
  ofstream out_cg(file_cg);

  out_cg << "<cg_molecule>\n";
  out_cg << "  <name>ppn</name> <!-- molecule name in cg representation -->\n";
  out_cg << "  <ident>propane</ident> <!-- molecule name in atomistic topology "
            "-->\n";
  out_cg << "  <topology> <!-- topology of one molecule -->\n";
  out_cg << "    <cg_beads>\n";
  out_cg << "      <cg_bead> <!-- definition of a coarse-grained bead -->\n";
  out_cg << "        <name>A1</name>\n";
  out_cg << "        <type>A</type>\n";
  out_cg << "        <mapping>A</mapping> <!-- reference to a map -->\n";
  out_cg << "        <!-- atoms belonging to this bead -->\n";
  out_cg << "        <beads>1:ppn:C1 1:ppn:H4 1:ppn:H5 1:ppn:H6</beads>\n";
  out_cg << "      </cg_bead>\n";
  out_cg << "      <!-- more bead definitions -->\n";
  out_cg << "    </cg_beads>\n";
  out_cg << "    <cg_bonded> <!-- bonded interactions -->\n";
  out_cg << "      <bond>\n";
  out_cg << "        <name>bond</name>\n";
  out_cg << "        <beads>\n";
  out_cg << "          A1 B1\n";
  out_cg << "          B1 A2\n";
  out_cg << "        </beads>\n";
  out_cg << "      </bond>\n";
  out_cg << "      <angle>\n";
  out_cg << "        <name>angle</name>\n";
  out_cg << "        <beads>\n";
  out_cg << "          A1 B1 A2\n";
  out_cg << "        </beads>\n";
  out_cg << "      </angle>\n";
  out_cg << "    </cg_bonded>\n";
  out_cg << "  </topology>\n";
  out_cg << "  <maps>\n";
  out_cg << "    <map> <!-- mapping A -->\n";
  out_cg << "      <name>A</name>\n";
  out_cg << "      <weights> 12 1 1 1 </weights>\n";
  out_cg << "    </map>\n";
  out_cg << "    <!-- more mapping definitions -->\n";
  out_cg << "  </maps>\n";
  out_cg << "</cg_molecule> <!-- end of the molecule -->\n";
  out_cg.close();

  Property cg_options;
  load_property_from_xml(cg_options, file_cg);

  string file_topology = "topology.xml";
  ofstream out_topology(file_topology);
  out_topology << "<topology>\n";
  out_topology << "  <!-- particle group name in H5MD file -->\n";
  out_topology << "  <h5md_particle_group name=\"atoms\" />\n";
  out_topology << "  <molecules>\n";
  out_topology
      << "  <!-- define molecule, number of beads, number of mols -->\n";
  out_topology << "    <molecule name=\"BUT\" nmols=\"4000\" nbeads=\"4\">\n";
  out_topology << "      <!-- composition of molecule, bead definition -->\n";
  out_topology
      << "      <bead name=\"C1\" type=\"C\" mass=\"15.035\" q=\"0.0\" />\n";
  out_topology
      << "      <bead name=\"C2\" type=\"C\" mass=\"14.028\" q=\"0.0\" />\n";
  out_topology
      << "      <bead name=\"C3\" type=\"C\" mass=\"14.028\" q=\"0.0\" />\n";
  out_topology
      << "      <bead name=\"C4\" type=\"C\" mass=\"15.035\" q=\"0.0\" />\n";
  out_topology << "    </molecule>\n";
  out_topology << "  </molecules>\n";
  out_topology << "  <!-- bonded terms -->\n";
  out_topology << "  <bonded>\n";
  out_topology << "    <bond>\n";
  out_topology << "      <name>bond1</name>\n";
  out_topology << "      <beads>\n";
  out_topology << "        BUT:C1 BUT:C2\n";
  out_topology << "      </beads>\n";
  out_topology << "    </bond>\n";
  out_topology << "    <bond>\n";
  out_topology << "      <name>bond2</name>\n";
  out_topology << "      <beads>\n";
  out_topology << "        BUT:C2 BUT:C3\n";
  out_topology << "      </beads>\n";
  out_topology << "    </bond>\n";
  out_topology << "    <angle>\n";
  out_topology << "      <name>angle1</name>\n";
  out_topology << "      <beads>\n";
  out_topology << "        BUT:C1 BUT:C2 BUT:C3\n";
  out_topology << "        BUT:C2 BUT:C3 BUT:C4\n";
  out_topology << "      </beads>\n";
  out_topology << "    </angle>\n";
  out_topology << "    <dihedral>\n";
  out_topology << "      <name>dihedral1</name>\n";
  out_topology << "      <beads>\n";
  out_topology << "        BUT:C1 BUT:C2 BUT:C3 BUT:C4\n";
  out_topology << "      </beads>\n";
  out_topology << "    </dihedral>\n";
  out_topology << "  </bonded>\n";
  out_topology << "</topology> \n";
  out_topology.close();

  Property topology_options;
  load_property_from_xml(topology_options, file_topology);
}

BOOST_AUTO_TEST_SUITE_END()
