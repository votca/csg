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

#define BOOST_TEST_MODULE cgmoleculedef_test
#include "../../include/votca/csg/cgmoleculedef.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(cgmoleculedef_test)

BOOST_AUTO_TEST_CASE(test_cgmoleculedef_constructor) {
  CGMoleculeDef cgmoleculedef;
}

BOOST_AUTO_TEST_CASE(test_load) {

  string file_name = "cg_molecule.xml";
  ofstream outfile(file_name);

  outfile << "<cg_molecule>" << endl;
  outfile << "  <name>ppn</name> <!-- molecule name in cg representation -->"
          << endl;
  outfile
      << "  <ident>propane</ident> <!-- molecule name in atomistic topology -->"
      << endl;
  outfile << "  <topology> <!-- topology of one molecule -->" << endl;
  outfile << "    <cg_beads>" << endl;
  outfile << "      <cg_bead> <!-- definition of a coarse-grained bead -->"
          << endl;
  outfile << "        <name>A1</name>" << endl;
  outfile << "        <type>A</type>" << endl;
  outfile << "        <mapping>A</mapping> <!-- reference to a map -->" << endl;
  outfile << "        <!-- atoms belonging to this bead -->" << endl;
  outfile << "        <beads>1:ppn:C1 1:ppn:H4 1:ppn:H5 1:ppn:H6</beads>"
          << endl;
  outfile << "      </cg_bead>" << endl;
  outfile << "      <!-- more bead definitions -->" << endl;
  outfile << "    </cg_beads>" << endl;
  outfile << "    <cg_bonded> <!-- bonded interactions -->" << endl;
  outfile << "      <bond>" << endl;
  outfile << "        <name>bond</name>" << endl;
  outfile << "        <beads>" << endl;
  outfile << "          A1 B1" << endl;
  outfile << "          B1 A2" << endl;
  outfile << "        </beads>" << endl;
  outfile << "      </bond>" << endl;
  outfile << "      <angle>" << endl;
  outfile << "        <name>angle</name>" << endl;
  outfile << "        <beads>" << endl;
  outfile << "          A1 B1 A2" << endl;
  outfile << "        </beads>" << endl;
  outfile << "      </angle>" << endl;
  outfile << "    </cg_bonded>" << endl;
  outfile << "  </topology>" << endl;
  outfile << "  <maps>" << endl;
  outfile << "    <map> <!-- mapping A -->" << endl;
  outfile << "      <name>A</name>" << endl;
  outfile << "      <weights> 12 1 1 1 </weights>" << endl;
  outfile << "    </map>" << endl;
  outfile << "    <!-- more mapping definitions -->" << endl;
  outfile << "  </maps>" << endl;
  outfile << "</cg_molecule> <!-- end of the molecule -->" << endl;

  outfile.close();

  CGMoleculeDef cgmoleculedef;
  cgmoleculedef.Load(file_name);

  BOOST_CHECK_EQUAL(cgmoleculedef.getCGType(), "ppn");
  BOOST_CHECK_EQUAL(cgmoleculedef.getAtomisticType(), "propane");
}

BOOST_AUTO_TEST_SUITE_END()
