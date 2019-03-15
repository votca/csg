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

#define BOOST_TEST_MODULE cgmoleculestencil_test
#include "../../include/votca/csg/cgmoleculestencil.h"
#include "../../include/votca/csg/csgtopology.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(cgmoleculestencil_test)

BOOST_AUTO_TEST_CASE(test_cgmoleculedef_constructor) {
  CGMoleculeStencil cgmoleculestencil("Propane", "AtomicPropane");
}

BOOST_AUTO_TEST_CASE(test_load) {
  /*
    string file_name = "cg_molecule.xml";
    ofstream outfile(file_name);

    outfile << "<cg_molecule>\n";
    outfile << "  <name>ppn</name> <!-- molecule name in cg representation
    -->\n"; outfile << "  <ident>propane</ident> <!-- molecule name in atomistic
    " "topology -->\n"; outfile << "  <topology> <!-- topology of one molecule
    -->\n"; outfile << "    <cg_beads>\n"; outfile << ""; outfile << " <cg_bead>
    <!-- definition of a coarse-grained bead -->\n"; outfile << "
    <name>A1</name>\n"; outfile << "        <type>A</type>\n"; outfile << "
    <mapping>A</mapping> <!-- reference to a map -->\n"; outfile << " <!-- atoms
    belonging to this bead -->\n"; outfile << "        <beads>1:ppn:C1 1:ppn:H4
    1:ppn:H5 1:ppn:H6</beads>\n"; outfile << "      </cg_bead>\n"; outfile <<
    ""; outfile << "	     <cg_bead>\n"; outfile << " <name>B1</name>\n";
    outfile << "        <type>B</type>\n";
    outfile << "        <symmetry>1</symmetry>\n";
    outfile << "        <mapping>B</mapping>\n";
    outfile << "        <beads> 1:ppn:C2 1:ppn:H7 1:ppn:H8  </beads>\n";
    outfile << "      </cg_bead>\n";
    outfile << "";
    outfile << "      <cg_bead>\n";
    outfile << "        <name>A2</name>\n";
    outfile << "        <type>A</type>\n";
    outfile << "        <symmetry>1</symmetry>\n";
    outfile << "        <mapping>A</mapping>\n";
    outfile << "        <beads> 1:ppn:C3 1:ppn:H9 1:ppn:H10 1:ppn:H11
    </beads>\n"; outfile << "      </cg_bead>\n"; outfile << ""; outfile << "
    </cg_beads>\n"; outfile << "    <cg_bonded> <!-- bonded interactions -->\n";
    outfile << "      <bond>\n";
    outfile << "        <name>bond</name>\n";
    outfile << "        <beads>\n";
    outfile << "          A1 B1\n";
    outfile << "          B1 A2\n";
    outfile << "        </beads>\n";
    outfile << "      </bond>\n";
    outfile << "      <angle>\n";
    outfile << "        <name>angle</name>\n";
    outfile << "        <beads>\n";
    outfile << "          A1 B1 A2\n";
    outfile << "        </beads>\n";
    outfile << "      </angle>\n";
    outfile << "    </cg_bonded>\n";
    outfile << "  </topology>\n";
    outfile << "  <maps>\n";
    outfile << "    <map> <!-- mapping A -->\n";
    outfile << "      <name>A</name>\n";
    outfile << "      <weights> 12 1 1 1 </weights>\n";
    outfile << "    </map>\n";
    outfile << "    <map>\n";
    outfile << "      <name>B</name>\n";
    outfile << "      <weights> 12 1 1 </weights>\n";
    outfile << "    </map>\n";
    outfile << "  </maps>\n";
    outfile << "</cg_molecule> <!-- end of the molecule -->\n";

    outfile.close();

    CGMoleculeStencil cgmoleculestencil;
    cgmoleculestencil.Load(file_name);

    BOOST_CHECK_EQUAL(cgmoleculestencil.getCGType(), "ppn");
    BOOST_CHECK_EQUAL(cgmoleculestencil.getAtomisticType(), "propane");

    CSG_Topology top;
    Molecule* cg_molecule = cgmoleculestencil.CreateMolecule(top);

    BOOST_CHECK_EQUAL(cg_molecule->getType(), "ppn");
    // Because it is the first cg_molecule created the id should be 0
    BOOST_CHECK_EQUAL(cg_molecule->getId(), 0);

    // ids   0 - 1 - 2
    // types A - B - A
    BOOST_CHECK_EQUAL(cg_molecule->BeadCount(), 3);

    // Ids of the beads should be 0, 1, and 2 because they are the first beads
    BOOST_CHECK(cg_molecule->BeadExist(0));
    BOOST_CHECK(cg_molecule->BeadExist(1));
    BOOST_CHECK(cg_molecule->BeadExist(2));

    // The bead types of the cg molecule should be A B and A
    BOOST_CHECK_EQUAL(cg_molecule->getBead(0)->getType(), "A");
    BOOST_CHECK_EQUAL(cg_molecule->getBead(1)->getType(), "B");
    BOOST_CHECK_EQUAL(cg_molecule->getBead(2)->getType(), "A");

    // Determine the molecule the beads are attached to should be 0
    BOOST_CHECK_EQUAL(cg_molecule->getBead(0)->getMoleculeId(), 0);
    BOOST_CHECK_EQUAL(cg_molecule->getBead(1)->getMoleculeId(), 0);
    BOOST_CHECK_EQUAL(cg_molecule->getBead(2)->getMoleculeId(), 0);*/
}

BOOST_AUTO_TEST_SUITE_END()
