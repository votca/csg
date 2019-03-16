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
#include "../../include/votca/csg/cgbeadstencil.h"
#include "../../include/votca/csg/cgmoleculestencil.h"
#include "../../include/votca/csg/csgtopology.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(cgmoleculestencil_test)

BOOST_AUTO_TEST_CASE(test_cgmoleculedef_constructor) {
  CGMoleculeStencil cgmoleculestencil("Propane", "AtomicPropane");
}

BOOST_AUTO_TEST_CASE(test_load) {

  CGBeadStencil bead_stencil1;
  bead_stencil1.atomic_subbeads_ =
      vector<string>{"1:ppn:C1", "1:ppn:H4", "1:ppn:H5", "1:ppn:H6"};
  bead_stencil1.cg_name_ = "A1";
  bead_stencil1.cg_symmetry_ = 1;
  bead_stencil1.cg_bead_type_ = "A";
  bead_stencil1.mapping_ = "A";
  bead_stencil1.subbead_weights_ = vector<double>{12.0, 1.0, 1.0, 1.0};

  CGBeadStencil bead_stencil2;
  bead_stencil2.atomic_subbeads_ =
      vector<string>{"1:ppn:C2", "1:ppn:H7", "1:ppn:H8"};
  bead_stencil2.cg_name_ = "B1";
  bead_stencil2.cg_symmetry_ = 1;
  bead_stencil2.cg_bead_type_ = "B";
  bead_stencil2.mapping_ = "B";
  bead_stencil2.subbead_weights_ = vector<double>{12.0, 1.0, 1.0};

  CGBeadStencil bead_stencil3;
  bead_stencil3.atomic_subbeads_ =
      vector<string>{"1:ppn:C3", "1:ppn:H9", "1:ppn:H10", "1:ppn:H11"};
  bead_stencil3.cg_name_ = "A2";
  bead_stencil3.cg_symmetry_ = 1;
  bead_stencil3.cg_bead_type_ = "A";
  bead_stencil3.mapping_ = "A";
  bead_stencil3.subbead_weights_ = vector<double>{12.0, 1.0, 1.0, 1.0};

  vector<CGBeadStencil> bead_stencils;
  bead_stencils.push_back(bead_stencil1);
  bead_stencils.push_back(bead_stencil2);
  bead_stencils.push_back(bead_stencil3);

  CGMoleculeStencil molecule_stencil("propane", "AtomicPropane");
  molecule_stencil.AddBeadStencil(bead_stencils);

  cout << "Testing getAtomicBeadNames A1" << endl;
  {
    vector<string> atomic_bead_names =
        molecule_stencil.getAtomicBeadNames("A1");

    // Check that the correct atomic beads are registered with cg-bead name "A1"
    map<string, bool> bead_found;
    bead_found["1:ppn:C1"] = false;
    bead_found["1:ppn:H4"] = false;
    bead_found["1:ppn:H5"] = false;
    bead_found["1:ppn:H6"] = false;
    for (string& atomic_bead_name : atomic_bead_names) {
      if (atomic_bead_name == "1:ppn:C1")
        bead_found["1:ppn:C1"] = true;
      else if (atomic_bead_name == "1:ppn:H4")
        bead_found["1:ppn:H4"] = true;
      else if (atomic_bead_name == "1:ppn:H5")
        bead_found["1:ppn:H5"] = true;
      else if (atomic_bead_name == "1:ppn:H6")
        bead_found["1:ppn:H6"] = true;
    }

    BOOST_CHECK(bead_found["1:ppn:C1"]);
    BOOST_CHECK(bead_found["1:ppn:H4"]);
    BOOST_CHECK(bead_found["1:ppn:H5"]);
    BOOST_CHECK(bead_found["1:ppn:H6"]);
    BOOST_CHECK_EQUAL(atomic_bead_names.size(), 4);
  }

  cout << "Testing getAtomicBeadNames B1" << endl;
  {
    vector<string> atomic_bead_names =
        molecule_stencil.getAtomicBeadNames("B1");

    // Check that the correct atomic beads are registered with cg-bead name "B1"
    map<string, bool> bead_found;
    bead_found["1:ppn:C2"] = false;
    bead_found["1:ppn:H7"] = false;
    bead_found["1:ppn:H8"] = false;
    for (string& atomic_bead_name : atomic_bead_names) {
      cout << atomic_bead_name << endl;
      if (atomic_bead_name == "1:ppn:C2")
        bead_found["1:ppn:C2"] = true;
      else if (atomic_bead_name == "1:ppn:H7")
        bead_found["1:ppn:H7"] = true;
      else if (atomic_bead_name == "1:ppn:H8")
        bead_found["1:ppn:H8"] = true;
    }

    BOOST_CHECK(bead_found["1:ppn:C2"]);
    BOOST_CHECK(bead_found["1:ppn:H7"]);
    BOOST_CHECK(bead_found["1:ppn:H8"]);
    BOOST_CHECK_EQUAL(atomic_bead_names.size(), 3);
  }

  cout << "Testing getAtomicBeadNames A2" << endl;
  {
    vector<string> atomic_bead_names =
        molecule_stencil.getAtomicBeadNames("A2");

    // Check that the correct atomic beads are registered with cg-bead name "A2"
    map<string, bool> bead_found;
    bead_found["1:ppn:C3"] = false;
    bead_found["1:ppn:H9"] = false;
    bead_found["1:ppn:H10"] = false;
    bead_found["1:ppn:H11"] = false;
    for (string& atomic_bead_name : atomic_bead_names) {
      if (atomic_bead_name == "1:ppn:C3")
        bead_found["1:ppn:C3"] = true;
      else if (atomic_bead_name == "1:ppn:H9")
        bead_found["1:ppn:H9"] = true;
      else if (atomic_bead_name == "1:ppn:H10")
        bead_found["1:ppn:H10"] = true;
      else if (atomic_bead_name == "1:ppn:H11")
        bead_found["1:ppn:H11"] = true;
    }

    BOOST_CHECK(bead_found["1:ppn:C3"]);
    BOOST_CHECK(bead_found["1:ppn:H9"]);
    BOOST_CHECK(bead_found["1:ppn:H10"]);
    BOOST_CHECK(bead_found["1:ppn:H11"]);
    BOOST_CHECK_EQUAL(atomic_bead_names.size(), 4);
  }

  // Returns the names of all the cg beads should be A1, B1 and A2
  cout << "Testing getCGBeadNames" << endl;
  {
    vector<string> cg_bead_names = molecule_stencil.getCGBeadNames();

    map<string, bool> bead_found;
    bead_found["A1"] = false;
    bead_found["B1"] = false;
    bead_found["A2"] = false;
    for (string& cg_bead_name : cg_bead_names) {
      if (cg_bead_name == "A1")
        bead_found["A1"] = true;
      else if (cg_bead_name == "B1")
        bead_found["B1"] = true;
      else if (cg_bead_name == "A2")
        bead_found["A2"] = true;
    }

    BOOST_CHECK(bead_found["A1"]);
    BOOST_CHECK(bead_found["B1"]);
    BOOST_CHECK(bead_found["A2"]);
    BOOST_CHECK_EQUAL(cg_bead_names.size(), 3);
  }

  /*
   * Given the the bead ids in an atomic molecule map them to the bead names.
   * This function makes some assumptions, it assumes that all the ids passed
   * in belong to the same molecule and it assumes that when the ids are sorted
   * from highest to lowest that they line up with the atomic bead ids passed
   * in.Note that bead names will come out in the order that the BeadStencils
   * were passed into the CGMoleculeStencil
   */
  cout << "Testing MapAtomicBeadIdsToAtomicBeadNames 1" << endl;
  {
    vector<int> atomic_bead_ids = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    unordered_map<int, string> atomic_name_and_id =
        molecule_stencil.MapAtomicBeadIdsToAtomicBeadNames(atomic_bead_ids);

    BOOST_CHECK_EQUAL(atomic_name_and_id[1], "1:ppn:C1");
    BOOST_CHECK_EQUAL(atomic_name_and_id[2], "1:ppn:H4");
    BOOST_CHECK_EQUAL(atomic_name_and_id[3], "1:ppn:H5");
    BOOST_CHECK_EQUAL(atomic_name_and_id[4], "1:ppn:H6");
    BOOST_CHECK_EQUAL(atomic_name_and_id[5], "1:ppn:C2");
    BOOST_CHECK_EQUAL(atomic_name_and_id[6], "1:ppn:H7");
    BOOST_CHECK_EQUAL(atomic_name_and_id[7], "1:ppn:H8");
    BOOST_CHECK_EQUAL(atomic_name_and_id[8], "1:ppn:C3");
    BOOST_CHECK_EQUAL(atomic_name_and_id[9], "1:ppn:H9");
    BOOST_CHECK_EQUAL(atomic_name_and_id[10], "1:ppn:H10");
    BOOST_CHECK_EQUAL(atomic_name_and_id[11], "1:ppn:H11");
    BOOST_CHECK_EQUAL(atomic_name_and_id.size(), 11);
  }

  cout << "Testing MapAtomicBeadIdsToAtomicBeadNames 2" << endl;
  {
    // The atomic beads should be sorted in order
    vector<int> atomic_bead_ids = {12, 16, 13, 15, 17, 18, 19, 20, 14, 22, 21};
    unordered_map<int, string> atomic_name_and_id =
        molecule_stencil.MapAtomicBeadIdsToAtomicBeadNames(atomic_bead_ids);

    BOOST_CHECK_EQUAL(atomic_name_and_id[12], "1:ppn:C1");
    BOOST_CHECK_EQUAL(atomic_name_and_id[13], "1:ppn:H4");
    BOOST_CHECK_EQUAL(atomic_name_and_id[14], "1:ppn:H5");
    BOOST_CHECK_EQUAL(atomic_name_and_id[15], "1:ppn:H6");
    BOOST_CHECK_EQUAL(atomic_name_and_id[16], "1:ppn:C2");
    BOOST_CHECK_EQUAL(atomic_name_and_id[17], "1:ppn:H7");
    BOOST_CHECK_EQUAL(atomic_name_and_id[18], "1:ppn:H8");
    BOOST_CHECK_EQUAL(atomic_name_and_id[19], "1:ppn:C3");
    BOOST_CHECK_EQUAL(atomic_name_and_id[20], "1:ppn:H9");
    BOOST_CHECK_EQUAL(atomic_name_and_id[21], "1:ppn:H10");
    BOOST_CHECK_EQUAL(atomic_name_and_id[22], "1:ppn:H11");
    BOOST_CHECK_EQUAL(atomic_name_and_id.size(), 11);
  }

  cout << "Testing MapCGBeadIdsToCGBeadNames 1" << endl;
  {
    vector<int> bead_ids = {1, 3, 2};
    unordered_map<int, string> bead_id_and_name =
        molecule_stencil.MapCGBeadIdsToCGBeadNames(bead_ids);
    BOOST_CHECK_EQUAL(bead_id_and_name[1], "A1");
    BOOST_CHECK_EQUAL(bead_id_and_name[2], "B1");
    BOOST_CHECK_EQUAL(bead_id_and_name[3], "A2");
  }

  cout << "Testing MapCGBeadIdsToCGBeadNames 2" << endl;
  {
    vector<int> bead_ids = {12, 10, 11};
    unordered_map<int, string> bead_id_and_name =
        molecule_stencil.MapCGBeadIdsToCGBeadNames(bead_ids);
    BOOST_CHECK_EQUAL(bead_id_and_name[10], "A1");
    BOOST_CHECK_EQUAL(bead_id_and_name[11], "B1");
    BOOST_CHECK_EQUAL(bead_id_and_name[12], "A2");
  }
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
