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

#define BOOST_TEST_MODULE beadmoleculebase_test
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <votca/csg/basebead.h>

#include "../../include/votca/csg/basemolecule.h"

using namespace std;
using namespace votca::csg;

class TestBead : public BaseBead {
 public:
  TestBead() : BaseBead(){};
};

BOOST_AUTO_TEST_SUITE(basemolecule_test)

BOOST_AUTO_TEST_CASE(test_basemolecule_constructor) {
  BaseMolecule<TestBead> base_molecule;
}

BOOST_AUTO_TEST_CASE(test_base_molecule_beadcount) {
  BaseMolecule<TestBead> base_molecule;
  BOOST_CHECK_EQUAL(base_molecule.BeadCount(), 0);
}

BOOST_AUTO_TEST_CASE(test_base_molecule_add_and_getbead) {
  BaseMolecule<TestBead> base_molecule;
  TestBead testbead;
  testbead.setId(2);
  base_molecule.AddBead(&testbead);
  BOOST_CHECK_EQUAL(base_molecule.BeadCount(), 1);
  auto testbead2 = base_molecule.getBead(2);
  BOOST_CHECK_EQUAL(testbead2->getId(), testbead.getId());
}

BOOST_AUTO_TEST_CASE(test_base_molecule_ConnectBeads) {
  BaseMolecule<TestBead> base_molecule;
  TestBead testbead1;
  testbead1.setId(1);
  testbead1.setName("Carbon");
  TestBead testbead2;
  testbead2.setId(2);
  testbead2.setName("Carbon");
  base_molecule.AddBead(&testbead1);
  base_molecule.AddBead(&testbead2);
  base_molecule.ConnectBeads(1, 2);
}

BOOST_AUTO_TEST_CASE(test_base_molecule_isSingleStructure) {
  BaseMolecule<TestBead> base_molecule;

  TestBead testbead1;
  testbead1.setName("Carbon");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Oxygen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  base_molecule.AddBead(&testbead1);
  base_molecule.AddBead(&testbead2);
  BOOST_CHECK(!base_molecule.isSingleStructure());

  // C - C
  base_molecule.ConnectBeads(1, 2);
  BOOST_CHECK(base_molecule.isSingleStructure());

  // C - C  O
  base_molecule.AddBead(&testbead3);
  BOOST_CHECK(!base_molecule.isSingleStructure());

  // C - C - O
  base_molecule.ConnectBeads(2, 3);
  BOOST_CHECK(base_molecule.isSingleStructure());

  // C - C - O  H - H
  base_molecule.AddBead(&testbead4);
  base_molecule.AddBead(&testbead5);
  base_molecule.ConnectBeads(4, 5);
  BOOST_CHECK(!base_molecule.isSingleStructure());
}

BOOST_AUTO_TEST_CASE(test_base_molecule_isStructureEquivalent) {
  BaseMolecule<TestBead> base_molecule1;
  BaseMolecule<TestBead> base_molecule2;

  // Beads for bead structure 1
  TestBead testbead1;
  testbead1.setName("Carbon");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Oxygen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  // Beads for bead structure 2
  TestBead testbead6;
  testbead6.setName("Carbon");
  testbead6.setId(6);

  TestBead testbead7;
  testbead7.setName("Carbon");
  testbead7.setId(7);

  TestBead testbead8;
  testbead8.setName("Oxygen");
  testbead8.setId(8);

  TestBead testbead9;
  testbead9.setName("Hydrogen");
  testbead9.setId(9);

  TestBead testbead10;
  testbead10.setName("Hydrogen");
  testbead10.setId(10);

  BOOST_CHECK(base_molecule1.isStructureEquivalent(base_molecule2));
  base_molecule1.AddBead(&testbead1);
  BOOST_CHECK(!base_molecule1.isStructureEquivalent(base_molecule2));
  base_molecule2.AddBead(&testbead6);
  BOOST_CHECK(base_molecule1.isStructureEquivalent(base_molecule2));

  base_molecule1.AddBead(&testbead2);
  base_molecule2.AddBead(&testbead7);

  base_molecule1.ConnectBeads(1, 2);
  BOOST_CHECK(!base_molecule1.isStructureEquivalent(base_molecule2));
  base_molecule2.ConnectBeads(6, 7);
  BOOST_CHECK(base_molecule1.isStructureEquivalent(base_molecule2));

  base_molecule1.AddBead(&testbead3);
  base_molecule1.AddBead(&testbead4);
  base_molecule1.AddBead(&testbead5);
  base_molecule1.ConnectBeads(2, 3);
  base_molecule1.ConnectBeads(4, 5);

  base_molecule2.AddBead(&testbead10);
  base_molecule2.AddBead(&testbead8);
  base_molecule2.AddBead(&testbead9);
  base_molecule2.ConnectBeads(7, 8);
  base_molecule2.ConnectBeads(9, 10);
  BOOST_CHECK(base_molecule1.isStructureEquivalent(base_molecule2));
}

BOOST_AUTO_TEST_CASE(test_base_molecule_getNeighBeads) {
  BaseMolecule<TestBead> base_molecule1;

  // Beads for bead structure 1
  // Make a methane molecule
  //
  //     H
  //     |
  // H - C - H
  //     |
  //     H
  //
  TestBead testbead1;
  testbead1.setName("Hydrogen");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Hydrogen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  // Make a Water molecule
  //
  // H - O - H
  //

  TestBead testbead6;
  testbead6.setName("Hydrogen");
  testbead6.setId(6);

  TestBead testbead7;
  testbead7.setName("Oxygen");
  testbead7.setId(7);

  TestBead testbead8;
  testbead8.setName("Hydrogen");
  testbead8.setId(8);

  base_molecule1.AddBead(&testbead1);
  base_molecule1.AddBead(&testbead2);
  base_molecule1.AddBead(&testbead3);
  base_molecule1.AddBead(&testbead4);
  base_molecule1.AddBead(&testbead5);
  base_molecule1.AddBead(&testbead6);
  base_molecule1.AddBead(&testbead7);
  base_molecule1.AddBead(&testbead8);

  // At this point non of the beads are connected so should return a vector of
  // size 0
  auto v1 = base_molecule1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v1.size(), 0);
  auto v2 = base_molecule1.getNeighBeads(2);
  BOOST_CHECK_EQUAL(v2.size(), 0);
  auto v3 = base_molecule1.getNeighBeads(3);
  BOOST_CHECK_EQUAL(v3.size(), 0);
  auto v4 = base_molecule1.getNeighBeads(4);
  BOOST_CHECK_EQUAL(v4.size(), 0);
  auto v5 = base_molecule1.getNeighBeads(5);
  BOOST_CHECK_EQUAL(v5.size(), 0);
  auto v6 = base_molecule1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v6.size(), 0);
  auto v7 = base_molecule1.getNeighBeads(7);
  BOOST_CHECK_EQUAL(v7.size(), 0);
  auto v8 = base_molecule1.getNeighBeads(8);
  BOOST_CHECK_EQUAL(v8.size(), 0);

  // Connect beads
  base_molecule1.ConnectBeads(1, 2);
  base_molecule1.ConnectBeads(3, 2);
  base_molecule1.ConnectBeads(4, 2);
  base_molecule1.ConnectBeads(5, 2);
  base_molecule1.ConnectBeads(6, 7);
  base_molecule1.ConnectBeads(7, 8);

  v1 = base_molecule1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v1.size(), 1);
  v2 = base_molecule1.getNeighBeads(2);
  BOOST_CHECK_EQUAL(v2.size(), 4);
  v3 = base_molecule1.getNeighBeads(3);
  BOOST_CHECK_EQUAL(v3.size(), 1);
  v4 = base_molecule1.getNeighBeads(4);
  BOOST_CHECK_EQUAL(v4.size(), 1);
  v5 = base_molecule1.getNeighBeads(5);
  BOOST_CHECK_EQUAL(v5.size(), 1);
  v6 = base_molecule1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v6.size(), 1);
  v7 = base_molecule1.getNeighBeads(7);
  BOOST_CHECK_EQUAL(v7.size(), 2);
  v8 = base_molecule1.getNeighBeads(8);
  BOOST_CHECK_EQUAL(v8.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_base_molecule_catchError) {

  {
    TestBead testbead1;
    testbead1.setName("Hydrogen");
    testbead1.setId(1);

    TestBead testbead2;
    testbead2.setName("Carbon");
    testbead2.setId(2);

    TestBead testbead3;
    testbead3.setName("Hydrogen");
    testbead3.setId(3);

    TestBead testbead4;
    testbead4.setName("Hydrogen");
    testbead4.setId(4);

    TestBead testbead5;
    testbead5.setName("Hydrogen");
    testbead5.setId(5);

    TestBead testbead6;
    testbead6.setName("Hydrogen");
    testbead6.setId(5);

    BaseMolecule<TestBead> base_molecule;
    base_molecule.AddBead(&testbead1);
    base_molecule.AddBead(&testbead2);
    base_molecule.AddBead(&testbead3);
    base_molecule.AddBead(&testbead4);
    base_molecule.AddBead(&testbead5);

    BOOST_CHECK_THROW(base_molecule.AddBead(&testbead6), invalid_argument);
    BOOST_CHECK_THROW(base_molecule.ConnectBeads(0, 1), invalid_argument);
    BOOST_CHECK_THROW(base_molecule.ConnectBeads(5, 6), invalid_argument);
    BOOST_CHECK_THROW(base_molecule.ConnectBeads(1, 1), invalid_argument);
  }
}
BOOST_AUTO_TEST_SUITE_END()
