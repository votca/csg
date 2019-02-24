/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE csg_topology_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <math.h>
#include <string>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>

#include "../../include/votca/csg/csgtopology.h"
#include "../../include/votca/csg/interaction.h"

using namespace std;
using namespace votca::tools;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(csg_topology_test)

BOOST_AUTO_TEST_CASE(constructors_test) { CSG_Topology top; }

BOOST_AUTO_TEST_CASE(box_test) {

  // Box takes a vector
  double x1 = 2.0;
  double y1 = 0.0;
  double z1 = 0.0;

  double x2 = 0.0;
  double y2 = 2.0;
  double z2 = 0.0;

  double x3 = 0.0;
  double y3 = 0.0;
  double z3 = 2.0;

  vec v1(x1, y1, z1);
  vec v2(x2, y2, z2);
  vec v3(x3, y3, z3);

  matrix box(v1, v2, v3);

  CSG_Topology top;
  top.setBox(box);

  auto vol = top.BoxVolume();

  BOOST_CHECK_CLOSE(vol, 8, 1e-5);
  auto box2 = top.getBox();

  BOOST_CHECK_EQUAL(box2.isClose(box, 1e-5), true);
}

BOOST_AUTO_TEST_CASE(simple_test) {

  CSG_Topology top;
  top.setStep(1);
  BOOST_CHECK_EQUAL(top.getStep(), 1);
  top.setTime(1.21);

  BOOST_CHECK_CLOSE(top.getTime(), 1.21, 1e-5);
}

/**
 * Test is creating a bead with the topology object and then ensuring that the
 * bead has the correct properties.
 **/
BOOST_AUTO_TEST_CASE(create_bead) {
  CSG_Topology top;
  // 1 - for spherical bead
  byte_t symmetry = 1;

  string bead_type = "type1";
  int bead_id = 1;
  int molecule_id = 1;
  int residue_id = 1;
  string residue_type = "Protein";
  double mass = 1.1;
  double charge = 0.3;

  auto bead_ptr = top.CreateBead(
      symmetry, bead_type, bead_id, molecule_id, residue_id, residue_type,
      basebead_constants::unassigned_element, mass, charge);

  BOOST_CHECK_CLOSE(bead_ptr->getQ(), 0.3, 1e-5);
  BOOST_CHECK_CLOSE(bead_ptr->getMass(), 1.1, 1e-5);
  BOOST_CHECK_EQUAL(bead_ptr->getResidueId(), residue_id);
  BOOST_CHECK_EQUAL(bead_ptr->getMoleculeId(), molecule_id);
  BOOST_CHECK_EQUAL(bead_ptr->getSymmetry(), symmetry);
  BOOST_CHECK(bead_ptr->getType() == bead_type);
  BOOST_CHECK_EQUAL(top.BeadCount(), 1);

  top.Cleanup();
}

/**
 * This test ensures that the interactions are stored correctly. Three beads
 * are created and two interactions are created connecting the three beads. The
 * interactions are then checked to be sure they hold the right information.
 **/
BOOST_AUTO_TEST_CASE(add_bonded_interation_test) {
  CSG_Topology top;
  // 1 - for spherical bead
  byte_t symmetry = 1;

  string bead_type = "type1";

  int bead_id = 0;
  int molecule_id = 1;
  int residue_id = 1;
  string residue_type = "Protein";
  double mass = 1.1;
  double charge = 0.3;

  // Create 3 beads
  auto bead_ptr = top.CreateBead(
      symmetry, bead_type, bead_id, molecule_id, residue_id, residue_type,
      basebead_constants::unassigned_element, mass, charge);
  bead_ptr->setId(0);

  bead_id = 1;
  auto bead_ptr2 = top.CreateBead(
      symmetry, bead_type, bead_id, molecule_id, residue_id, residue_type,
      basebead_constants::unassigned_element, mass, charge);
  bead_ptr2->setId(1);

  bead_id = 2;
  string bead_type3 = "bead_test3";
  auto bead_ptr3 = top.CreateBead(
      symmetry, bead_type3, bead_id, molecule_id, residue_id, residue_type,
      basebead_constants::unassigned_element, mass, charge);
  bead_ptr3->setId(2);

  BOOST_CHECK_EQUAL(top.BeadCount(), 3);

  // Create two bonded interactions
  int bond_id = 0;
  string interaction_group = "BONDS";
  Interaction* bond1 = top.CreateInteraction(
      Interaction::interaction_type::bond, interaction_group, bond_id,
      molecule_id, vector<int>{0, 1});

  ++bond_id;
  // bond1->setGroup(interaction_group);
  // auto bond2 = new IBond(1, 2);
  Interaction* bond2 = top.CreateInteraction(
      Interaction::interaction_type::bond, interaction_group, bond_id,
      molecule_id, vector<int>{1, 2});
  // bond2->setGroup(interaction_group);

  // top.AddBondedInteraction(bond1);
  // top.AddBondedInteraction(bond2);

  const vector<unique_ptr<Interaction>>& interaction_container =
      top.BondedInteractions();
  BOOST_CHECK_EQUAL(interaction_container.size(), 2);

  cout << "interaction name " << interaction_container.at(0)->getName() << endl;
  cout << "interaction name " << interaction_container.at(1)->getName() << endl;
  BOOST_CHECK_EQUAL(interaction_container.at(0)->getBeadId(0), 0);
  BOOST_CHECK_EQUAL(interaction_container.at(0)->getBeadId(1), 1);
  BOOST_CHECK_EQUAL(interaction_container.at(1)->getBeadId(0), 1);
  BOOST_CHECK_EQUAL(interaction_container.at(1)->getBeadId(1), 2);

  top.Cleanup();
}

BOOST_AUTO_TEST_SUITE_END()
