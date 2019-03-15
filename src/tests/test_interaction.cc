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

#define BOOST_TEST_MODULE interaction_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/csgtopology.h>
#include <votca/csg/interaction.h>
#include <votca/csg/molecule.h>
#include <votca/tools/constants.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(interaction_test)

BOOST_AUTO_TEST_CASE(test_interaction_setters_getters) {

  CSG_Topology top;

  string bead_type = "type1";
  byte_t symmetry = 1;
  int bead_id = 1;
  int molecule_id = 1;
  int residue_id = 1;
  string residue_type = "Protein";
  double mass = 1.1;
  double charge = 0.3;
  top.CreateBead(symmetry, bead_type, bead_id, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);
  ++bead_id;
  top.CreateBead(symmetry, bead_type, bead_id, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);
  ++bead_id;
  top.CreateBead(symmetry, bead_type, bead_id, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);
  ++bead_id;
  top.CreateBead(symmetry, bead_type, bead_id, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);

  int bond_id = 1;
  Interaction* bond1 = top.CreateInteraction(
      InteractionType::bond, "BONDS", bond_id, molecule_id, vector<int>{1, 2});
  bond1->setGroupId(1);

  string label = bond1->getLabel();
  cout << label << endl;
  bool correct_label =
      label.compare("molecule id 1:group name BONDS:group id 1:index 1") == 0;
  BOOST_CHECK(correct_label);
  BOOST_CHECK_EQUAL(bond1->getGroupId(), 1);
  BOOST_CHECK_EQUAL(bond1->getIndex(), 1);
  BOOST_CHECK_EQUAL(bond1->getMolecule(), 1);
  BOOST_CHECK_EQUAL(bond1->BeadCount(), 2);
  BOOST_CHECK_EQUAL(bond1->getBeadId(0), 1);
  BOOST_CHECK_EQUAL(bond1->getBeadId(1), 2);
  string groupName = bond1->getGroup();
  correct_label = groupName.compare("BONDS") == 0;

  BOOST_CHECK(correct_label);

  ++bond_id;
  Interaction* angle1 =
      top.CreateInteraction(InteractionType::angle, "ANGLES", bond_id,
                            molecule_id, vector<int>{1, 2, 3});
  angle1->setGroupId(1);

  label = angle1->getLabel();
  cout << label << endl;
  correct_label =
      label.compare("molecule id 1:group name ANGLES:group id 1:index 2") == 0;
  BOOST_CHECK(correct_label);
  BOOST_CHECK_EQUAL(angle1->getGroupId(), 1);
  BOOST_CHECK_EQUAL(angle1->getIndex(), 2);
  BOOST_CHECK_EQUAL(angle1->getMolecule(), 1);
  BOOST_CHECK_EQUAL(angle1->BeadCount(), 3);
  BOOST_CHECK_EQUAL(angle1->getBeadId(0), 1);
  BOOST_CHECK_EQUAL(angle1->getBeadId(1), 2);
  BOOST_CHECK_EQUAL(angle1->getBeadId(2), 3);
  groupName = angle1->getGroup();
  correct_label = groupName.compare("ANGLES") == 0;

  ++bond_id;
  Interaction* dihedral1 =
      top.CreateInteraction(InteractionType::dihedral, "DIHEDRAL", bond_id,
                            molecule_id, vector<int>{1, 2, 3, 4});
  dihedral1->setGroupId(1);

  label = dihedral1->getLabel();
  cout << label << endl;
  correct_label =
      label.compare("molecule id 1:group name DIHEDRAL:group id 1:index 3") ==
      0;
  BOOST_CHECK(correct_label);
  BOOST_CHECK_EQUAL(dihedral1->getGroupId(), 1);
  BOOST_CHECK_EQUAL(dihedral1->getIndex(), 3);
  BOOST_CHECK_EQUAL(dihedral1->getMolecule(), 1);
  BOOST_CHECK_EQUAL(dihedral1->BeadCount(), 4);
  BOOST_CHECK_EQUAL(dihedral1->getBeadId(0), 1);
  BOOST_CHECK_EQUAL(dihedral1->getBeadId(1), 2);
  BOOST_CHECK_EQUAL(dihedral1->getBeadId(2), 3);
  BOOST_CHECK_EQUAL(dihedral1->getBeadId(3), 4);
  groupName = dihedral1->getGroup();
}

BOOST_AUTO_TEST_SUITE_END()
