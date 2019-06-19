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

#define BOOST_TEST_MODULE topologytypecontainer_test
#include "../../include/votca/csg/topologytypecontainer.h"
#include <boost/test/unit_test.hpp>
#include <string>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(topologytypecontainer_test)

BOOST_AUTO_TEST_CASE(test_topologytypecontainer_constructor) {
  TopologyTypeContainer type_cont;
}

BOOST_AUTO_TEST_CASE(test_molecule_methods) {

  TopologyTypeContainer type_cont;
  string molecule_type = "propane";
  BOOST_CHECK(!type_cont.MoleculeTypeExist(molecule_type));
  BOOST_CHECK_EQUAL(type_cont.MoleculeTypeCount(), 0);

  type_cont.AddMoleculeType(molecule_type);
  BOOST_CHECK(type_cont.MoleculeTypeExist(molecule_type));
  BOOST_CHECK_EQUAL(type_cont.MoleculeTypeCount(), 1);

  // Adding the same type should do nothing
  type_cont.AddMoleculeType(molecule_type);
  BOOST_CHECK_EQUAL(type_cont.MoleculeTypeCount(), 1);

  string molecule_type2 = "rubrene";
  type_cont.AddMoleculeType(molecule_type2);
  BOOST_CHECK(type_cont.MoleculeTypeExist(molecule_type2));
  BOOST_CHECK_EQUAL(type_cont.MoleculeTypeCount(), 2);

  vector<string> molecule_types = type_cont.getMoleculeTypes();

  sort(molecule_types.begin(), molecule_types.end());
  BOOST_CHECK_EQUAL(molecule_types.size(), 2);
  BOOST_CHECK_EQUAL(molecule_types.at(0), molecule_type);
  BOOST_CHECK_EQUAL(molecule_types.at(1), molecule_type2);
}

BOOST_AUTO_TEST_CASE(test_residue_methods) {

  TopologyTypeContainer type_cont;
  string residue_type = "arginine";
  BOOST_CHECK(!type_cont.ResidueTypeExist(residue_type));
  BOOST_CHECK_EQUAL(type_cont.ResidueTypeCount(), 0);

  type_cont.AddResidueType(residue_type);
  BOOST_CHECK(type_cont.ResidueTypeExist(residue_type));
  BOOST_CHECK_EQUAL(type_cont.ResidueTypeCount(), 1);

  // Adding the same type should do nothing
  type_cont.AddResidueType(residue_type);
  BOOST_CHECK_EQUAL(type_cont.ResidueTypeCount(), 1);

  string residue_type2 = "histidine";
  type_cont.AddResidueType(residue_type2);
  BOOST_CHECK(type_cont.ResidueTypeExist(residue_type2));
  BOOST_CHECK_EQUAL(type_cont.ResidueTypeCount(), 2);

  vector<string> residue_types = type_cont.getResidueTypes();

  sort(residue_types.begin(), residue_types.end());
  BOOST_CHECK_EQUAL(residue_types.size(), 2);
  BOOST_CHECK_EQUAL(residue_types.at(0), residue_type);
  BOOST_CHECK_EQUAL(residue_types.at(1), residue_type2);
}

BOOST_AUTO_TEST_CASE(test_bead_methods) {

  TopologyTypeContainer type_cont;
  string bead_type = "CH3";
  BOOST_CHECK(!type_cont.BeadTypeExist(bead_type));
  BOOST_CHECK_EQUAL(type_cont.BeadTypeCount(), 0);

  type_cont.AddBeadType(bead_type);
  BOOST_CHECK(type_cont.BeadTypeExist(bead_type));
  BOOST_CHECK_EQUAL(type_cont.BeadTypeCount(), 1);

  // Adding the same type should do nothing
  type_cont.AddBeadType(bead_type);
  BOOST_CHECK_EQUAL(type_cont.BeadTypeCount(), 1);

  string bead_type2 = "CH2";
  type_cont.AddBeadType(bead_type2);
  BOOST_CHECK(type_cont.BeadTypeExist(bead_type2));
  BOOST_CHECK_EQUAL(type_cont.BeadTypeCount(), 2);

  vector<string> bead_types = type_cont.getBeadTypes();

  sort(bead_types.begin(), bead_types.end());
  BOOST_CHECK_EQUAL(bead_types.size(), 2);
  BOOST_CHECK_EQUAL(bead_types.at(0), bead_type2);
  BOOST_CHECK_EQUAL(bead_types.at(1), bead_type);
}

BOOST_AUTO_TEST_SUITE_END()
