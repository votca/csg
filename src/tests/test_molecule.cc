/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE beadstructure_test

#include <boost/test/unit_test.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <votca/csg/topology.h>
#include <votca/csg/molecule.h>
#include <votca/csg/bead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/interaction.h>

using namespace std;
using namespace votca::csg;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(molecule_test)

BOOST_AUTO_TEST_CASE(test_molecule_constructor) {
  Topology top;
  string molecule_name = "Methane";
  auto mol = top.CreateMolecule(molecule_name);
}

BOOST_AUTO_TEST_CASE(test_molecule_base_methods){

  Topology top;
  string molecule_name = "Methane";
  auto mol = top.CreateMolecule(molecule_name);
  BOOST_CHECK_EQUAL(mol->getId(),0);
  bool name_same = !(molecule_name.compare(mol->getName()));
  BOOST_CHECK(name_same);
  
}

BOOST_AUTO_TEST_CASE(test_molecule_bead_methods){

  Topology top;
  string molecule_name = "Methane";
  auto mol = top.CreateMolecule(molecule_name);
 
  BOOST_CHECK_EQUAL(mol->BeadCount(),0);

  // Check that errors are thrown if bead does not exist in the molecule
  BOOST_CHECK_THROW(mol->getBead(0),invalid_argument);
  auto beads_with_same_name = mol->getIdsOfBeadsWithName("C1");
  BOOST_CHECK_EQUAL(beads_with_same_name.size(),0);
  BOOST_CHECK_THROW(mol->getBeadName(0),invalid_argument);
   
  // Create bead
  string bead_type_name = "1";
  auto b_type = top.GetOrCreateBeadType(bead_type_name);

  int symmetry = 1;
  string name = "C1";
  int resnr = 0;
  double mass = 1.21;
  double charge = -0.87;
  auto b = top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

  mol->AddBead(b);
  
  auto b2 = mol->getBead(0);
  BOOST_CHECK_EQUAL(b2->getId(),b->getId());
  beads_with_same_name =  mol->getIdsOfBeadsWithName("C1");
  BOOST_CHECK_EQUAL(beads_with_same_name.size(),1);
  BOOST_CHECK_EQUAL(beads_with_same_name.at(0),0);
  BOOST_CHECK_EQUAL(mol->BeadCount(),1);
  
  bool string_eql = !(mol->getBeadName(0).compare("C1"));
  BOOST_CHECK(string_eql);
  
}

BOOST_AUTO_TEST_CASE(test_molecule_interaction){

  Topology top;
  string molecule_name = "Methane";
  auto mol = top.CreateMolecule(molecule_name);

  string bead_type_name2 = "1";
  auto b_type = top.GetOrCreateBeadType(bead_type_name2);

  int symmetry = 1;
  string name = "C1";
  int resnr = 0;
  double mass = 1.21;
  double charge = -0.87;
  auto b = top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

  mol->AddBead(b);
 
  string name2 = "H1";
  string bead_type_name = "2";
  mass = 0.1;
  charge = 0.4;
  auto b_type2 = top.GetOrCreateBeadType(bead_type_name2);
  auto b2 = top.CreateBead(symmetry,name2,b_type2,resnr,mass,charge);
  mol->AddBead(b2);

  string name3 = "H2";
  b2 = top.CreateBead(symmetry,name3,b_type2,resnr,mass,charge);
  mol->AddBead(b2);
  string name4 = "H3";
  b2 = top.CreateBead(symmetry,name4,b_type2,resnr,mass,charge);
  mol->AddBead(b2);
  string name5 = "H4";
  b2 = top.CreateBead(symmetry,name5,b_type2,resnr,mass,charge);
  mol->AddBead(b2);
  BOOST_CHECK_EQUAL(mol->BeadCount(),5);
  
  auto ibond = shared_ptr<IBond>(new IBond(0,1));  
  mol->AddInteraction(ibond);

  auto ibond2 = shared_ptr<IBond>(new IBond(0,2));  
  mol->AddInteraction(ibond2);

  auto ibond3 = shared_ptr<IBond>(new IBond(0,3));  
  mol->AddInteraction(ibond3);

  auto ibond4 = shared_ptr<IBond>(new IBond(0,4));  
  mol->AddInteraction(ibond4);

  auto bond_interactions = mol->Interactions();
  BOOST_CHECK_EQUAL(bond_interactions.size(),4);

}

BOOST_AUTO_TEST_CASE(test_comparison_of_molecules){

  Topology top;
  string molecule_name = "Methane";
  auto mol = top.CreateMolecule(molecule_name);

  string molecule_name2 = "Methane2";
  auto mol2 = top.CreateMolecule(molecule_name2);

  string bead_type_name = "1";
  auto b_type = top.GetOrCreateBeadType(bead_type_name);

  int symmetry = 1;
  double mass = 1.21;
  double charge = -0.87;

  int resnr = 0;
  string name_C1 = "C";
  auto b1_C = top.CreateBead(symmetry,name_C1,b_type,resnr,mass,charge);
  b1_C->setId(0);
  
  mol->AddBead(b1_C);

  // Adding the hydrogens 
  string bead_type_name2 = "2";
  mass = 0.1;
  charge = 0.4;
  auto b_type2 = top.GetOrCreateBeadType(bead_type_name2);

  string name2 = "H";
  auto b1_H = top.CreateBead(symmetry,name2,b_type2,resnr,mass,charge);
  b1_H->setId(1);
  mol->AddBead(b1_H);
  string name3 = "H";
  auto b2_H = top.CreateBead(symmetry,name3,b_type2,resnr,mass,charge);
  b2_H->setId(2);
  mol->AddBead(b2_H);
  string name4 = "H";
  auto b3_H = top.CreateBead(symmetry,name4,b_type2,resnr,mass,charge);
  b3_H->setId(3);
  mol->AddBead(b3_H);
  string name5 = "H";
  auto b4_H = top.CreateBead(symmetry,name5,b_type2,resnr,mass,charge);
  b4_H->setId(4);
  mol->AddBead(b4_H);

  mol->ConnectBeads(0,1);
  mol->ConnectBeads(0,2);
  mol->ConnectBeads(0,3);
  mol->ConnectBeads(0,4);

  //
  //      H1
  //      |
  // H2 - C0 - H4
  //      |
  //      H5
  //


  BOOST_CHECK_EQUAL(mol->BeadCount(),5);

  // Create the second molecule 
  mass = 1.21;
  charge = -0.87;
  int resnr2 = 1;
  string name_C2 = "C";
  auto b2_C = top.CreateBead(symmetry,name_C2,b_type,resnr2,mass,charge);
  b2_C->setId(5);
  mol2->AddBead(b2_C);

  mass = 0.1;
  charge = 0.4;
  string name6 = "H"; 
  auto b5_H = top.CreateBead(symmetry,name6,b_type2,resnr2,mass,charge);
  b5_H->setId(6);
  mol2->AddBead(b5_H);
  string name7 = "H";
  auto b6_H = top.CreateBead(symmetry,name7,b_type2,resnr2,mass,charge);
  b5_H->setId(7);
  mol2->AddBead(b6_H);
  string name8 = "H";
  auto b7_H = top.CreateBead(symmetry,name8,b_type2,resnr2,mass,charge);
  b5_H->setId(8);
  mol2->AddBead(b7_H);

  mol2->ConnectBeads(5,6);
  mol2->ConnectBeads(5,7);
  mol2->ConnectBeads(5,8);

  //
  //      H1              H7
  //      |               |
  // H2 - C0 - H4    H6 - C5 - H8
  //      |
  //      H3
  //
  BOOST_CHECK(!mol->isStructureEquivalent(*mol2));

  string name9 = "H";
  auto b8_H = top.CreateBead(symmetry,name9,b_type2,resnr2,mass,charge);
  b8_H->setId(9);
  mol2->AddBead(b8_H);
  mol2->ConnectBeads(5,9);

  //
  //      H1              H7
  //      |               |
  // H2 - C0 - H4    H6 - C5 - H8
  //      |               |
  //      H3              H9
  //
  BOOST_CHECK(mol->isStructureEquivalent(*mol2));
  
}

  
BOOST_AUTO_TEST_SUITE_END()
