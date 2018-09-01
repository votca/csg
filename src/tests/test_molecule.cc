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
  Molecule * mol = top.CreateMolecule(molecule_name);
}

BOOST_AUTO_TEST_CASE(test_molecule_base_methods){

  Topology top;
  string molecule_name = "Methane";
  Molecule * mol = top.CreateMolecule(molecule_name);
  BOOST_CHECK_EQUAL(mol->getId(),0);
  bool name_same = !(molecule_name.compare(mol->getName()));
  BOOST_CHECK(name_same);
  
}

BOOST_AUTO_TEST_CASE(test_molecule_bead_methods){

  Topology top;
  string molecule_name = "Methane";
  Molecule * mol = top.CreateMolecule(molecule_name);
 
  BOOST_CHECK_EQUAL(mol->BeadCount(),0);

  // Check that errors are thrown if bead does not exist in the molecule
  BOOST_CHECK_THROW(mol->getBead(0),invalid_argument);
  auto beads_with_same_name = mol->getIdsOfBeadsWithName("C1");
  BOOST_CHECK_EQUAL(beads_with_same_name.size(),0);
  BOOST_CHECK_THROW(mol->getBeadName(0),invalid_argument);
   
  // Create bead
  string bead_type_name = "1";
  BeadType * b_type = top.GetOrCreateBeadType(bead_type_name);

  int symmetry = 1;
  string name = "C1";
  int resnr = 0;
  double mass = 1.21;
  double charge = -0.87;
  Bead * b = top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

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
  Molecule * mol = top.CreateMolecule(molecule_name);

  string bead_type_name = "1";
  BeadType * b_type = top.GetOrCreateBeadType(bead_type_name);

  int symmetry = 1;
  string name = "C1";
  int resnr = 0;
  double mass = 1.21;
  double charge = -0.87;
  Bead * b = top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

  mol->AddBead(b);
 
  string name2 = "H1";
  Bead * b2 = top.CreateBead(symmetry,name2,b_type,resnr,mass,charge);
  mol->AddBead(b2);

  BOOST_CHECK_EQUAL(mol->BeadCount(),2);
  
  IBond bond_interaction(0,1);  
  mol->AddInteraction(&bond_interaction);

  vector<Interaction *> bond_interactions = mol->Interactions();
  BOOST_CHECK_EQUAL(bond_interactions.size(),1);

}

  
BOOST_AUTO_TEST_SUITE_END()
