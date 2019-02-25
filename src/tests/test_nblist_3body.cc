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

#define BOOST_TEST_MODULE nblist_3body_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <votca/csg/bead.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/csgtopology.h>
#include <votca/csg/nblist_3body.h>
#include <votca/tools/constants.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(nblist_3body_test)

BOOST_AUTO_TEST_CASE(test_nblist_3body_constructor) { NBList_3Body nb; }

BOOST_AUTO_TEST_CASE(test_nblist_3body_generate_list) {
  NBList_3Body *nb;
  nb = new NBList_3Body();

  nb->setCutoff(2.0);

  CSG_Topology top;

  matrix m;
  m.ZeroMatrix();
  m[0][0] = 5.0;
  m[1][1] = 5.0;
  m[2][2] = 5.0;

  top.setBox(m);

  vec pos;

  Molecule *mol;
  int molecule_id = 1;

  mol = top.CreateMolecule(molecule_id,
                           topology_constants::unassigned_molecule_type);

  int symmetry = 1;
  string bead_type = "CG";
  int bead_id0 = 0;
  string residue_type = "Protein";
  int residue_id = 0;
  double mass = 1.0;
  double charge = -1.0;
  Bead *b;
  b = top.CreateBead(symmetry, bead_type, bead_id0, molecule_id, residue_id,
                     residue_type, topology_constants::unassigned_element, mass,
                     charge);
  pos[0] = 0.0;
  pos[1] = 0.0;
  pos[2] = 0.0;
  b->setPos(pos);
  mol->AddBead(b);
  b->setMoleculeId(mol->getId());

  int bead_id1 = 1;
  symmetry = 1;
  bead_type = "CG";
  residue_id = 0;
  mass = 2.0;
  charge = -2.0;
  b = top.CreateBead(symmetry, bead_type, bead_id1, molecule_id, residue_id,
                     residue_type, topology_constants::unassigned_element, mass,
                     charge);
  mol->AddBead(b);
  b->setMoleculeId(mol->getId());
  pos[0] = 1.0;
  pos[1] = 0.0;
  pos[2] = 0.0;
  b->setPos(pos);

  int bead_id2 = 2;
  symmetry = 1;
  bead_type = "CG";
  residue_id = 0;
  mass = 3.0;
  charge = -3.0;
  b = top.CreateBead(symmetry, bead_type, bead_id2, molecule_id, residue_id,
                     residue_type, topology_constants::unassigned_element, mass,
                     charge);
  mol->AddBead(b);
  b->setMoleculeId(mol->getId());
  pos[0] = 1.0;
  pos[1] = 1.0;
  pos[2] = 0.0;
  b->setPos(pos);

  BeadList beads;
  beads.Generate(top, "CG");

  nb->Generate(beads, true);

  BOOST_CHECK_EQUAL(nb->size(), 3);

  for (NBList_3Body::iterator triple_iter = nb->begin();
       triple_iter != nb->end(); ++triple_iter) {

    // Basically check that all 3 ids { 0, 1, 2 } are present in triple_iter
    map<int, int> bead_number_and_bead_ids;
    bead_number_and_bead_ids[1] = ((*triple_iter)->bead1()->getId());
    bead_number_and_bead_ids[2] = ((*triple_iter)->bead2()->getId());
    bead_number_and_bead_ids[3] = ((*triple_iter)->bead3()->getId());

    vector<bool> found_beads(3, false);
    for (const pair<const int, int> &bead_number_and_bead_id :
         bead_number_and_bead_ids) {
      if (bead_number_and_bead_id.second == 0) {
        found_beads.at(0) = true;
      } else if (bead_number_and_bead_id.second == 1) {
        found_beads.at(1) = true;
      } else if (bead_number_and_bead_id.second == 2) {
        found_beads.at(2) = true;
      }
    }

    for (const bool &found : found_beads) {
      BOOST_CHECK(found);
    }

    // Check that the distance between beads correct
    if (bead_number_and_bead_ids[1] == 0 && bead_number_and_bead_ids[2] == 1) {
      // For ids 0 and 1
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 1 &&
               bead_number_and_bead_ids[2] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 0 &&
               bead_number_and_bead_ids[2] == 2) {
      // For ids 0 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 2 &&
               bead_number_and_bead_ids[2] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 1 &&
               bead_number_and_bead_ids[2] == 2) {
      // For ids 1 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 2 &&
               bead_number_and_bead_ids[2] == 1) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
    } else {
      throw runtime_error("Failed to trigger test something is off!");
    }
    if (bead_number_and_bead_ids[1] == 0 && bead_number_and_bead_ids[3] == 1) {
      // For ids 0 and 1
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 1 &&
               bead_number_and_bead_ids[3] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 0 &&
               bead_number_and_bead_ids[3] == 2) {
      // For ids 0 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 2 &&
               bead_number_and_bead_ids[3] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 1 &&
               bead_number_and_bead_ids[3] == 2) {
      // For ids 1 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[1] == 2 &&
               bead_number_and_bead_ids[3] == 1) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
    } else {
      throw runtime_error("Failed to trigger test something is off!");
    }
    if (bead_number_and_bead_ids[2] == 0 && bead_number_and_bead_ids[3] == 1) {
      // For ids 0 and 1
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[2] == 1 &&
               bead_number_and_bead_ids[3] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[2] == 0 &&
               bead_number_and_bead_ids[3] == 2) {
      // For ids 0 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[2] == 2 &&
               bead_number_and_bead_ids[3] == 0) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.414214, 1e-4);
    } else if (bead_number_and_bead_ids[2] == 1 &&
               bead_number_and_bead_ids[3] == 2) {
      // For ids 1 and 2
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);
    } else if (bead_number_and_bead_ids[2] == 2 &&
               bead_number_and_bead_ids[3] == 1) {
      BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);
    } else {
      throw runtime_error("Failed to trigger test something is off!");
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
