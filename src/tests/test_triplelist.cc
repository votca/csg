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

#define BOOST_TEST_MODULE triplelist_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/beadtriple.h>
#include <votca/csg/topology.h>
#include <votca/csg/triplelist.h>
#include <votca/tools/constants.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(triplelist_test)

BOOST_AUTO_TEST_CASE(triplelist_constructor) {
  TripleList<Bead *, BeadTriple> triplelist;
}

BOOST_AUTO_TEST_CASE(triplelist_add_triple) {
  TripleList<Bead *, BeadTriple> triplelist;

  Topology top;

  string bead_type = "CG";

  byte_t symmetry = 1;
  int molecule_id = 1;
  int bead_id1 = 1;
  string residue_type = "Residue";
  int residue_id = 0;
  double mass = 1.0;
  double charge = -1.0;

  top.CreateBead(symmetry, bead_type, bead_id1, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);

  int bead_id2 = 2;
  symmetry = 1;
  residue_id = 0;
  mass = 2.0;
  charge = -2.0;

  top.CreateBead(symmetry, bead_type, bead_id2, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);

  int bead_id3 = 3;
  symmetry = 1;
  residue_id = 0;
  mass = 3.0;
  charge = -3.0;

  top.CreateBead(symmetry, bead_type, bead_id3, molecule_id, residue_id,
                 residue_type, topology_constants::unassigned_element, mass,
                 charge);

  Eigen::Vector3d dist12(0.1, 0.2, 0.3);
  Eigen::Vector3d dist13(0.2, 0.4, 0.3);
  Eigen::Vector3d dist23(0.1, 0.2, 0.0);

  BeadTriple *testtriple =
      new BeadTriple(top.getBead(bead_id1), top.getBead(bead_id2),
                     top.getBead(bead_id3), dist12, dist13, dist23);

  // Because we added the testriple to the triplelist there is no need for us
  // to delete it, it will be deleted by triplelist container
  triplelist.AddTriple(testtriple);

  BeadTriple *triplefront, *tripleback;

  triplefront = triplelist.front();
  BOOST_CHECK_CLOSE(triplefront->bead1()->getMass(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(triplefront->bead2()->getMass(), 2.0, 1e-5);
  BOOST_CHECK_CLOSE(triplefront->bead3()->getMass(), 3.0, 1e-5);
  BOOST_CHECK_EQUAL(triplefront->bead1()->getResidueId(), 0);
  BOOST_CHECK_EQUAL(triplefront->bead2()->getResidueId(), 0);
  BOOST_CHECK_EQUAL(triplefront->bead3()->getResidueId(), 0);
  BOOST_CHECK_EQUAL(triplelist.size(), 1);

  tripleback = triplelist.back();
  BOOST_CHECK_CLOSE(tripleback->bead1()->getMass(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(tripleback->bead2()->getMass(), 2.0, 1e-5);
  BOOST_CHECK_CLOSE(tripleback->bead3()->getMass(), 3.0, 1e-5);
  BOOST_CHECK_EQUAL(tripleback->bead1()->getResidueId(), 0);
  BOOST_CHECK_EQUAL(tripleback->bead2()->getResidueId(), 0);
  BOOST_CHECK_EQUAL(tripleback->bead3()->getResidueId(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
