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

#define BOOST_TEST_MODULE bead_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <math.h>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/topology.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
class Molecule;
}  // namespace csg
}  // namespace votca

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(bead_test)

BOOST_AUTO_TEST_CASE(test_bead_constructor) {

  Topology top;

  string bead_type_name = "C1";
  int symmetry = 1;
  string name = "dummy";
  int resnr = 0;
  string residue_name = "DNA";
  double mass = 1.21;
  double charge = -0.87;

  top.CreateBead<Bead>(symmetry, name, bead_type_name, resnr, residue_name,
                       mass, charge);
}

BOOST_AUTO_TEST_CASE(test_bead_getters) {

  Topology top;

  string bead_type_name = "C1";

  int symmetry = 1;
  string name = "dummy";
  int resnr = 0;
  string residue_name = "DNA";
  double mass = 1.21;
  double charge = -0.87;

  Bead* b = top.CreateBead<Bead>(symmetry, name, bead_type_name, resnr,
                                 residue_name, mass, charge);

  BOOST_CHECK_CLOSE(b->getMass(), mass, 1e-5);
  BOOST_CHECK_CLOSE(b->getQ(), charge, 1e-5);
  BOOST_CHECK_EQUAL(b->getId(), 0);
  BOOST_CHECK_EQUAL(b->getName(), name);
  BOOST_CHECK_EQUAL(b->getResidueNumber(), resnr);
  BOOST_CHECK_EQUAL(b->getResidueName(), residue_name);
  BOOST_CHECK_EQUAL(b->getSymmetry(), symmetry);
}

BOOST_AUTO_TEST_CASE(test_bead_setters) {

  Topology top;

  string bead_type_name = "C1";

  int symmetry = 1;
  string name = "dummy";
  int resnr = 0;
  string residue_name = "DNA";
  double mass = 1.21;
  double charge = -0.87;

  Bead* b = top.CreateBead<Bead>(symmetry, name, bead_type_name, resnr,
                                 residue_name, mass, charge);

  double newMass = 9.4;
  double newCharge = 2.6;

  b->setM(newMass);
  b->setQ(newCharge);

  vec xyz(0.1, 0.2, 0.3);
  b->setPos(xyz);

  vec xyz_vel(-2.0, 0.32, 32.0);
  b->setVel(xyz_vel);

  string molecule_name = "TestMol";
  Molecule* mol = top.CreateMolecule(molecule_name);

  b->setMolecule(mol);

  BOOST_CHECK_CLOSE(b->getMass(), newMass, 1e-5);
  BOOST_CHECK_CLOSE(b->getQ(), newCharge, 1e-5);

  auto new_xyz = b->getPos();
  BOOST_CHECK(new_xyz.isClose(xyz, 3));
  auto new_xyz_vel = b->getVel();
  BOOST_CHECK(new_xyz_vel.isClose(xyz_vel, 3));

  auto mol_new = b->getMolecule();
  bool same = !(molecule_name.compare(mol_new->getName()));
  BOOST_CHECK(same);
}

BOOST_AUTO_TEST_SUITE_END()
