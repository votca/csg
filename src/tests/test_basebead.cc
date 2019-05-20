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

#define BOOST_TEST_MODULE basebead_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <string>
#include <votca/csg/basebead.h>
#include <votca/csg/molecule.h>
#include <votca/csg/topology.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

class TestBead : public BaseBead {
 public:
  TestBead() : BaseBead(){};
};

BOOST_AUTO_TEST_SUITE(basebead_test)

BOOST_AUTO_TEST_CASE(test_basebead_constructor) { TestBead basebead; }

BOOST_AUTO_TEST_CASE(test_basebead_getters_setters) {

  TestBead basebead;
  BOOST_CHECK_CLOSE(basebead.getMass(), 0.0, 1e-5);
  BOOST_CHECK(!basebead.HasPos());

  basebead.setId(0);
  BOOST_CHECK_EQUAL(basebead.getId(), 0);

  basebead.setType("Bead1");
  string name = "Bead1";
  BOOST_CHECK(name == basebead.getType());

  basebead.setMass(1.0);
  BOOST_CHECK_CLOSE(basebead.getMass(), 1.0, 1e-5);

  Eigen::Vector3d xyz(-1.3, 2.9, 9.2);
  basebead.setPos(xyz);
  BOOST_CHECK(basebead.HasPos());
  Eigen::Vector3d xyz2 = basebead.getPos();
  Eigen::Vector3d xyz_ref = xyz;

  BOOST_CHECK_EQUAL(xyz2.isApprox(xyz_ref, 1e-5), true);

  Eigen::Vector3d xyz3 = basebead.Pos();

  BOOST_CHECK_EQUAL(xyz3.isApprox(xyz_ref, 1e-5), true);

  Topology top;
  int molecule_id = 1;
  basebead.setMoleculeId(molecule_id);
  BOOST_CHECK_EQUAL(basebead.getMoleculeId(), molecule_id);
}
BOOST_AUTO_TEST_SUITE_END()
