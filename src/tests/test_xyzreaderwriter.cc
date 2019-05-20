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

#define BOOST_TEST_MODULE xyzreaderwriter_test
#include "../../include/votca/csg/topology.h"
#include "../../include/votca/csg/topologyreader.h"
#include "../../include/votca/csg/trajectorywriter.h"
#include <boost/any.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <fstream>
#include <stddef.h>
#include <string>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
namespace votca {
namespace csg {
class Bead;
}  // namespace csg
}  // namespace votca

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(xyzreaderwriter_test)

BOOST_AUTO_TEST_CASE(test_topologyreader) {

  //////////////////////////////////////////////////////////////////////////
  // Create pdb test file
  //////////////////////////////////////////////////////////////////////////
  // Distances are in angstroms by default in .xyz format
  ofstream outfile("Molecule1.xyz");
  outfile << "12" << endl;
  outfile << "benzene example" << endl;
  outfile << "C        0.00000        1.40272        0.00000" << endl;
  outfile << "H        0.00000        2.49029        0.00000" << endl;
  outfile << "C       -1.21479        0.70136        0.00000" << endl;
  outfile << "H       -2.15666        1.24515        0.00000" << endl;
  outfile << "C       -1.21479       -0.70136        0.00000" << endl;
  outfile << "H       -2.15666       -1.24515        0.00000" << endl;
  outfile << "C        0.00000       -1.40272        0.00000" << endl;
  outfile << "H        0.00000       -2.49029        0.00000" << endl;
  outfile << "C        1.21479       -0.70136        0.00000" << endl;
  outfile << "H        2.15666       -1.24515        0.00000" << endl;
  outfile << "C        1.21479        0.70136        0.00000" << endl;
  outfile << "H        2.15666        1.24515        0.00000" << endl;

  outfile.close();
  //////////////////////////////////////////////////////////////////////////

  Elements ele;

  Topology top;
  TopologyReader::RegisterPlugins();
  string str = "Molecule1.xyz";
  unique_ptr<TopologyReader> reader = TopReaderFactory().Create(str);

  boost::any any_ptr(&top);
  reader->ReadTopology(str, any_ptr);
  BOOST_CHECK_EQUAL(reader != NULL, true);
  BOOST_CHECK_EQUAL(top.BeadCount(), 12);

  vector<int> residue_id =
      vector<int>(12, topology_constants::unassigned_residue_id);

  vector<string> bead_type = {"C", "H", "C", "H", "C", "H",
                              "C", "H", "C", "H", "C", "H"};

  vector<string> bead_element = {"C", "H", "C", "H", "C", "H",
                                 "C", "H", "C", "H", "C", "H"};

  vector<string> residue_type =
      vector<string>(12, topology_constants::unassigned_residue_type);

  vector<double> x = {0.00000, 0.00000, -1.21479, -2.15666, -1.21479, -2.15666,
                      0.00000, 0.00000, 1.21479,  2.15666,  1.21479,  2.15666};

  vector<double> y = {1.40272,  2.49029,  0.70136,  1.24515,
                      -0.70136, -1.24515, -1.40272, -2.49029,
                      -0.70136, -1.24515, 0.70136,  1.24515};

  vector<double> z = vector<double>(12, 0.0);

  // convert to nm
  for (size_t index = 0; index < 12; ++index) {
    x.at(index) = x.at(index) / 10;
    y.at(index) = y.at(index) / 10;
  }

  Bead* bd;
  Eigen::Vector3d v;
  for (int i = 0; i < 12; i++) {
    bd = top.getBead(i);
    BOOST_CHECK_EQUAL(bd->getId(), i);
    BOOST_CHECK_EQUAL(bd->getResidueId(), residue_id.at(i));
    BOOST_CHECK_EQUAL(bd->getType(), bead_type.at(i));
    BOOST_CHECK_EQUAL(bd->getElement(), bead_element.at(i));
    BOOST_CHECK_EQUAL(bd->getResidueType(), residue_type.at(i));
    v = bd->getPos();
    BOOST_CHECK_CLOSE(bd->getQ(), 0, 1e-5);
    BOOST_CHECK_CLOSE(v.x(), x.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.y(), y.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.z(), z.at(i), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(test_trajectorywriter) {

  //////////////////////////////////////////////////////////////////////////
  // Create pdb test file
  //////////////////////////////////////////////////////////////////////////
  // Distances are in angstroms by default in .xyz format
  ofstream outfile("Molecule1.xyz");
  outfile << "12" << endl;
  outfile << "benzene example" << endl;
  outfile << "C        0.00000        1.40272        0.00000" << endl;
  outfile << "H        0.00000        2.49029        0.00000" << endl;
  outfile << "C       -1.21479        0.70136        0.00000" << endl;
  outfile << "H       -2.15666        1.24515        0.00000" << endl;
  outfile << "C       -1.21479       -0.70136        0.00000" << endl;
  outfile << "H       -2.15666       -1.24515        0.00000" << endl;
  outfile << "C        0.00000       -1.40272        0.00000" << endl;
  outfile << "H        0.00000       -2.49029        0.00000" << endl;
  outfile << "C        1.21479       -0.70136        0.00000" << endl;
  outfile << "H        2.15666       -1.24515        0.00000" << endl;
  outfile << "C        1.21479        0.70136        0.00000" << endl;
  outfile << "H        2.15666        1.24515        0.00000" << endl;

  outfile.close();
  //////////////////////////////////////////////////////////////////////////

  Elements ele;

  Topology top;
  TopologyReader::RegisterPlugins();
  string str = "Molecule1.xyz";
  unique_ptr<TopologyReader> reader = TopReaderFactory().Create(str);

  boost::any any_ptr(&top);
  reader->ReadTopology(str, any_ptr);
  BOOST_CHECK_EQUAL(reader != NULL, true);
  BOOST_CHECK_EQUAL(top.BeadCount(), 12);

  TrajectoryWriter::RegisterPlugins();
  string str2 = "Molecule2.xyz";
  unique_ptr<TrajectoryWriter> writer = TrjWriterFactory().Create(str2);
  writer->Open(str2);
  writer->Write(&top);
  writer->Close();

  Topology top2;
  unique_ptr<TopologyReader> reader2 = TopReaderFactory().Create(str2);
  boost::any any_ptr2(&top2);
  reader2->ReadTopology(str2, any_ptr2);
  BOOST_CHECK_EQUAL(reader2 != NULL, true);
  BOOST_CHECK_EQUAL(top2.BeadCount(), 12);

  vector<int> residue_id =
      vector<int>(12, topology_constants::unassigned_residue_id);

  vector<string> bead_type = {"C", "H", "C", "H", "C", "H",
                              "C", "H", "C", "H", "C", "H"};

  vector<string> bead_element = {"C", "H", "C", "H", "C", "H",
                                 "C", "H", "C", "H", "C", "H"};

  vector<string> residue_type =
      vector<string>(12, topology_constants::unassigned_residue_type);

  vector<double> x = {0.00000, 0.00000, -1.21479, -2.15666, -1.21479, -2.15666,
                      0.00000, 0.00000, 1.21479,  2.15666,  1.21479,  2.15666};

  vector<double> y = {1.40272,  2.49029,  0.70136,  1.24515,
                      -0.70136, -1.24515, -1.40272, -2.49029,
                      -0.70136, -1.24515, 0.70136,  1.24515};

  vector<double> z = vector<double>(12, 0.0);

  // convert to nm
  for (size_t index = 0; index < 12; ++index) {
    x.at(index) = x.at(index) / 10;
    y.at(index) = y.at(index) / 10;
  }
  Bead* bd;
  Eigen::Vector3d v;
  for (int i = 0; i < 10; i++) {
    bd = top2.getBead(i);
    BOOST_CHECK_EQUAL(bd->getId(), i);
    BOOST_CHECK_EQUAL(bd->getResidueId(), residue_id.at(i));
    BOOST_CHECK_EQUAL(bd->getType(), bead_type.at(i));
    BOOST_CHECK_EQUAL(bd->getElement(), bead_element.at(i));
    BOOST_CHECK_EQUAL(bd->getResidueType(), residue_type.at(i));
    v = bd->getPos();
    BOOST_CHECK_CLOSE(bd->getQ(), 0, 1e-5);
    BOOST_CHECK_CLOSE(v.x(), x.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.y(), y.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.z(), z.at(i), 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
