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

#define BOOST_TEST_MODULE tabulatedpotential_test
#include <boost/test/unit_test.hpp>

#include "../csg_boltzmann/tabulatedpotential.h"
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include <votca/tools/constants.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

Eigen::VectorXd getColumnFromFile(string file_name, int column) {
  vector<double> data;
  ifstream file;
  file.open(file_name);
  string line;
  if (file.is_open()) {
    while (getline(file, line)) {
      string word;
      istringstream ss(line);
      for (int i = 0; i < column; ++i) {
        ss >> word;
      }
      data.push_back(stod(word));
    }
    file.close();
  }
  Eigen::VectorXd result =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(data.data(), data.size());
  return result;
}

BOOST_AUTO_TEST_SUITE(tabulatedpotential_test)

BOOST_AUTO_TEST_CASE(test_tabulatedpotential_constructor) {
  TabulatedPotential tabulatedpotential;
  BOOST_CHECK_EQUAL(tabulatedpotential.getTemperature(), 300);
  auto iterations_to_smooth = tabulatedpotential.getSmoothIterations();
  BOOST_CHECK_EQUAL(iterations_to_smooth.first, 0);
  BOOST_CHECK_EQUAL(iterations_to_smooth.second, 0);
}

BOOST_AUTO_TEST_CASE(test_register) {
  TabulatedPotential tablulatedpotential;
  map<string, AnalysisTool *> commands;
  tablulatedpotential.Register(commands);
}

BOOST_AUTO_TEST_CASE(test_command) {

  CSG_Topology top;
  BondedStatistics bonded_statistics;
  string interaction_group = "interaction";
  vector<string> interactions;
  interactions.push_back("file_interactions.txt");
  // Setup BondedStatistics Object
  {

    Eigen::Matrix3d box = 20 * Eigen::Matrix3d::Identity();
    top.setBox(box);

    const BoundaryCondition *bc = top.getBoundaryCondition();
    BOOST_CHECK_EQUAL(bc->getBoxType(), BoundaryCondition::typeOrthorhombic);
    // Create three beads
    byte_t symmetry = 1;

    string bead_type = "H2";

    double mass = 0.9;
    double charge = 0.0;

    // Create a bunch of H2 molecules, each bead is considered a molecule they
    // are placed on a regular grid so the table properties can be compared
    // consistently
    cout << "Creating Beads" << endl;
    int number_of_H2_molecules = 0;
    int bead_id = 0;
    int residue_id = 0;
    for (double x = 2.0; x < (box(0, 0) - 2.0); x += 4.0) {
      for (double y = 2.0; y < (box(1, 1) - 2.0); y += 3.0) {
        for (double z = 2.0; z < (box(2, 2) - 2.0); z += 4.0) {
          ++residue_id;

          Eigen::Vector3d bead_pos(x, y, z);
          Bead *bead_ptr = top.CreateBead(
              symmetry, bead_type, bead_id, number_of_H2_molecules, residue_id,
              topology_constants::unassigned_residue_type,
              topology_constants::unassigned_element, mass, charge);
          bead_ptr->setPos(bead_pos);
          ++number_of_H2_molecules;
          ++bead_id;
        }
      }
    }

    cout << "Creating interactions" << endl;
    int molecule_id = 0;
    // Every molecule interacts with every other molecule
    int interaction_id = 0;
    for (int index = 0; index < number_of_H2_molecules; ++index) {
      for (int index2 = index + 1; index2 < number_of_H2_molecules; ++index2) {
        string interaction_group2 = interaction_group;
        interaction_group2 += to_string(index) + "_" + to_string(index2);
        cout << interaction_group2 << endl;
        auto bond = top.CreateInteraction(
            InteractionType::bond, interaction_group2, interaction_id,
            molecule_id, vector<int>{index, index2});
        interactions.push_back("molecule id " + to_string(molecule_id) +
                               ":group name " + interaction_group +
                               to_string(index) + "_" + to_string(index2) +
                               ":group id " + to_string(bond->getGroupId()) +
                               ":index " + to_string(interaction_id));
        cout << "bond index " << bond->getIndex() << " " << bond->getGroup()
             << " " << bond->getGroupId() << " "
             << InteractionTypeToString(bond->getType()) << endl;
      }
    }

    cout << "Evalulating topology" << endl;
    bonded_statistics.BeginCG(&top, nullptr);
    cout << "Eval Config" << endl;
    bonded_statistics.EvalConfiguration(&top, nullptr);
  }  // End of setup
  cout << "Creating DataCollection" << endl;
  DataCollection<double> &bonded_values = bonded_statistics.BondedValues();
  cout << "bonded_values after pulling out of statistics "
       << bonded_values.size() << endl;

  TabulatedPotential tabulatedpotential;
  map<string, AnalysisTool *> commands;
  tabulatedpotential.Register(commands);

  // Test 1
  cout << "Running Test 1" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    string command = "tab";
    cout << "Command with args" << endl;
    cout << "bonded statistics " << bonded_statistics.BondedValues().size()
         << endl;
    auto iter = bonded_values.begin();
    while (iter != bonded_values.end()) {
      cout << (*iter)->getName() << endl;
      ++iter;
    }
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    cout << "Command with interactions" << endl;
    tabulatedpotential.Command(bonded_statistics, command, interactions);
    cout << "Get columns from file" << endl;
    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);
    Eigen::VectorXd column3 = getColumnFromFile(interactions.at(0), 3);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 3.0, 5.86, 8.73, 11.59, 14.46;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 5.16, 1.84, 0.00, 0.36, 5.70;
    Eigen::VectorXd col3_ref(5);
    col3_ref << 1.16, 0.90, 0.26, -1.00, -1.87;

    cout << column1 << endl;
    cout << col1_ref << endl;
    cout << column2 << endl;
    cout << column3 << endl;
    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column3.isApprox(col3_ref, 1e-2), true);

  }  // End of Test 1

  // Test 2
  cout << "Running Test 2" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "smooth_pdf", "2"};
    string command = "tab";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);
    Eigen::VectorXd column3 = getColumnFromFile(interactions.at(0), 3);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 3.0, 5.86, 8.73, 11.59, 14.46;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 0.79, 0.24, 0.00, 0.14, 0.77;
    Eigen::VectorXd col3_ref(5);
    col3_ref << 0.193, 0.1385, 0.018, -0.134, -0.220;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column3.isApprox(col3_ref, 1e-2), true);

  }  // End of Test 2

  // Test 3
  cout << "Running Test 3" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "smooth_pdf", "2"};
    vector<string> arguments3{"set", "smooth_pot", "1"};
    string command = "tab";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, arguments3);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);
    Eigen::VectorXd column3 = getColumnFromFile(interactions.at(0), 3);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 3.0, 5.86, 8.73, 11.59, 14.46;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 0.58, 0.36, 0.26, 0.33, 0.59;
    Eigen::VectorXd col3_ref(5);
    col3_ref << 0.0777308, 0.0567485, 0.00544064, -0.0573399, -0.0897949;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column3.isApprox(col3_ref, 1e-2), true);

  }  // End of Test 3

  // Test 4
  cout << "Running Test 4" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "periodic", "1"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 3.0, 5.86, 8.73, 11.59, 14.46;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 0.0284148, 0.059643, 0.124631, 0.108033, 0.0284148;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);

  }  // End of Test 4

  // Test 5
  cout << "Running Test 5" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "periodic", "1"};
    vector<string> arguments3{"set", "normalize", "0"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, arguments3);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 3.0, 5.86, 8.73, 11.59, 14.46;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 404.0, 848.0, 1772.0, 1536.0, 404.0;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);

  }  // End of Test 5

  // Test 6
  cout << "Running Test 6" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "periodic", "1"};
    vector<string> arguments3{"set", "extend", "0"};
    vector<string> arguments4{"set", "auto", "0"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, arguments3);
    tabulatedpotential.Command(bonded_statistics, command, arguments4);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 0.0, 0.25, 0.5, 0.75, 1;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 1508.0, 1164.0, 436.0, 1452.0, 1508.0;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);
  }  // End of Test 6

  // Test 7
  cout << "Running Test 7" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "scale", "bond"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 0.0, 0.25, 0.5, 0.75, 1;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 19372.0, 18624.0, 1744.0, 2581.33, 19372.0;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);

  }  // End of Test 7

  // Test 8
  cout << "Running Test 8" << endl;
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set", "scale", "angle"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics, command, arguments);
    tabulatedpotential.Command(bonded_statistics, command, arguments2);
    tabulatedpotential.Command(bonded_statistics, command, interactions);

    Eigen::VectorXd column1 = getColumnFromFile(interactions.at(0), 1);
    Eigen::VectorXd column2 = getColumnFromFile(interactions.at(0), 2);

    Eigen::VectorXd col1_ref(5);
    col1_ref << 0.0, 0.25, 0.5, 0.75, 1;
    Eigen::VectorXd col2_ref(5);
    col2_ref << 5593.78, 4704.86, 909.42, 2130.16, 5593.78;

    BOOST_CHECK_EQUAL(column1.isApprox(col1_ref, 1e-2), true);
    BOOST_CHECK_EQUAL(column2.isApprox(col2_ref, 1e-2), true);

  }  // End of Test 8

  top.Cleanup();
}

BOOST_AUTO_TEST_SUITE_END()
