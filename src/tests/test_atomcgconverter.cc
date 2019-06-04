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

#define BOOST_TEST_MODULE atomcgconverter_test
#include "../../include/votca/csg/atomcgconverter.h"
#include "../../include/votca/csg/molecule.h"
#include "../../include/votca/csg/topology.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/tools/constants.h>
#include <votca/tools/property.h>
using namespace std;
using namespace votca::csg;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(atomcgconverter_test)

BOOST_AUTO_TEST_CASE(test_atomcgconverter_constructor) {
  AtomCGConverter converter(vector<string> empty);
}

BOOST_AUTO_TEST_CASE(test_load) {

  Topology atom_top;
  atom_top.setStep(0);
  // Create two propane molecules

  // Propane atoms will be positioned as
  //
  //      H5   H7   H9
  //      |    |    |
  // H4 - C1 - C2 - C3 - H11
  //      |    |    |
  //      H6   H8   H10
  //
  Eigen::Vector3d pos_c1(1.0, 1.0, 0.0);
  Eigen::Vector3d pos_h5(1.0, 2.0, 0.0);
  Eigen::Vector3d pos_h6(1.0, 0.0, 0.0);
  Eigen::Vector3d pos_h4(0.0, 1.0, 0.0);
  Eigen::Vector3d pos_c2(2.0, 1.0, 0.0);
  Eigen::Vector3d pos_h7(2.0, 2.0, 0.0);
  Eigen::Vector3d pos_h8(2.0, 0.0, 0.0);
  Eigen::Vector3d pos_c3(3.0, 1.0, 0.0);
  Eigen::Vector3d pos_h11(4.0, 1.0, 0.0);
  Eigen::Vector3d pos_h9(3.0, 2.0, 0.0);
  Eigen::Vector3d pos_h10(3.0, 0.0, 0.0);

  Eigen::Vector3d translate(0.0, 5.0, 0.0);
  // Only the H4 and H11 atoms will be given velocities and forces
  //
  // Velocities
  //
  // <- H4     H11 ->
  //
  Eigen::Vector3d vel_h4(-0.5, 0.0, 0.0);
  Eigen::Vector3d vel_h11(0.6, 0.0, 0.0);

  // Forces
  //
  // <- H4     <- H11
  //
  Eigen::Vector3d force_h4(-0.3, 0.0, 0.0);
  Eigen::Vector3d force_h11(-0.7, 0.0, 0.0);

  // Used to initialize all other forces and velocities
  Eigen::Vector3d initialize(0.0, 0.0, 0.0);

  byte_t atom_bead_sym = 1;
  string mol_type = "propane";
  int bead_id = 1;
  Molecule* propane_mol;
  for (int mol_id = 1; mol_id < 3; ++mol_id) {

    propane_mol = &atom_top.CreateMolecule(mol_id, mol_type);

    // cout << "Address original of molecule " << endl;
    // cout << "id " << mol_id << " address " << propane_mol << endl;
    Bead& C1 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "C", 12, 0.0);
    C1.setPos(pos_c1);
    C1.setVel(initialize);
    C1.setF(initialize);
    propane_mol->AddBead(C1);

    ++bead_id;
    Bead& C2 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "C", 12, 0.0);
    C2.setPos(pos_c2);
    C2.setVel(initialize);
    C2.setF(initialize);
    propane_mol->AddBead(C2);

    ++bead_id;
    Bead& C3 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "C", 12, 0.0);
    C3.setPos(pos_c3);
    C3.setVel(initialize);
    C3.setF(initialize);
    propane_mol->AddBead(C3);

    ++bead_id;
    Bead& H4 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H4.setPos(pos_h4);
    H4.setVel(vel_h4);
    H4.setF(force_h4);
    propane_mol->AddBead(H4);

    ++bead_id;
    Bead& H5 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H5.setPos(pos_h5);
    H5.setVel(initialize);
    H5.setF(initialize);
    propane_mol->AddBead(H5);

    ++bead_id;
    Bead& H6 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H6.setPos(pos_h6);
    H6.setVel(initialize);
    H6.setF(initialize);
    propane_mol->AddBead(H6);

    ++bead_id;
    Bead& H7 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H7.setPos(pos_h7);
    H7.setVel(initialize);
    H7.setF(initialize);
    propane_mol->AddBead(H7);

    ++bead_id;
    Bead& H8 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H8.setPos(pos_h8);
    H8.setVel(initialize);
    H8.setF(initialize);
    propane_mol->AddBead(H8);

    ++bead_id;
    Bead& H9 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                   topology_constants::unassigned_residue_id,
                                   topology_constants::unassigned_residue_type,
                                   "H", 12, 0.0);
    H9.setPos(pos_h9);
    H9.setVel(initialize);
    H9.setF(initialize);
    propane_mol->AddBead(H9);

    ++bead_id;
    Bead& H10 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                    topology_constants::unassigned_residue_id,
                                    topology_constants::unassigned_residue_type,
                                    "H", 12, 0.0);
    H10.setPos(pos_h10);
    H10.setVel(initialize);
    H10.setF(initialize);
    propane_mol->AddBead(H10);

    ++bead_id;
    Bead& H11 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                    topology_constants::unassigned_residue_id,
                                    topology_constants::unassigned_residue_type,
                                    "H", 12, 0.0);
    H11.setPos(pos_h11);
    H11.setVel(vel_h11);
    H11.setF(force_h11);
    // propane_mol->AddBead(H11);
    atom_top.getMolecule(mol_id).AddBead(H11);
    ++bead_id;

    pos_c1 += translate;
    pos_h5 += translate;
    pos_h6 += translate;
    pos_h4 += translate;
    pos_c2 += translate;
    pos_h7 += translate;
    pos_h8 += translate;
    pos_c3 += translate;
    pos_h11 += translate;
    pos_h9 += translate;
    pos_h10 += translate;
  }

  BOOST_CHECK_EQUAL(atom_top.BeadCount(), 22);
  BOOST_CHECK_EQUAL(atom_top.MoleculeCount(), 2);

  BOOST_CHECK_EQUAL(atom_top.getMolecule(1).BeadCount(), 11);
  BOOST_CHECK_EQUAL(atom_top.getMolecule(2).BeadCount(), 11);

  string file_cg = "cg_molecule.xml";
  ofstream out_cg(file_cg);

  out_cg << "<cg_molecule>\n";
  out_cg << "  <name>ppn</name> <!-- molecule name in cg representation -->\n";
  out_cg << "  <ident>propane</ident> \n";
  out_cg << "  <topology> <!-- topology of one molecule -->\n";
  out_cg << "   <cg_beads>\n";
  out_cg << "      <cg_bead>\n";
  out_cg << "        <name>b1</name>\n";
  out_cg << "        <type>CH3</type>\n";
  out_cg << "        <symmetry>1</symmetry>\n";
  out_cg << "        <mapping>A</mapping>\n";
  out_cg << "        <beads> 1:ppn:C1 1:ppn:H4 1:ppn:H5 1:ppn:H6 </beads>\n";
  out_cg << "      </cg_bead>\n";
  out_cg << "      <cg_bead>\n";
  out_cg << "        <name>b2</name>\n";
  out_cg << "        <type>CH2</type>\n";
  out_cg << "        <symmetry>1</symmetry>\n";
  out_cg << "        <mapping>B</mapping>\n";
  out_cg << "        <beads> 1:ppn:C2 1:ppn:H7 1:ppn:H8  </beads>\n";
  out_cg << "      </cg_bead>\n";
  out_cg << "      <cg_bead>\n";
  out_cg << "        <name>b3</name>\n";
  out_cg << "        <type>CH3</type>\n";
  out_cg << "        <symmetry>1</symmetry>\n";
  out_cg << "        <mapping>A</mapping>\n";
  out_cg << "        <beads> 1:ppn:C3 1:ppn:H9 1:ppn:H10 1:ppn:H11 </beads>\n";
  out_cg << "      </cg_bead>\n";
  out_cg << "    </cg_beads>\n";
  out_cg << "    <cg_bonded> <!-- bonded interactions -->\n";
  out_cg << "      <bond>\n";
  out_cg << "        <name>bond</name>\n";
  out_cg << "        <beads>\n";
  out_cg << "          b1 b2\n";
  out_cg << "          b2 b3\n";
  out_cg << "        </beads>\n";
  out_cg << "      </bond>\n";
  out_cg << "      <angle>\n";
  out_cg << "        <name>angle</name>\n";
  out_cg << "        <beads>\n";
  out_cg << "          b1 b2 b3\n";
  out_cg << "        </beads>\n";
  out_cg << "      </angle>\n";
  out_cg << "    </cg_bonded>\n";
  out_cg << "  </topology>\n";
  out_cg << "  <maps>\n";
  out_cg << "    <map> <!-- mapping A -->\n";
  out_cg << "      <name>A</name>\n";
  out_cg << "      <weights> 12 1 1 1 </weights>\n";
  out_cg << "    </map>\n";
  out_cg << "    <map>\n";
  out_cg << "      <name>B</name>\n";
  out_cg << "      <weights> 12 1 1 </weights>\n";
  out_cg << "     </map>  \n";
  out_cg << "  </maps>\n";
  out_cg << "</cg_molecule> <!-- end of the molecule -->\n";
  out_cg.close();

  // Create converter
  AtomCGConverter converter;
  converter.LoadMoleculeStencil(file_cg);

  BOOST_CHECK_EQUAL(converter.getCGMoleculeType("propane"), string("ppn"));
  BOOST_CHECK_EQUAL(converter.getAtomisticMoleculeType("ppn"),
                    string("propane"));

  // Create an empty topology
  Topology cg_top;

  cout << "Converting top" << endl;
  cg_top = converter.Convert(atom_top);

  cout << "Conversion success" << endl;
  cg_top.setStep(0);

  vector<int> molecule_ids = cg_top.getMoleculeIds();
  sort(molecule_ids.begin(), molecule_ids.end());
  BOOST_CHECK_EQUAL(molecule_ids.at(0), 1);
  BOOST_CHECK_EQUAL(molecule_ids.at(1), 2);

  vector<int> cg_bead_ids = cg_top.getBeadIds();
  sort(cg_bead_ids.begin(), cg_bead_ids.end());
  BOOST_CHECK_EQUAL(cg_bead_ids.size(), 6);

  vector<int> bead_ids_mol_1 =
      cg_top.getMolecule(molecule_ids.at(0)).getBeadIds();
  BOOST_CHECK_EQUAL(bead_ids_mol_1.size(), 3);
  vector<int> bead_ids_mol_2 =
      cg_top.getMolecule(molecule_ids.at(1)).getBeadIds();
  BOOST_CHECK_EQUAL(bead_ids_mol_2.size(), 3);
  sort(bead_ids_mol_1.begin(), bead_ids_mol_1.end());
  sort(bead_ids_mol_2.begin(), bead_ids_mol_2.end());

  vector<string> bead_types = cg_top.getBeadTypes();
  // Should be CH3 and CH2
  BOOST_CHECK_EQUAL(bead_types.size(), 2);
  sort(bead_types.begin(), bead_types.end());
  BOOST_CHECK_EQUAL(bead_types.at(0), "CH2");
  BOOST_CHECK_EQUAL(bead_types.at(1), "CH3");

  Bead& bead = cg_top.getBead(bead_ids_mol_1.at(0));
  BOOST_CHECK_EQUAL(bead.getType(), "CH3");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 1.1333333333333333, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 1.0, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);

  bead = cg_top.getBead(bead_ids_mol_1.at(1));
  BOOST_CHECK_EQUAL(bead.getType(), "CH2");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 1.0714285714285714, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 1.857142857142857, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);

  bead = cg_top.getBead(bead_ids_mol_1.at(2));
  BOOST_CHECK_EQUAL(bead.getType(), "CH3");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 2.2666666666666671, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 0.2, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);

  bead = cg_top.getBead(bead_ids_mol_2.at(0));
  BOOST_CHECK_EQUAL(bead.getType(), "CH3");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 1.1333333333333333, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 6.0, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);

  bead = cg_top.getBead(bead_ids_mol_2.at(1));
  BOOST_CHECK_EQUAL(bead.getType(), "CH2");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 1.0714285714285714, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 6.85714, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);

  bead = cg_top.getBead(bead_ids_mol_2.at(2));
  BOOST_CHECK_EQUAL(bead.getType(), "CH3");
  cout << bead.getPos() << endl;
  BOOST_CHECK_CLOSE(bead.getPos().x(), 2.2666666666666671, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().y(), 5.2, 1E-4);
  BOOST_CHECK_CLOSE(bead.getPos().z(), 0.0, 1E-4);
}

BOOST_AUTO_TEST_SUITE_END()
