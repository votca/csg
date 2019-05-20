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

#define BOOST_TEST_MODULE beadmap_test

#include "../../include/votca/csg/bead.h"
#include "../../include/votca/csg/beadmap.h"
#include "../../include/votca/csg/orthorhombicbox.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/tools/constants.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

class TestBead : public Bead {
 public:
  TestBead(){};
  void setSymmetry(byte_t sym) { symmetry_ = sym; }
  void setElement(string element) { element_symbol_.setName(element); }
};

class TestBeadSphere : public Map_Sphere {
 public:
  TestBeadSphere(){};
};

class TestBeadEllipsoid : public Map_Ellipsoid {
 public:
  TestBeadEllipsoid(){};
};
BOOST_AUTO_TEST_SUITE(beadmap_test)

BOOST_AUTO_TEST_CASE(test_beadmap_constructor) {
  TestBeadSphere test_sphere;
  TestBeadEllipsoid test_ellip;
}

/** @brief Will test that the spherical mapping works
 *
 * Takes 3 beads and maps thier positions, velocities, forces and masses onto a
 * single coarse grained bead.
 */
BOOST_AUTO_TEST_CASE(test_bead_sphere_apply) {

  vector<string> subbeads = {"1:ppn:C2", "1:ppn:H7", "1:ppn:H8"};
  vector<double> weights = {12.0, 1.0, 1.0};
  vector<double> ds;
  TestBeadSphere test_sphere;
  test_sphere.InitializeBeadMap(subbeads, weights, ds);

  vector<string> subbeads_names = test_sphere.getAtomicBeadNames();
  vector<bool> found_bead(3, false);
  BOOST_CHECK_EQUAL(subbeads_names.size(), 3);
  for (string& subbead : subbeads_names) {
    if (subbead == subbeads.at(0)) {
      found_bead.at(0) = true;
    } else if (subbead == subbeads.at(1)) {
      found_bead.at(1) = true;
    } else if (subbead == subbeads.at(2)) {
      found_bead.at(2) = true;
    }
  }

  // Check that all the beads were found
  for (bool found : found_bead) {
    BOOST_CHECK(found);
  }

  Eigen::Matrix3d box;
  box << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0;
  OrthorhombicBox boundaries;
  boundaries.setBox(box);

  TestBead beadC2;

  Eigen::Vector3d vel(0.0, 1.3, -0.1);
  beadC2.setVel(vel);
  Eigen::Vector3d pos(3.4, 2.4, 4.5);
  beadC2.setPos(pos);
  Eigen::Vector3d force(0.3, -0.4, 0.2);
  beadC2.setF(force);
  beadC2.setMass(12.0);
  beadC2.setType("C");
  beadC2.setId(2);

  TestBead beadH7;
  Eigen::Vector3d velH7(0.3, 0.3, 0.1);
  beadH7.setVel(velH7);
  Eigen::Vector3d posH7(3.5, 3.4, 4.5);
  beadH7.setPos(posH7);
  Eigen::Vector3d forceH7(-0.1, -0.2, -0.1);
  beadH7.setF(forceH7);
  beadH7.setMass(1.0);
  beadH7.setType("H");
  beadH7.setId(7);

  TestBead beadH8;
  Eigen::Vector3d velH8(-0.2, 0.4, 0.1);
  beadH8.setVel(velH8);
  Eigen::Vector3d posH8(3.5, 2.4, 2.8);
  beadH8.setPos(posH8);
  Eigen::Vector3d forceH8(0.1, -0.2, -0.1);
  beadH8.setF(forceH8);
  beadH8.setMass(1.0);
  beadH8.setType("H");
  beadH8.setId(8);

  map<string, const Bead*> atomic_beads;
  atomic_beads[subbeads.at(0)] = &beadC2;
  atomic_beads[subbeads.at(1)] = &beadH7;
  atomic_beads[subbeads.at(2)] = &beadH8;

  TestBead cg_bead;
  Eigen::Vector3d initialize(0.0, 0.0, 0.0);
  cg_bead.setType("CH2");
  cg_bead.setF(initialize);
  cg_bead.setVel(initialize);
  cg_bead.setPos(initialize);
  cg_bead.setU(initialize);
  cg_bead.setV(initialize);
  cg_bead.setW(initialize);
  cg_bead.setMass(0.0);
  cg_bead.setSymmetry(1);  // For sphere

  test_sphere.UpdateCGBead(&boundaries, atomic_beads, cg_bead);

  cout << cg_bead.getF() << endl;
  cout << cg_bead.getVel() << endl;
  cout << cg_bead.getPos() << endl;
  cout << cg_bead.getU() << endl;
  cout << cg_bead.getV() << endl;
  cout << cg_bead.getW() << endl;
  cout << cg_bead.getMass() << endl;

  BOOST_CHECK_CLOSE(cg_bead.getF().x(), 0.3, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().y(), -0.8, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().z(), 0.0, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getVel().x(), 0.00714285714285714, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().y(), 1.1642857142857141, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().z(), -0.071428571428571425, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getPos().x(), 3.4142857142857141, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().y(), 2.4714285714285711, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().z(), 4.3785714285714281, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getMass(), 14.0, 1E-5);
}

/** @brief Will test that the ellipsoidal mapping works
 *
 * Takes 3 beads and maps their positions, velocities, forces, orientation and
 * masses onto a single coarse grained bead.
 */
BOOST_AUTO_TEST_CASE(test_bead_ellipsoid_apply) {

  vector<string> subbeads = {"1:ppn:C2", "1:ppn:H7", "1:ppn:H8"};
  vector<double> weights = {12.0, 1.0, 1.0};
  vector<double> ds;
  TestBeadEllipsoid test_ellip;
  test_ellip.InitializeBeadMap(subbeads, weights, ds);

  vector<string> subbeads_names = test_ellip.getAtomicBeadNames();
  vector<bool> found_bead(3, false);
  BOOST_CHECK_EQUAL(subbeads_names.size(), 3);
  for (string& subbead : subbeads_names) {
    if (subbead == subbeads.at(0)) {
      found_bead.at(0) = true;
    } else if (subbead == subbeads.at(1)) {
      found_bead.at(1) = true;
    } else if (subbead == subbeads.at(2)) {
      found_bead.at(2) = true;
    }
  }

  // Check that all the beads were found
  for (bool found : found_bead) {
    BOOST_CHECK(found);
  }

  Eigen::Matrix3d box;
  box << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0;
  OrthorhombicBox boundaries;
  boundaries.setBox(box);

  TestBead beadC2;
  Eigen::Vector3d vel(0.0, 1.3, -0.1);
  beadC2.setVel(vel);
  Eigen::Vector3d pos(3.4, 2.4, 4.5);
  beadC2.setPos(pos);
  Eigen::Vector3d force(0.3, -0.4, 0.2);
  beadC2.setF(force);
  beadC2.setMass(12.0);

  TestBead beadH7;
  Eigen::Vector3d velH7(0.3, 0.3, 0.1);
  beadH7.setVel(velH7);
  Eigen::Vector3d posH7(3.5, 3.4, 4.5);
  beadH7.setPos(posH7);
  Eigen::Vector3d forceH7(-0.1, -0.2, -0.1);
  beadH7.setF(forceH7);
  beadH7.setMass(1.0);

  TestBead beadH8;
  Eigen::Vector3d velH8(-0.2, 0.4, 0.1);
  beadH8.setVel(velH8);
  Eigen::Vector3d posH8(3.5, 2.4, 2.8);
  beadH8.setPos(posH8);
  Eigen::Vector3d forceH8(0.1, -0.2, -0.1);
  beadH8.setF(forceH8);
  beadH8.setMass(1.0);

  map<string, const Bead*> atomic_beads;
  atomic_beads[subbeads.at(0)] = &beadC2;
  atomic_beads[subbeads.at(1)] = &beadH7;
  atomic_beads[subbeads.at(2)] = &beadH8;

  TestBead cg_bead;
  cg_bead.setSymmetry(3);  // For Ellipsoid
  Eigen::Vector3d initialize(0.0, 0.0, 0.0);
  cg_bead.setF(initialize);
  cg_bead.setVel(initialize);
  cg_bead.setPos(initialize);
  cg_bead.setU(initialize);
  cg_bead.setV(initialize);
  cg_bead.setW(initialize);
  cg_bead.setMass(0.0);

  test_ellip.UpdateCGBead(&boundaries, atomic_beads, cg_bead);

  cout << cg_bead.getF() << endl;
  cout << cg_bead.getVel() << endl;
  cout << cg_bead.getPos() << endl;
  cout << cg_bead.getU() << endl;
  cout << cg_bead.getV() << endl;
  cout << cg_bead.getW() << endl;
  cout << cg_bead.getMass() << endl;

  BOOST_CHECK_CLOSE(cg_bead.getF().x(), 0.3, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().y(), -0.8, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().z(), 0.0, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getVel().x(), 0.00714285714285714, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().y(), 1.1642857142857141, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().z(), -0.071428571428571425, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getPos().x(), 3.4142857142857141, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().y(), 2.4714285714285711, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().z(), 4.3785714285714281, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getU().x(), -0.50871738647660514, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getU().y(), 0.84141289136147968, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getU().z(), -0.18229362837775606, 1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getV().x(), 0.0995037, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getV().y(), 0.995037, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getV().z(), 0.0, 1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getW().x(), 0.29377573551808389, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getW().y(), -0.029377573551808416, 1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getW().z(), -0.95542282545112811, 1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getMass(), 14.0, 1E-5);
}

/**
 * @brief Assuming that we have a mapping file of the following format
 *
 * <cg_molecule>
 *   <name>ppn</name> <!-- molecule name in cg representation -->
 *   <ident>propane</ident> <!-- molecule name in atomistic topology -->
 *   <topology> <!-- topology of one molecule -->
 *     <cg_beads>
 *       <cg_bead> <!-- definition of a coarse-grained bead -->
 *         <name>A1</name>
 *         <type>A</type>
 *         <mapping>A</mapping> <!-- reference to a map -->
 *         <!-- atoms belonging to this bead -->
 *         <beads>1:ppn:C1 1:ppn:H4 1:ppn:H5 1:ppn:H6</beads>
 *       </cg_bead>
 * 	     <cg_bead>
 *         <name>B1</name>
 *         <type>B</type>
 *         <symmetry>1</symmetry>
 *         <mapping>B</mapping>
 *         <beads> 1:ppn:C2 1:ppn:H7 1:ppn:H8  </beads>
 *       </cg_bead>
 *       <cg_bead>
 *         <name>A2</name>
 *         <type>A</type>
 *         <symmetry>1</symmetry>
 *         <mapping>A</mapping>
 *         <beads> 1:ppn:C3 1:ppn:H9 1:ppn:H10 1:ppn:H11 </beads>
 *       </cg_bead>
 *     </cg_beads>
 *     <cg_bonded> <!-- bonded interactions -->
 *       <bond>
 *         <name>bond</name>
 *         <beads>
 *           A1 B1
 *           B1 A2
 *         </beads>
 *       </bond>
 *       <angle>
 *         <name>angle</name>
 *         <beads>
 *           A1 B1 A2
 *         </beads>
 *       </angle>
 *     </cg_bonded>
 *   </topology>
 *   <maps>
 *     <map> <!-- mapping A -->
 *       <name>A</name>
 *       <weights> 12 1 1 1 </weights>
 *     </map>
 *     <map>
 *       <name>B</name>
 *       <weights> 12 1 1 </weights>
 *     </map>
 *   </maps>
 * </cg_molecule> <!-- end of the molecule -->
 *
 */
BOOST_AUTO_TEST_CASE(test_atomtocgmapper_apply) {

  CGBeadStencil bead_stencil1;
  bead_stencil1.cg_name_ = "A1";
  bead_stencil1.cg_bead_type_ = "A";
  bead_stencil1.cg_symmetry_ = 1;
  bead_stencil1.mapping_ = "A";
  bead_stencil1.atomic_subbeads_ =
      vector<string>{"1:ppn:C1", "1:ppn:H4", "1:ppn:H5", "1:ppn:H6"};
  bead_stencil1.subbead_weights_ = vector<double>{12, 1, 1, 1};

  CGBeadStencil bead_stencil2;
  bead_stencil2.cg_name_ = "B1";
  bead_stencil2.cg_bead_type_ = "B";
  bead_stencil2.cg_symmetry_ = 1;
  bead_stencil2.mapping_ = "B";
  bead_stencil2.atomic_subbeads_ =
      vector<string>{"1:ppn:C2", "1:ppn:H7", "1:ppn:H8"};
  bead_stencil2.subbead_weights_ = vector<double>{12, 1, 1};

  CGBeadStencil bead_stencil3;
  bead_stencil3.cg_name_ = "A2";
  bead_stencil3.cg_bead_type_ = "A";
  bead_stencil3.cg_symmetry_ = 1;
  bead_stencil3.mapping_ = "A";
  bead_stencil3.atomic_subbeads_ =
      vector<string>{"1:ppn:C3", "1:ppn:H9", "1:ppn:H10", "1:ppn:H11"};
  bead_stencil3.subbead_weights_ = vector<double>{12, 1, 1, 1};

  unordered_map<string, CGBeadStencil> all_molecule_bead_stencils;
  all_molecule_bead_stencils[bead_stencil1.cg_name_] = bead_stencil1;
  all_molecule_bead_stencils[bead_stencil2.cg_name_] = bead_stencil2;
  all_molecule_bead_stencils[bead_stencil3.cg_name_] = bead_stencil3;

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

  Topology atom_top;
  int mol_id = 0;
  string atomic_mol_type = "propane";
  byte_t atom_bead_sym = 1;
  Molecule* atom_mol = &atom_top.CreateMolecule(mol_id, atomic_mol_type);
  int bead_id = 1;
  Bead& C1 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "C", 12, 0.0);
  C1.setPos(pos_c1);
  C1.setVel(initialize);
  C1.setF(initialize);
  atom_mol->AddBead(C1);

  bead_id = 2;
  Bead& C2 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "C", 12, 0.0);
  C2.setPos(pos_c2);
  C2.setVel(initialize);
  C2.setF(initialize);
  atom_mol->AddBead(C2);

  bead_id = 3;
  Bead& C3 = atom_top.CreateBead(atom_bead_sym, "C", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "C", 12, 0.0);
  C3.setPos(pos_c3);
  C3.setVel(initialize);
  C3.setF(initialize);
  atom_mol->AddBead(C3);

  bead_id = 4;
  Bead& H4 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H4.setPos(pos_h4);
  H4.setVel(vel_h4);
  H4.setF(force_h4);
  atom_mol->AddBead(H4);

  bead_id = 5;
  Bead& H5 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H5.setPos(pos_h5);
  H5.setVel(initialize);
  H5.setF(initialize);
  atom_mol->AddBead(H5);

  bead_id = 6;
  Bead& H6 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H6.setPos(pos_h6);
  H6.setVel(initialize);
  H6.setF(initialize);
  atom_mol->AddBead(H6);

  bead_id = 7;
  Bead& H7 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H7.setPos(pos_h7);
  H7.setVel(initialize);
  H7.setF(initialize);
  atom_mol->AddBead(H7);

  bead_id = 8;
  Bead& H8 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H8.setPos(pos_h8);
  H8.setVel(initialize);
  H8.setF(initialize);
  atom_mol->AddBead(H8);

  bead_id = 9;
  Bead& H9 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                 topology_constants::unassigned_residue_id,
                                 topology_constants::unassigned_residue_type,
                                 "H", 12, 0.0);
  H9.setPos(pos_h9);
  H9.setVel(initialize);
  H9.setF(initialize);
  atom_mol->AddBead(H9);

  bead_id = 10;
  Bead& H10 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                  topology_constants::unassigned_residue_id,
                                  topology_constants::unassigned_residue_type,
                                  "H", 12, 0.0);
  H10.setPos(pos_h10);
  H10.setVel(initialize);
  H10.setF(initialize);
  atom_mol->AddBead(H10);

  bead_id = 11;
  Bead& H11 = atom_top.CreateBead(atom_bead_sym, "H", bead_id, mol_id,
                                  topology_constants::unassigned_residue_id,
                                  topology_constants::unassigned_residue_type,
                                  "H", 12, 0.0);
  H11.setPos(pos_h11);
  H11.setVel(vel_h11);
  H11.setF(force_h11);
  atom_mol->AddBead(H11);

  Topology cg_top;
  string cg_mol_type = "propane";
  Molecule* cg_mol = &cg_top.CreateMolecule(mol_id, cg_mol_type);
  int cg_bead_id = 1;
  Bead& A1 = cg_top.CreateBead(
      bead_stencil1.cg_symmetry_, bead_stencil1.cg_bead_type_, cg_bead_id,
      mol_id, topology_constants::unassigned_residue_id,
      topology_constants::unassigned_residue_type,
      topology_constants::unassigned_element, 0.0, 0.0);
  A1.setPos(initialize);
  A1.setVel(initialize);
  A1.setF(initialize);
  cg_mol->AddBead(A1);

  cg_bead_id = 2;
  Bead& B1 = cg_top.CreateBead(
      bead_stencil2.cg_symmetry_, bead_stencil2.cg_bead_type_, cg_bead_id,
      mol_id, topology_constants::unassigned_residue_id,
      topology_constants::unassigned_residue_type,
      topology_constants::unassigned_element, 0.0, 0.0);
  B1.setPos(initialize);
  B1.setVel(initialize);
  B1.setF(initialize);
  cg_mol->AddBead(B1);

  cg_bead_id = 3;
  Bead& A2 = cg_top.CreateBead(
      bead_stencil3.cg_symmetry_, bead_stencil3.cg_bead_type_, cg_bead_id,
      mol_id, topology_constants::unassigned_residue_id,
      topology_constants::unassigned_residue_type,
      topology_constants::unassigned_element, 0.0, 0.0);
  A2.setPos(initialize);
  A2.setVel(initialize);
  A2.setF(initialize);
  cg_mol->AddBead(A2);

  AtomToCGMoleculeMapper mapper(atomic_mol_type, cg_mol_type);

  vector<string> bead_stencil_order = {"A1", "B1", "A2"};
  mapper.InitializeMoleculeMap(all_molecule_bead_stencils, bead_stencil_order);

  pair<int, map<int, vector<pair<string, int>>>>
      cgmolid_cgbeadid_atomname_and_id;

  vector<pair<string, int>> atoms_in_A1;
  atoms_in_A1.push_back(
      make_pair(bead_stencil1.atomic_subbeads_.at(0), C1.getId()));
  atoms_in_A1.push_back(
      make_pair(bead_stencil1.atomic_subbeads_.at(1), H4.getId()));
  atoms_in_A1.push_back(
      make_pair(bead_stencil1.atomic_subbeads_.at(2), H5.getId()));
  atoms_in_A1.push_back(
      make_pair(bead_stencil1.atomic_subbeads_.at(3), H6.getId()));

  vector<pair<string, int>> atoms_in_B1;
  atoms_in_B1.push_back(
      make_pair(bead_stencil2.atomic_subbeads_.at(0), C2.getId()));
  atoms_in_B1.push_back(
      make_pair(bead_stencil2.atomic_subbeads_.at(1), H7.getId()));
  atoms_in_B1.push_back(
      make_pair(bead_stencil2.atomic_subbeads_.at(2), H8.getId()));

  vector<pair<string, int>> atoms_in_A2;
  atoms_in_A2.push_back(
      make_pair(bead_stencil3.atomic_subbeads_.at(0), C3.getId()));
  atoms_in_A2.push_back(
      make_pair(bead_stencil3.atomic_subbeads_.at(1), H9.getId()));
  atoms_in_A2.push_back(
      make_pair(bead_stencil3.atomic_subbeads_.at(2), H10.getId()));
  atoms_in_A2.push_back(
      make_pair(bead_stencil3.atomic_subbeads_.at(3), H11.getId()));

  cgmolid_cgbeadid_atomname_and_id.first = mol_id;
  cgmolid_cgbeadid_atomname_and_id.second[A1.getId()] = atoms_in_A1;
  cgmolid_cgbeadid_atomname_and_id.second[B1.getId()] = atoms_in_B1;
  cgmolid_cgbeadid_atomname_and_id.second[A2.getId()] = atoms_in_A2;

  // This should update the cg beads with the appropriate atomistic bead pos,
  // velocity etc
  mapper.UpdateCGMolecule(atom_top, cg_top, cgmolid_cgbeadid_atomname_and_id);

  BOOST_CHECK_CLOSE(A1.getPos().x(), 0.933333, 1E-4);
  BOOST_CHECK_CLOSE(A1.getPos().y(), 1.0, 1E-4);
  BOOST_CHECK_CLOSE(A1.getPos().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getPos().x(), 2.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getPos().y(), 1.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getPos().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getPos().x(), 3.0666666666666673, 1E-4);
  BOOST_CHECK_CLOSE(A2.getPos().y(), 1.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getPos().z(), 0.0, 1E-4);

  BOOST_CHECK_CLOSE(A1.getVel().x(), -0.033333333333333333, 1E-4);
  BOOST_CHECK_CLOSE(A1.getVel().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A1.getVel().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getVel().x(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getVel().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getVel().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getVel().x(), 0.04, 1E-4);
  BOOST_CHECK_CLOSE(A2.getVel().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getVel().z(), 0.0, 1E-4);

  BOOST_CHECK_CLOSE(A1.getF().x(), -0.3, 1E-4);
  BOOST_CHECK_CLOSE(A1.getF().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A1.getF().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getF().x(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getF().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getF().z(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getF().x(), -0.7, 1E-4);
  BOOST_CHECK_CLOSE(A2.getF().y(), 0.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getF().z(), 0.0, 1E-4);

  BOOST_CHECK_CLOSE(A1.getMass(), 48.0, 1E-4);
  BOOST_CHECK_CLOSE(B1.getMass(), 36.0, 1E-4);
  BOOST_CHECK_CLOSE(A2.getMass(), 48.0, 1E-4);
}

BOOST_AUTO_TEST_SUITE_END()
