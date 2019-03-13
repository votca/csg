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

using namespace std;
using namespace votca::csg;

class TestBead : public Bead {
  public:
    TestBead() {};
    void setSymmetry(byte_t sym){
      symmetry_ = sym;
    }
};

class TestBeadSphere : public Map_Sphere {
  public:
    TestBeadSphere() {};
};

class TestBeadEllipsoid : public Map_Ellipsoid {
  public:
    TestBeadEllipsoid() {};
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
  vector<double> weights = { 12.0, 1.0, 1.0 }; 
  vector<double> ds;
  TestBeadSphere test_sphere;
  test_sphere.Initialize( subbeads,weights,ds);

  vector<string> subbeads_names = test_sphere.getAtomicBeadNames();
  vector<bool> found_bead(3,false);
  BOOST_CHECK_EQUAL(subbeads_names.size(),3);
  for( string & subbead : subbeads_names ){
    if(subbead == subbeads.at(0) ){
      found_bead.at(0)  = true;
    }else if(subbead == subbeads.at(1)){
      found_bead.at(1)  = true;
    }else if(subbead == subbeads.at(2)){
      found_bead.at(2)  = true;
    }
  }

  // Check that all the beads were found
  for ( bool found : found_bead ){
    BOOST_CHECK(found);
  }

  vec row1(10.0,0.0,0.0);
  vec row2(0.0,10.0,0.0);
  vec row3(0.0,0.0,10.0);
  matrix box(row1,row2,row3);
  OrthorhombicBox boundaries;
  boundaries.setBox(box);
  
  TestBead beadC2;
  vec vel(0.0,1.3,-0.1);
  beadC2.setVel(vel);
  vec pos(3.4,2.4,4.5);
  beadC2.setPos(pos);
  vec force(0.3,-0.4,0.2);
  beadC2.setF(force);
  beadC2.setMass(12.0);

  TestBead beadH7;
  vec velH7(0.3,0.3,0.1);
  beadH7.setVel(velH7);
  vec posH7(3.5,3.4,4.5);
  beadH7.setPos(posH7);
  vec forceH7(-0.1,-0.2,-0.1);
  beadH7.setF(forceH7);
  beadH7.setMass(1.0);

  TestBead beadH8;
  vec velH8(-0.2,0.4,0.1);
  beadH8.setVel(velH8);
  vec posH8(3.5,2.4,2.8);
  beadH8.setPos(posH8);
  vec forceH8(0.1,-0.2,-0.1);
  beadH8.setF(forceH8);
  beadH8.setMass(1.0);

  map<string,Bead*> atomic_beads;
  atomic_beads[subbeads.at(0)] = &beadC2;
  atomic_beads[subbeads.at(1)] = &beadH7;
  atomic_beads[subbeads.at(2)] = &beadH8;

  TestBead cg_bead;
  vec initialize(0.0);
  cg_bead.setF(initialize);
  cg_bead.setVel(initialize);
  cg_bead.setPos(initialize);
  cg_bead.setU(initialize);
  cg_bead.setV(initialize);
  cg_bead.setW(initialize);
  cg_bead.setMass(0.0);
  cg_bead.setSymmetry(1); // For sphere
  
  test_sphere.Apply(&boundaries,atomic_beads,&cg_bead);

  cout << cg_bead.getF() << endl; 
  cout << cg_bead.getVel() << endl; 
  cout << cg_bead.getPos() << endl; 
  cout << cg_bead.getU() << endl; 
  cout << cg_bead.getV() << endl; 
  cout << cg_bead.getW() << endl; 
  cout << cg_bead.getMass() << endl;

  BOOST_CHECK_CLOSE(cg_bead.getF().getX(), 0.3, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().getY(), -0.8, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().getZ(), 0.0, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getVel().getX(),0.08571428,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().getY(),1.71428571,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().getZ(),0.08571428571,1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getPos().getX(),8.91428571,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().getY(),7.02857142,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().getZ(),10.11428571,1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getMass(),14.0,1E-5); 
}

/** @brief Will test that the ellipsoidal mapping works
 *
 * Takes 3 beads and maps their positions, velocities, forces, orientation and
 * masses onto a single coarse grained bead.
 */
BOOST_AUTO_TEST_CASE(test_bead_ellipsoid_apply) {

  vector<string> subbeads = {"1:ppn:C2", "1:ppn:H7", "1:ppn:H8"};
  vector<double> weights = { 12.0, 1.0, 1.0 }; 
  vector<double> ds;
  TestBeadEllipsoid test_ellip;
  test_ellip.Initialize( subbeads,weights,ds);

  vector<string> subbeads_names = test_ellip.getAtomicBeadNames();
  vector<bool> found_bead(3,false);
  BOOST_CHECK_EQUAL(subbeads_names.size(),3);
  for( string & subbead : subbeads_names ){
    if(subbead == subbeads.at(0) ){
      found_bead.at(0)  = true;
    }else if(subbead == subbeads.at(1)){
      found_bead.at(1)  = true;
    }else if(subbead == subbeads.at(2)){
      found_bead.at(2)  = true;
    }
  }

  // Check that all the beads were found
  for ( bool found : found_bead ){
    BOOST_CHECK(found);
  }

  vec row1(10.0,0.0,0.0);
  vec row2(0.0,10.0,0.0);
  vec row3(0.0,0.0,10.0);
  matrix box(row1,row2,row3);
  OrthorhombicBox boundaries;
  boundaries.setBox(box);
  
  TestBead beadC2;
  vec vel(0.0,1.3,-0.1);
  beadC2.setVel(vel);
  vec pos(3.4,2.4,4.5);
  beadC2.setPos(pos);
  vec force(0.3,-0.4,0.2);
  beadC2.setF(force);
  beadC2.setMass(12.0);

  TestBead beadH7;
  vec velH7(0.3,0.3,0.1);
  beadH7.setVel(velH7);
  vec posH7(3.5,3.4,4.5);
  beadH7.setPos(posH7);
  vec forceH7(-0.1,-0.2,-0.1);
  beadH7.setF(forceH7);
  beadH7.setMass(1.0);

  TestBead beadH8;
  vec velH8(-0.2,0.4,0.1);
  beadH8.setVel(velH8);
  vec posH8(3.5,2.4,2.8);
  beadH8.setPos(posH8);
  vec forceH8(0.1,-0.2,-0.1);
  beadH8.setF(forceH8);
  beadH8.setMass(1.0);

  map<string,Bead*> atomic_beads;
  atomic_beads[subbeads.at(0)] = &beadC2;
  atomic_beads[subbeads.at(1)] = &beadH7;
  atomic_beads[subbeads.at(2)] = &beadH8;

  TestBead cg_bead;
  cg_bead.setSymmetry(3); // For Ellipsoid
  vec initialize(0.0);
  cg_bead.setF(initialize);
  cg_bead.setVel(initialize);
  cg_bead.setPos(initialize);
  cg_bead.setU(initialize);
  cg_bead.setV(initialize);
  cg_bead.setW(initialize);
  cg_bead.setMass(0.0);
  
  test_ellip.Apply(&boundaries,atomic_beads,&cg_bead);

  cout << cg_bead.getF() << endl; 
  cout << cg_bead.getVel() << endl; 
  cout << cg_bead.getPos() << endl; 
  cout << cg_bead.getU() << endl; 
  cout << cg_bead.getV() << endl; 
  cout << cg_bead.getW() << endl; 
  cout << cg_bead.getMass() << endl;

  BOOST_CHECK_CLOSE(cg_bead.getF().getX(), 0.3, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().getY(), -0.8, 1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getF().getZ(), 0.0, 1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getVel().getX(),0.08571428,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().getY(),1.71428571,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getVel().getZ(),0.08571428571,1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getPos().getX(),8.91428571,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().getY(),7.02857142,1E-5);
  BOOST_CHECK_CLOSE(cg_bead.getPos().getZ(),10.11428571,1E-5);

  BOOST_CHECK_CLOSE(cg_bead.getU().getX(),-0.729109,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getU().getY(),0.676912,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getU().getZ(),-0.100947695,1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getV().getX(),0.0995037,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getV().getY(),0.995037,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getV().getZ(),0.0,1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getW().getX(),0.1256767,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getW().getY(),-0.012567673,1E-4);
  BOOST_CHECK_CLOSE(cg_bead.getW().getZ(),-0.991992,1E-4);

  BOOST_CHECK_CLOSE(cg_bead.getMass(),14.0,1E-5); 
}

BOOST_AUTO_TEST_CASE(test_atomtocgmapper_apply) {


/*
  string file_name = "cg_molecule.xml";
  ofstream outfile(file_name);

  outfile << "<cg_molecule>\n";
  outfile << "  <name>ppn</name> <!-- molecule name in cg representation -->\n";
  outfile << "  <ident>propane</ident> <!-- molecule name in atomistic "
             "topology -->\n";
  outfile << "  <topology> <!-- topology of one molecule -->\n";
  outfile << "    <cg_beads>\n";
  outfile << "";
  outfile << "      <cg_bead> <!-- definition of a coarse-grained bead -->\n";
  outfile << "        <name>A1</name>\n";
  outfile << "        <type>A</type>\n";
  outfile << "        <mapping>A</mapping> <!-- reference to a map -->\n";
  outfile << "        <!-- atoms belonging to this bead -->\n";
  outfile << "        <beads>1:ppn:C1 1:ppn:H4 1:ppn:H5 1:ppn:H6</beads>\n";
  outfile << "      </cg_bead>\n";
  outfile << "";
  outfile << "	     <cg_bead>\n";
  outfile << "        <name>B1</name>\n";
  outfile << "        <type>B</type>\n";
  outfile << "        <symmetry>1</symmetry>\n";
  outfile << "        <mapping>B</mapping>\n";
  outfile << "        <beads> 1:ppn:C2 1:ppn:H7 1:ppn:H8  </beads>\n";
  outfile << "      </cg_bead>\n";
  outfile << "";
  outfile << "      <cg_bead>\n";
  outfile << "        <name>A2</name>\n";
  outfile << "        <type>A</type>\n";
  outfile << "        <symmetry>1</symmetry>\n";
  outfile << "        <mapping>A</mapping>\n";
  outfile << "        <beads> 1:ppn:C3 1:ppn:H9 1:ppn:H10 1:ppn:H11 </beads>\n";
  outfile << "      </cg_bead>\n";
  outfile << "";
  outfile << "    </cg_beads>\n";
  outfile << "    <cg_bonded> <!-- bonded interactions -->\n";
  outfile << "      <bond>\n";
  outfile << "        <name>bond</name>\n";
  outfile << "        <beads>\n";
  outfile << "          A1 B1\n";
  outfile << "          B1 A2\n";
  outfile << "        </beads>\n";
  outfile << "      </bond>\n";
  outfile << "      <angle>\n";
  outfile << "        <name>angle</name>\n";
  outfile << "        <beads>\n";
  outfile << "          A1 B1 A2\n";
  outfile << "        </beads>\n";
  outfile << "      </angle>\n";
  outfile << "    </cg_bonded>\n";
  outfile << "  </topology>\n";
  outfile << "  <maps>\n";
  outfile << "    <map> <!-- mapping A -->\n";
  outfile << "      <name>A</name>\n";
  outfile << "      <weights> 12 1 1 1 </weights>\n";
  outfile << "    </map>\n";
  outfile << "    <map>\n";
  outfile << "      <name>B</name>\n";
  outfile << "      <weights> 12 1 1 </weights>\n";
  outfile << "    </map>\n";
  outfile << "  </maps>\n";
  outfile << "</cg_molecule> <!-- end of the molecule -->\n";

  outfile.close();

  CGMoleculeStencil cgmoleculestencil;
  cgmoleculestencil.Load(file_name);

  BOOST_CHECK_EQUAL(cgmoleculestencil.getCGType(), "ppn");
  BOOST_CHECK_EQUAL(cgmoleculestencil.getAtomisticType(), "propane");

  CSG_Topology top;
  Molecule* cg_molecule = cgmoleculestencil.CreateMolecule(top);

  BOOST_CHECK_EQUAL(cg_molecule->getType(), "ppn");
  // Because it is the first cg_molecule created the id should be 0
  BOOST_CHECK_EQUAL(cg_molecule->getId(), 0);

  // ids   0 - 1 - 2
  // types A - B - A
  BOOST_CHECK_EQUAL(cg_molecule->BeadCount(), 3);

  // Ids of the beads should be 0, 1, and 2 because they are the first beads
  BOOST_CHECK(cg_molecule->BeadExist(0));
  BOOST_CHECK(cg_molecule->BeadExist(1));
  BOOST_CHECK(cg_molecule->BeadExist(2));

  // The bead types of the cg molecule should be A B and A
  BOOST_CHECK_EQUAL(cg_molecule->getBead(0)->getType(), "A");
  BOOST_CHECK_EQUAL(cg_molecule->getBead(1)->getType(), "B");
  BOOST_CHECK_EQUAL(cg_molecule->getBead(2)->getType(), "A");

  // Determine the molecule the beads are attached to should be 0
  BOOST_CHECK_EQUAL(cg_molecule->getBead(0)->getMoleculeId(), 0);
  BOOST_CHECK_EQUAL(cg_molecule->getBead(1)->getMoleculeId(), 0);
  BOOST_CHECK_EQUAL(cg_molecule->getBead(2)->getMoleculeId(), 0);*/
}

BOOST_AUTO_TEST_SUITE_END()
