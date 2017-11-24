/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#include <vector>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/tools/getline.h>
#include "pdbtopologyreader.h"

namespace votca { namespace csg {
    using namespace boost;
using namespace std;

bool PDBTopologyReader::ReadTopology(string file,  Topology &top)
{
    top.Cleanup();

    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topology file: " + file);

    string line;
    // Two column vector for storing all bonds
    // 1 - id of first atom
    // 2 - id of second atom
    vector<vector<int>> bond_pairs;
    // Store pointers to every bead
    // WARNING we are assuming in the bead_vec that the indices of the beads
    //         correspond to the order in which they are read in. As in the first
    //         bead read in will be at index 0, etc...
    vector<Bead *> bead_vec;
    ////////////////////////////////////////////////////////////////////////////////
    // Read in information from .pdb file
    ////////////////////////////////////////////////////////////////////////////////
    int bead_count = 0 ;
    while ( std::getline(_fl, line) ){
        if( wildcmp("CRYST1*",line.c_str())){
             string a, b, c, alpha, beta, gamma;
             try {
                                            // 1 -  6       Record name    "CRYST1"
               a    =string(line,(7-1),9);  // 7 - 15       Real(9.3)      a           (Angstroms)
               b    =string(line,(16-1),9); //16 - 24       Real(9.3)      b           (Angstroms)
               c    =string(line,(25-1),9); //25 - 33       Real(9.3)      c           (Angstroms)
               alpha=string(line,(34-1),7); //34 - 40       Real(7.2)      alpha       (degrees)
               beta =string(line,(41-1),7); //41 - 47       Real(7.2)      beta        (degrees)
               gamma=string(line,(48-1),7); //48 - 54       Real(7.2)      gamma       (degrees)
                                            //56 - 66       LString        Space group
                                            //67 - 70       Integer        Z value
            } catch (std::out_of_range& err) {
              throw std::runtime_error("Misformated pdb file in CRYST1 line");
            }
            boost::algorithm::trim(a);
            boost::algorithm::trim(b);
            boost::algorithm::trim(c);
            boost::algorithm::trim(alpha);
            boost::algorithm::trim(beta);
            boost::algorithm::trim(gamma);
	        if ((!wildcmp("90*",alpha.c_str()))||(!wildcmp("90*",alpha.c_str()))||(!wildcmp("90*",alpha.c_str()))){
	         throw std::runtime_error("Non cubical box in pdb file not implemented, yet!");
            }
            double aa = boost::lexical_cast<double>(a)/10.0;
            double bb = boost::lexical_cast<double>(b)/10.0;
            double cc = boost::lexical_cast<double>(c)/10.0;
	        top.setBox(matrix(vec(aa, 0 , 0 ),
	                          vec(0 , bb, 0 ),
                              vec(0 , 0 , cc)));

	    }
        if( wildcmp("CONECT*",line.c_str())){
            vector<string> bonded_atms;
            string atm1;
            // Keep track of the number of bonds
            int num_bonds = 0;
            try {
                // If the CONECT keyword is found then there must be at least
                // two atom identifiers, more than that is optional. 
                
                stringstream ss;
                ss << string(line.substr(6));
                // 1 -  6       Record name    "CONECT"
                //11 -  7       Real(5)        atm1           (ID)
                ss >> atm1; 
                string temp_atm;          
                ss >> temp_atm;
                bonded_atms.push_back(temp_atm);
                num_bonds = 1;
                cout << atm1 << " " << temp_atm << " ";
                ss >> temp_atm;
                // Here we have taken a less rigorous approach to the .pdb files
                // we do not care at this point how large the ids of the atoms are
                // they can be greater than 99,999 with this approach. 
                while(ss){
                    cout << temp_atm << " ";
                    bonded_atms.push_back(temp_atm);
                    num_bonds++;
                    ss >> temp_atm;
                }    
                cout << endl;

            } catch (std::out_of_range& err) {
                throw std::runtime_error("Misformated pdb file in CONECT line\n"+line);
            }

            vector<int> row(2);
            boost::algorithm::trim(atm1);
            int at1   = boost::lexical_cast<int>(atm1);
            row.at(0) = at1;

            for(auto bonded_atm=bonded_atms.begin();bonded_atm!=bonded_atms.end();bonded_atm++){
                int at2 = boost::lexical_cast<int>(*bonded_atm);
                row.at(1)= at2 ;
                // Because every bond will be counted twice in a .pdb file
                // we will only add bonds where the id (atm1) is less than the bonded_atm
                if(at1<at2){
                    cerr << "Adding bond pair " << at1 << " " << at2 << endl;
                    bond_pairs.push_back(row);
                }
            }

	    }

        if( wildcmp("ATOM*",line.c_str()) || wildcmp("HETATM*",line.c_str())){
            
            //      according to PDB format
	        string x,y,z, resNum, resName, atName;
            string charge;
            //string atNum;
            try {
	      /* Some pdb don't include all this, read only what we really need*/
	      /* leave this here in case we need more later*/
              //string recType    (line,( 1-1),6); // str       ,  "ATOM", "HETATM"
              //atNum    =    string(line,( 7-1),6); // int       , Atom serial number
              atName   =    string(line,(13-1),4); // str       , Atom name
              //string atAltLoc   (line,(17-1),1); // char      , Alternate location indicator
              resName  =    string(line,(18-1),3); // str       , Residue name
              //string chainID    (line,(22-1),1); // char      , Chain identifier
              resNum   =    string(line,(23-1),4); // int       , Residue sequence number
              //string atICode    (line,(27-1),1); // char      , Code for insertion of res
              x        =    string(line,(31-1),8); // float 8.3 , x
              y        =    string(line,(39-1),8); // float 8.3 , y
              z        =    string(line,(47-1),8); // float 8.3 , z
              //string atOccup    (line,(55-1),6); // float 6.2 , Occupancy
              //string atTFactor  (line,(61-1),6); // float 6.2 , Temperature factor
              //string segID      (line,(73-1),4); // str       , Segment identifier
              //elem_sym =  string(line,(77-1),2); // str       , Element symbol
              charge   =    string(line,(79-1),2); // str       , Charge on the atom

            } catch (std::out_of_range& err) {
              string err_msg = "Misformated pdb file in atom line # "+ boost::lexical_cast<string>(bead_count)
                               +"\n the correct pdb file format requires 80 characters in width. Furthermore, "
                               +"\n to read the topology in from a .pdb file the following attributes must be "
                               +"\n specified:                                                                "
                               +"\n Atom Name, Residue Name, Residue Number, x, y, z, charge (optional)     \n";
              throw std::runtime_error(err_msg);
            }
            boost::algorithm::trim(atName);
            boost::algorithm::trim(resName);
            boost::algorithm::trim(resNum);
            boost::algorithm::trim(x);
            boost::algorithm::trim(y);
            boost::algorithm::trim(z);
            boost::algorithm::trim(charge);
            
            bead_count++;

            int resnr;
            try {
                resnr = boost::lexical_cast<int>(resNum);
            } catch(bad_lexical_cast &) {
                throw std::runtime_error("Cannot convert resNum='"+ resNum+"' to int, that usallly means: misformated pdb file");
            }
            if (resnr < 1)
                throw std::runtime_error("Misformated pdb file, resnr has to be > 0");
            //TODO: fix the case that resnr is not in ascending order
            if(resnr > top.ResidueCount()) {
                while ((resnr-1)>top.ResidueCount()){ //pdb resnr should start with 1 but accept sloppy files
                    top.CreateResidue("DUMMY"); // create dummy residue, hopefully it will never show
                    cout << "Warning: residue numbers not continous, create DUMMY residue with nr " << top.ResidueCount() << endl;
                }
                top.CreateResidue(resName);
            }
            // This is not correct, but still better than no type at all!
            BeadType *type = top.GetOrCreateBeadType(atName);

            // Determine if the charge has been provided in the .pdb file or if we will
            // be assuming it is 0
            double ch=0;
            if(charge!=""){
                ch = boost::lexical_cast<double>(charge);
            }
            // CreateBead takes 6 parameters in the following order
            // 1 - symmetry of the bead (1-indicates sphere, 3-indicates ellipsoidal)
            // 2 - name of the bead     (string)
            // 3 - bead type            (BeadType *)
            // 4 - residue number       (int)
            // 5 - mass                 (double)
            // 6 - charge               (double)
            //
            // res -1 as internal number starts with 0
            Bead *b;
            b = top.CreateBead(1, atName, type, resnr-1,this->getMass(atName),ch);

            // convert to nm from A
            b->setPos(vec(
                        boost::lexical_cast<double>(x)/10.0,
                        boost::lexical_cast<double>(y)/10.0,
                        boost::lexical_cast<double>(z)/10.0
                        ));

            bead_vec.push_back(b);

        }

        if (( line == "ENDMDL" ) || ( line == "END" ) || ( _fl.eof())){
            break;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Sort data and determine atom structure, connect with top (molecules, bonds)
    ////////////////////////////////////////////////////////////////////////////////

    // Now we need to add the bond pairs
    // Step 1 create a vector each element is associated with an atom 
    //        and indicates what molecule the atom is attached too.
    // WARNING We are assuming the atom ids are contiguous with no gaps
    vector<int> atm_molecule(bead_count,-1);
    // Step 2 Create a vector of molecules
    // Each index of the vector will represent a molecule
    // Each index will contain a link list that contains the ids of the atms in the molecule
    vector<list<int>> molecule_atms;
    // Keep track of the number of molecules we have created through an index
    int mol_index = 0;

    // Step 3 Cycle through all bonds
    for(auto row=bond_pairs.begin();row!=bond_pairs.end();row++){
        
        // Step 4 Check to see if either atm referred to in the bond is alread
        //        attached to a molecule
        cerr << "Getting atoms from row" << endl;
        int atm_id1  = row->at(0);
        int atm_id2  = row->at(1);
        int mol_atm1;
        int mol_atm2;
        cerr << "Sorting bond with atoms " << atm_id1 << " " << atm_id2 << endl; 
        try {
            mol_atm1 = atm_molecule.at(atm_id1-1);
            mol_atm2 = atm_molecule.at(atm_id2-1);
        } catch(const std::out_of_range & err) {
            string err_msg = "One of the atoms in the bond: "+boost::lexical_cast<string>(atm_id1)+ 
                             "  "+boost::lexical_cast<string>(atm_id2)+" does not exist\n"+
                             "Keep in mind that the atom with the largest id has a value of "+
                             boost::lexical_cast<string>(atm_molecule.size());
            throw std::runtime_error(err_msg);
        }
        // This means neither atom is attached to a molecule
        if(mol_atm1==-1 && mol_atm2==-1){
            // We are going to create a new row for a new molecule
            list<int> atms_in_mol;
            atms_in_mol.push_back(atm_id1);
            atms_in_mol.push_back(atm_id2);
            molecule_atms.push_back(atms_in_mol);
            // Associate atm1 and atm2 with the molecule index
            atm_molecule.at(atm_id1-1)=mol_index;
            atm_molecule.at(atm_id2-1)=mol_index;
            
            // Increment the molecule index
            cerr << "Creating molecule " << mol_index << " with atms " << atm_id1 << " " << atm_id2 << endl;
            mol_index++;
        // This means only atm2 is attached to a molecule
        }else if(mol_atm1==-1){
            // Add atm1 to the molecule that contains atm2
            molecule_atms.at(mol_atm2).push_back(atm_id1);
            // Associate atm1 with the molecule it is now part of
            atm_molecule.at(atm_id1-1)=mol_atm2;
            cerr << "Adding atom " << atm_id1 << " to molecule id " << mol_atm2 << endl;
        // This means only atm1 is attached to a molecule
        }else if(mol_atm2==-1){
            // Add atm2 to the molecule that contains atm1
            molecule_atms.at(mol_atm1).push_back(atm_id2);
            // Associate atm1 with the molecule it is now part of
            atm_molecule.at(atm_id2-1)=mol_atm1;
            cerr << "Adding atom " << atm_id2 << " to molecule id " << mol_atm1 << endl;

        }else if(mol_atm1!=mol_atm2){
        // This means both atm1 and atm2 are attached to a molecule     
        // But if they are already attached to the same molecule there is 
        // nothing else to be done. 
            int chosen_mol;
            int obsolete_mol;
            // We will merge the atms to the molecule with the smallest index
            if(mol_atm1<mol_atm2){
                chosen_mol   = mol_atm1;
                obsolete_mol = mol_atm2;
            }else if(mol_atm2<mol_atm1){
                chosen_mol   = mol_atm2;
                obsolete_mol = mol_atm1;
            }
            
            cerr << "Joining molecules " << chosen_mol << " and " << obsolete_mol << endl;
            // Grab the atoms from the obsolete molecule
            list<int> obs_mol_atms = molecule_atms.at(obsolete_mol);
            // We will clear out the atms 
            molecule_atms.at(obsolete_mol).clear();
            // Now we will proceed to cycle through the atms that were in the now
            // obsolete molecule
            for(auto atm_temp=obs_mol_atms.begin();atm_temp!=obs_mol_atms.end();atm_temp++){
                // Add the atms from the obsolete molecule to the chosen molecule

                molecule_atms.at(chosen_mol).push_back(*atm_temp);
                // Make sure the newly added atoms are now pointing at the chosen molecule
                atm_molecule.at(*atm_temp-1)=chosen_mol;
            }
        }
    }
 
    int i=0;
    for(auto lis=molecule_atms.begin();lis!=molecule_atms.end();lis++){
        cout << "Molecule " << i << endl;
        cout << "Atoms: ";
        for(auto atm_ind=lis->begin();atm_ind!=lis->end();atm_ind++){
            cout << *atm_ind << " ";
        }
        cout << endl;
    }
    cout << endl;     
    // Now that we know which interactions belong to which molecules we can:
    // 1 Add the molecules
    // 2 Add the bond interactions 

    // Molecule Vector
    vector<Molecule *> mol_vec;   
    // Used to keep track of molecules as they are re-indexed. Re-indexing
    // is neccessary because in the process of sorting the atoms into molecules
    // , some molecules were joined together. The atoms from one of the molecules
    // are all moved to the other molecule. This left a molecule with no atoms in
    // it. Hence, we do not want to record these empty molecules. 
    vector<int> mol_vec_new_ind;   
    cerr << "Creating residues" << endl;
    // Cycle through the molecules
    for(int ind=0, mol_ind=0; mol_ind < molecule_atms.size();ind++){
        string mol_name = "PDB Molecule "+boost::lexical_cast<string>(mol_ind); 
        // Before adding a molecule ensure that it is not an empty molecule
        // It must contain atoms to be a valid molecule object, we will re-index
        // the molecules starting at 0
        if(molecule_atms.at(ind).size()>0){
            cerr << "Creating molecule with name " << mol_name << endl;
            Molecule *mi = top.CreateMolecule(mol_name);
            mol_vec.push_back(mi);
            // Add all the atoms to the appropriate molecule object
            list<int> atm_list = molecule_atms.at(ind);
            cerr << "Adding atoms to the molecule " << endl;
            for(auto atm_temp = atm_list.begin();atm_temp!=atm_list.end();atm_temp++){
                string residuename = "DUM";
                cerr << *atm_temp << " ";
                mi->AddBead(bead_vec.at(*atm_temp-1),residuename);
            }
            cerr << endl;
            mol_vec_new_ind.push_back(mol_ind);
            mol_ind++;
        }else{
            // If the molecule has no atoms given an index of -1
            mol_vec_new_ind.push_back(-1);
        }
    }
    cout << "Adding bonds" << endl;
    int bond_indx = 0;
    // Cyle through the bonds and add them to the appropriate molecule
    for(auto bond_pair=bond_pairs.begin();bond_pair!=bond_pairs.end();bond_pair++){
        // Should be able to just look at one of the atoms the bond is attached too
        // because the other will also be attached to the same molecule. 
        int atm_id1 = bond_pair->at(0);
        int atm_id2 = bond_pair->at(1);
        cout << "Atom 1 " << atm_id1 << " Atom 2 " << atm_id2 << endl;
        int mol_ind  = atm_molecule.at(atm_id1-1);
        cout << "Molecule index " << mol_ind << endl;
        Molecule *mi = mol_vec.at(mol_ind);      
        int new_mol_ind = mol_vec_new_ind.at(mol_ind);  

        // Grab the id of the bead associated with the atom
        // It may be the case that the atom id's and bead id's are different
        int bead_id1 = bead_vec.at(atm_id1-1)->getId();
        int bead_id2 = bead_vec.at(atm_id2-1)->getId();
        Interaction *ic = new IBond(bead_id1,bead_id2);
        ic->setGroup("BONDS");
        ic->setIndex(bond_indx);
        bond_indx++;
        ic->setMolecule(new_mol_ind);  
        top.AddBondedInteraction(ic);
        mi->AddInteraction(ic); 
    }
    cout << "RebuildExclusions" << endl;
    // Finally we want to build an exclusion matrix 
    top.RebuildExclusions();
    _fl.close();
    return true;
}

}}

