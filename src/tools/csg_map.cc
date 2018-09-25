/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectorywriter.h>

using namespace std;
using namespace votca::csg;

class CsgMapApp
    : public CsgApplication
{
public:
    string ProgramName() { return "csg_map"; }
    void HelpText(ostream &out) {
      out << "Convert a reference atomistic trajectory or configuration into a coarse-grained one \n"
          << "based on a mapping xml-file. The mapping can be applied to either an entire trajectory \n"
          << "or a selected set of frames only (see options).\n"
	  << "Examples:\n"
	  << "* csg_map --top FA-topol.tpr --trj FA-traj.trr --out CG-traj.xtc --cg cg-map.xml\n"
	  << "* csg_map --top FA-topol.tpr --trj FA-conf.gro --out CG-conf.gro --cg cg-map.xml\n"
	  << "* csg_map --top FA-topol.tpr --trj FA-traj.xtc --out FA-history.dlph --no-map\n"
	  << "* csg_map --top FA-field.dlpf --trj FA-history.dlph --out CG-history.dlph --cg cg-map.xml\n"
	  << "* csg_map --top .dlpf --trj .dlph --out .dlph --cg cg-map.xml  convert HISTORY to HISTORY_CGV\n";
    }

    bool DoTrajectory() { return true;}
    bool DoMapping() { return true;}

    void Initialize() {
        CsgApplication::Initialize();
        AddProgramOptions()
	  ("out", boost::program_options::value<string>(),"  output file for coarse-grained trajectory")
	  ("vel", "  Write mapped velocities (if available)")
	  ("force", "  Write mapped forces (if available)")
	  ("hybrid", "  Create hybrid trajectory containing both atomistic and coarse-grained");
    }

    bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("trj", "no trajectory file specified");
        CheckRequired("out", "need to specify output trajectory");
        return true;
    }

    void BeginEvaluate(Topology *top, Topology *top_ref);
void EvalConfiguration(Topology *top, Topology *top_ref) {
        if (!_do_hybrid) {
            // simply write the topology mapped by csgapplication class
            if (_do_vel) top->SetHasVel(true);
            if (_do_force) top->SetHasForce(true);
            _writer->Write(top);
        } else {
            // we want to combine atomistic and coarse-grained into one topology
            Topology *hybtol = new Topology();

            hybtol->setBox(top->getBox());
            hybtol->setTime(top->getTime());
            hybtol->setStep(top->getStep());

            // copy all residues from both
            for(auto res : top_ref->Residues() ) hybtol->CreateResidue((res)->getName());
            for(auto res : top->Residues() ) hybtol->CreateResidue((res)->getName());

            // copy all molecules and beads
          
            for(auto mol : top_ref->Molecules() ){
                auto mi = hybtol->CreateMolecule(mol->getName());
                for (int i = 0; i < mol->BeadCount(); i++) {
                    // copy atomistic beads of molecule
                    int beadid = mol->getBead(i)->getId();

                    auto bi = dynamic_pointer_cast<Bead>(mol->getBead(i));
                    auto type = hybtol->GetOrCreateBeadType(bi->getType()->getName());
                    auto bn = hybtol->CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr(), bi->getM(), bi->getQ());
                    bn->setOptions(bi->Options());
                    bn->setPos(bi->getPos());
                    if (bi->HasVel()) bn->setVel(bi->getVel());
                    if (bi->HasF()) bn->setF(bi->getF());

                    mi->AddBead(hybtol->Beads()[beadid]);

                }

                if (mi->getId() < top->MoleculeCount()) {
                    // copy cg beads of molecule
                    auto cgmol = top->Molecules()[mi->getId()];
                    for (int i = 0; i < cgmol->BeadCount(); i++) {
                        auto bi = dynamic_pointer_cast<Bead>(cgmol->getBead(i));
                        // todo: this is a bit dirty as a cg bead will always have the resid of its first parent
                        auto bparent = dynamic_pointer_cast<Bead>(mol->getBead(0));
                        auto type = hybtol->GetOrCreateBeadType(bi->getType()->getName());
                        auto bn = hybtol->CreateBead(bi->getSymmetry(), bi->getName(), type, bparent->getResnr(), bi->getM(), bi->getQ());
                        bn->setOptions(bi->Options());
                        bn->setPos(bi->getPos());
                        if (bi->HasVel()) bn->setVel(bi->getVel());
                        mi->AddBead(bi);
                    }
                }
                
            }
           hybtol->setBox(top_ref->getBox());

            _writer->Write(hybtol);
        }
    }

    void EndEvaluate() {
        _writer->Close();
        delete _writer;
    }

protected:
    TrajectoryWriter *_writer;
    bool _do_hybrid;
    bool _do_vel;
    bool _do_force;

};

void CsgMapApp::BeginEvaluate(Topology *top, Topology *top_atom) {
    string out = OptionsMap()["out"].as<string > ();
    cout << "writing coarse-grained trajectory to " << out << endl;
    _writer = TrjWriterFactory().Create(out);
    if (_writer == NULL)
        throw runtime_error("output format not supported: " + out);

    _do_hybrid = false;
    if(OptionsMap().count("hybrid")){
        if (!_do_mapping)
            throw runtime_error("options hybrid and no-map not compatible");
        cout << "Doing hybrid mapping..." << endl;
        _do_hybrid = true;
    }

    _do_vel = false;
    if(OptionsMap().count("vel")){
        _do_vel = true;
    }

    _do_force = false;
    if(OptionsMap().count("force")){
        _do_force = true;
    }

    _writer->Open(out);
};

int main(int argc, char **argv)
{
    CsgMapApp app;
    return app.Exec(argc, argv);
}

