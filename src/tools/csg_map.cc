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

#include "../../include/votca/csg/topology.h"
#include <fstream>
#include <stddef.h>
#include <stdexcept>
#include <string>
#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectorywriter.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

class CsgMapApp : public CsgApplication {
 public:
  string ProgramName() { return "csg_map"; }
  void HelpText(ostream &out) {
    out << "Convert a reference atomistic trajectory or configuration into a "
           "coarse-grained one \n"
        << "based on a mapping xml-file. The mapping can be applied to either "
           "an entire trajectory \n"
        << "or a selected set of frames only (see options).\n"
        << "Examples:\n"
        << "* csg_map --top FA-topol.tpr --trj FA-traj.trr --out CG-traj.xtc "
           "--cg cg-map.xml\n"
        << "* csg_map --top FA-topol.tpr --trj FA-conf.gro --out CG-conf.gro "
           "--cg cg-map.xml\n"
        << "* csg_map --top FA-topol.tpr --trj FA-traj.xtc --out "
           "FA-history.dlph --no-map\n"
        << "* csg_map --top FA-field.dlpf --trj FA-history.dlph --out "
           "CG-history.dlph --cg cg-map.xml\n"
        << "* csg_map --top .dlpf --trj .dlph --out .dlph --cg cg-map.xml  "
           "convert HISTORY to HISTORY_CGV\n";
  }

  bool DoTrajectory() { return true; }
  bool DoMapping() { return true; }

  void Initialize() {
    CsgApplication::Initialize();
    AddProgramOptions()("out", boost::program_options::value<string>(),
                        "  output file for coarse-grained trajectory")(
        "vel", "  Write mapped velocities (if available)")(
        "force", "  Write mapped forces (if available)")(
        "hybrid",
        "  Create hybrid trajectory containing both atomistic and "
        "coarse-grained");
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
      hybtol->Copy(*top);
      _writer->Write(hybtol);
    }
  }

  void EndEvaluate() { _writer->Close(); }

 protected:
  unique_ptr<TrajectoryWriter> _writer;
  bool _do_hybrid;
  bool _do_vel;
  bool _do_force;
};

void CsgMapApp::BeginEvaluate(Topology *top, Topology *top_atom) {
  string out = OptionsMap()["out"].as<string>();
  cout << "writing coarse-grained trajectory to " << out << endl;
  _writer = TrjWriterFactory().Create(out);
  if (_writer == NULL)
    throw runtime_error("output format not supported: " + out);

  _do_hybrid = false;
  if (OptionsMap().count("hybrid")) {
    if (!_do_mapping)
      throw runtime_error("options hybrid and no-map not compatible");
    cout << "Doing hybrid mapping..." << endl;
    _do_hybrid = true;
  }

  _do_vel = false;
  if (OptionsMap().count("vel")) {
    _do_vel = true;
  }

  _do_force = false;
  if (OptionsMap().count("force")) {
    _do_force = true;
  }

  _writer->Open(out);
};

int main(int argc, char **argv) {
  CsgMapApp app;
  return app.Exec(argc, argv);
}
