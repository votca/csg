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
#include <votca/tools/histogramnew.h>
#include <votca/tools/tokenizer.h>
#include <votca/csg/csgapplication.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;


class CsgDensityApp
    : public CsgApplication
{
    string ProgramName() { return "csg_density"; }
    void HelpText(ostream &out) { 
        out << "Calculates the mass density distribution along a box axis or radial density profile from reference point";
    }

    // some program options are added here
    void Initialize();
    
    // we want to process a trajectory
    bool DoTrajectory() {return true;}
    bool DoMapping() {return true;}
    bool DoMappingDefault(void) { return false; }
   
    // write out results in EndEvaluate
    void EndEvaluate();
    void BeginEvaluate(Topology *top, Topology *top_atom);
    void EvalConfiguration(Topology *top, Topology *top_ref);

    bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("out", "no output topology specified");
        CheckRequired("trj", "no trajectory file specified");
        return true;
    };


protected:
    string _filter, _out;
    HistogramNew _dist;
    string _dens_type;
    double _rmax;
    int _nbin;
    double _scale;
    double _step;
    int _frames;
    int _nblock;
    int _block_length;
    vec _ref;
    vec _axis;
    string _axisname;
    string _molname;
    double _area;
    void WriteDensity(int nframes,const string &suffix="");
};


int main(int argc, char** argv)
{
    CsgDensityApp app;
    return app.Exec(argc, argv);
}

void CsgDensityApp::BeginEvaluate(Topology *top, Topology *top_atom) {

     matrix box = top->getBox();
     vec a = box.getCol(0);
     vec b = box.getCol(1);
     vec c = box.getCol(2);
    
    _dist.setPeriodic(true);
    _axis=vec(0,0,0);
    _area=0;
    if(_axisname=="x") {
      _axis.setX(1);
      _rmax = abs(a);
      _area= abs(b^c);
    }
    else if(_axisname=="y") {
      _axis.setY(1);
       _rmax = abs(b);
      _area= abs(a^c);
    }
    else if(_axisname=="z") {
      _axis.setZ(1);
      _rmax = abs(c);
      _area= abs(a^b);
    }
    else if(_axisname=="r") {
      _dist.setPeriodic(false);
      _rmax = min(min(abs(a/2), abs(b/2)), abs(c/2));
    } else {
      throw std::runtime_error("unknown axis type");
    }

    if(OptionsMap().count("rmax"))
      _rmax = OptionsMap()["rmax"].as<double>();
    
    if(OptionsMap().count("block-length")){
       _block_length=OptionsMap()["block-length"].as<int>();
    } else {
      _block_length=0;
    }


    if (_axisname=="r") {
      if(!OptionsMap().count("ref"))
         _ref = a/2+b/2+c/2;
       cout << "Using referece point: " << _ref << endl;
    } 
    else if(OptionsMap().count("ref"))
        throw std::runtime_error("reference center can only be used in case of spherical density");
    
    _nbin=(int)floor(_rmax/_step);
    _dist.Initialize(0, _rmax, _nbin);
    
    cout << "rmax: " << _rmax << endl;
    cout << "axis: " << _axisname << endl;
    cout << "Bins: " << _nbin << endl;
    _frames=0;
    _nblock=0;
}  

void CsgDensityApp::EvalConfiguration(Topology *top, Topology *top_ref)
{
    // loop over all molecules
    bool did_something = false;
    for(auto mol : top->Molecules() ){
        if(!wildcmp(_molname.c_str(),mol->getName().c_str())) 
            continue;
        int N = mol->BeadCount();
            for(int i=0; i<N; i++) {
                auto b = dynamic_pointer_cast<Bead>(mol->getBead(i));        
                if(!wildcmp(_filter.c_str(), b->getName().c_str()))
                    continue;
                double r;
                if (_axisname=="r") {
                    r = abs(top->BCShortestConnection(_ref, b->getPos()));
                } else {
                    r = b->getPos() *  _axis;
                }
                if (_dens_type=="mass") {
                  _dist.Process(r, b->getM());
                } else {
                  _dist.Process(r, 1.0);
                }
                did_something = true;
            }
    }
    _frames++;
    if (!did_something) throw std::runtime_error("No molecule in selection");
    if(_block_length != 0) {
        if((_nframes % _block_length)==0) {
            _nblock++;
	    string suffix = string("_") + boost::lexical_cast<string>(_nblock);
	    WriteDensity(_block_length,suffix);
	    _dist.Clear();
	}
    }
}


// output everything when processing frames is done
void CsgDensityApp::WriteDensity(int nframes, const string &suffix)
{
  if (_axisname=="r") {
    _dist.data().y() = _scale/(nframes*_rmax/(double)_nbin *4*M_PI) * 
            _dist.data().y().cwiseQuotient(_dist.data().x().cwiseAbs2());
                  
  } else {
    _dist.data().y() = _scale/((double)nframes * _area * _rmax/ (double)_nbin ) *_dist.data().y();
  }
  _dist.data().Save(_out + suffix);    
}

void CsgDensityApp::EndEvaluate()
{
  if(_block_length == 0) 
    WriteDensity(_frames);
}

// add our user program options
void CsgDensityApp::Initialize()
{
    CsgApplication::Initialize();
    // add program option to pick molecule
    AddProgramOptions("Specific options:")
             ("type", boost::program_options::value<string>(&_dens_type)->default_value("mass"), "density type: mass or number")
             ("axis", boost::program_options::value<string>(&_axisname)->default_value("r"), "[x|y|z|r] density axis (r=spherical)")
             ("step", boost::program_options::value<double>(&_step)->default_value(0.01), "spacing of density")
             ("block-length", boost::program_options::value<int>(), "  write blocks of this length, the averages are cleared after every write")
             ("out", boost::program_options::value<string>(&_out), "Output file")
             ("rmax", boost::program_options::value<double>(), "rmax (default for [r] =min of all box vectors/2, else l )")
             ("scale", boost::program_options::value<double>(&_scale)->default_value(1.0), "scale factor for the density")
             ("molname", boost::program_options::value<string>(&_molname)->default_value("*"), "molname")
             ("filter", boost::program_options::value<string>(&_filter)->default_value("*"), "filter bead names")
             ("ref", boost::program_options::value<vec>(&_ref), "reference zero point");
}

