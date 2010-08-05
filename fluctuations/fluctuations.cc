/*                                                                                                                                                    
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)                                                                                   
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

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <boost/program_options.hpp>
#include <votca/csg/csgapplication.h>
#include <votca/tools/average.h>
#include <votca/tools/tokenizer.h>
#include <votca/csg/cgengine.h>
#include <math.h>

//using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgFluctuations
    : public CsgApplication
{
    string ProgramName() { return "fluctuations"; }
    void HelpText(ostream &out) { 
        out << "calculate density fluctuations in a subvolume";
    }

    // some program options are added here

    void Initialize() {
        CsgApplication::Initialize();
        // add program option to pick molecule
        AddProgramOptions("Fluctuation options")
                ("filter", boost::program_options::value<string> (&_filter)->default_value("*"), "filter molecule names")
                ("rmax", boost::program_options::value<double>(), "maximal radial distance from refmol to be considerd")
                ("rmin", boost::program_options::value<double>(&_rmin)->default_value(0.0), "minimal radial distance from refmol to be considerd")
                ("refmol", boost::program_options::value<string>(&_refmol)->default_value(""), "Reference molecule")
                ("nbin", boost::program_options::value<int>(&_nbins)->default_value(100), "Number of bins")
                ("outfile", boost::program_options::value<string>(&_outfilename)->default_value("fluctuations.dat"), "Output file");


    }
     bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("rmax");
        return true;
    }
    
    // we want to process a trajectory
    bool DoTrajectory() {return true;}
    bool DoMapping() { return true;}

    void BeginEvaluate(Topology *top, Topology *top_atom) {
        _filter = OptionsMap()["filter"].as<string>();
        _refmol = OptionsMap()["refmol"].as<string>();
        _rmin = OptionsMap()["rmin"].as<double>();
        _rmax = OptionsMap()["rmax"].as<double>();
        _nbins = OptionsMap()["nbin"].as<int>();
        _outfilename = OptionsMap()["outfile"].as<string>();
        _nframes=0;

        _N_avg = new double [_nbins];
        _N_sq_avg = new double [_nbins];
        N = new int [_nbins];
        for (int i = 0; i < _nbins; i++) {
            _N_avg[i] = 0;
            _N_sq_avg[i] = 0;

        }

        cout << "Calculating fluctions for " << _rmin << "<r<" << _rmax;
        cout << "using " << _nbins << " bins" <<endl;
         _outfile.open(_outfilename.c_str());
        if(!_outfile)
                throw runtime_error("cannot open hist_u.xvg for output");

    }

    // write out results in EndEvaluate
    void EndEvaluate();
    // do calculation in this function
    void EvalConfiguration(Topology *top, Topology *top_ref);

protected:
    // number of particles in dV
    int _nbins;
    double *_N_avg;
    // sqare
    double *_N_sq_avg;
    int * N;
    string _filter;
    string _refmol;
    double _rmax;
    double _rmin;
    vec  _ref;
    int _nframes;
    string _outfilename;
    ofstream _outfile;

};

int main(int argc, char** argv)
{
    CsgFluctuations app;

    return app.Exec(argc, argv);
}

void CsgFluctuations::EvalConfiguration(Topology *conf, Topology*conf_atom = 0)
{
     vec eR;
     double r;
     int rbin;

     
            
    if (_refmol != ""){
            for(BeadContainer::iterator iter = conf->Beads().begin();
            iter!=conf->Beads().end();++iter) {
                Bead *bead = *iter;
                if(wildcmp(_refmol.c_str(), bead->getName().c_str())){
                    _ref = bead->getPos();
                    cout << " Solute pos " << _ref << endl;
                }

            }
        }

     for (int i=0; i< _nbins; i++){
         N[i]=0;
     }

     /* check how many molecules are in each bin*/
    for (BeadContainer::iterator iter = conf->Beads().begin();
            iter != conf->Beads().end(); ++iter) {
        Bead *bead = *iter;
        if (!wildcmp(_filter.c_str(), bead->getName().c_str())) continue;


        eR = bead->getPos() - _ref;
        r = abs(eR);

        if (r > _rmin && r < _rmax) {
            rbin = (int)_nbins*(double)((r - _rmin)/(_rmax - _rmin));
            N[rbin]++;
        }
    }

    /* update averages*/
    for (int i=0; i< _nbins; i++){
        _N_avg[i] += N[i];
        _N_sq_avg[i] += N[i]*N[i];
    }

    _nframes++;
}


// output everything when processing frames is done
void CsgFluctuations::EndEvaluate()
{
    cout << "Writing results to " << _outfilename << endl;
    _outfile << "# radius number_fluct avg_number" << endl;
    
    for (int i=0; i< _nbins; i++){
        _N_avg[i] /=_nframes;
        _N_sq_avg[i] /= _nframes;
    }
    for (int i=0; i< _nbins; i++){
        _outfile << _rmin+i*(_rmax-_rmin)/_nbins << " ";
        _outfile << (_N_sq_avg[i]-_N_avg[i]*_N_avg[i])/_N_avg[i] << " "; //fluctuation
        _outfile << _N_avg[i] << endl;
     }
}



// add our user program options


