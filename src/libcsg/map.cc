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

#include <votca/csg/map.h>
#include <votca/csg/topology.h>
#include <iostream>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <numeric>
#include <boost/lexical_cast.hpp>

namespace votca { namespace csg {

using namespace std;

Map::~Map()
{
    vector<BeadMap *>::iterator iter;
    for(iter=_maps.begin();iter!=_maps.end();++iter)
        delete (*iter);    
    _maps.clear();
}

void Map::Apply()
{
    for(auto bm : _maps) bm->Apply();
}

void Map_Sphere::Initialize(shared_ptr<Molecule> in, shared_ptr<Bead> out, Property *opts_bead, Property *opts_map) {
    BeadMap::Initialize(in, out, opts_bead, opts_map);
    
    vector<string> beads;
    vector<double> weights;
    vector<double> fweights;

    // get the beads
    string s(_opts_bead->get("beads").value());
    Tokenizer tok_beads(s, " \n\t");
    tok_beads.ToVector(beads);

    // get vector of weights
    Tokenizer tok_weights(_opts_map->get("weights").value(), " \n\t");
    tok_weights.ConvertToVector<double>(weights);

    // check weather weights and # beads matches
    if (beads.size() != weights.size())
        throw runtime_error(string("number of subbeads in " +
                opts_bead->get("name").as<string>()
                + " and number of weights in map "
                + opts_map->get("name").as<string>() + " do not match"));

    // normalize the weights
    double norm = 1./ std::accumulate(weights.begin(), weights.end(), 0.);
    
    transform(weights.begin(), weights.end(), weights.begin(), bind2nd(multiplies<double>(), norm));
    // get the d vector if exists or initialize same as weights
    vector<double> d;
    if(_opts_map->exists("d")) {
        Tokenizer tok_weights(_opts_map->get("d").value(), " \n\t");
        tok_weights.ConvertToVector(d);
        // normalize d coefficients
        norm = 1./std::accumulate(d.begin(), d.end(), 0.);
        transform(d.begin(), d.end(), d.begin(), bind2nd(multiplies<double>(), norm));
    } else {
        // initialize force-weights with weights
        d.resize(weights.size());
        copy(weights.begin(), weights.end(), d.begin());
    }

    // check weather number of d coeffs is correct
    if (beads.size() != d.size()) {
        throw runtime_error(string("number of subbeads in " +
            opts_bead->get("name").as<string>()
            + " and number of d-coefficients in map "
            + opts_map->get("name").as<string>() + " do not match"));
    }

    fweights.resize(weights.size());
    // calculate force weights by d_i/w_i
    for(size_t i=0; i<weights.size(); ++i) {
        if(weights[i] == 0 && d[i]!=0) {
            throw runtime_error(
                "A d coefficient is nonzero while weights is zero in mapping "
                + opts_map->get("name").as<string>());
        }
        if(weights[i] != 0)
            fweights[i] = d[i] / weights[i];
        else
            fweights[i] = 0;        
    }

    for (size_t i = 0; i < beads.size(); ++i) {
      vector<string> bead_name_components;
      Tokenizer tok_bead_name_components(beads[i], ":");
      tok_bead_name_components.ToVector(bead_name_components);
      string beadname = bead_name_components.at(2);
      cout << "Size of vector of string of beads " << beads.size() << endl;
      //vector<int> bead_ids = in->getIdsOfBeadsWithName(beads[i]);
      vector<int> bead_ids = in->getIdsOfBeadsWithName(beadname);
      for(auto v : beads){
        cout << v << endl;
      }
      if (bead_ids.size() == 0){
        //string err = "mapping error: bead " + beads[i] + " does not exist";
        string err = "mapping error: bead " + beadname + " does not exist";
        throw std::runtime_error(err);
      }else if(bead_ids.size()>1){
        string err = "impossible to resolve correct bead from name " + 
          beads[i] + " as more than one bead has the name.";
      }
      auto basebead = in->getBead(bead_ids.at(0));
      auto bead = dynamic_pointer_cast<Bead>(basebead);
      AddElem(bead, weights[i], fweights[i]);
    }
}

void Map_Sphere::Apply()
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    bool bPos, bVel, bF;
    bPos=bVel=bF=false;
    _out->ParentBeads().clear();

    // the following is needed for pbc treatment
    auto top = _out->getParent();
    double max_dist = 0.5*top->ShortestBoxSize();
    vec r0 = vec(0,0,0);
    string name0;
    int id0=0;
    if(_matrix.size() > 0) {
        if(_matrix.front()._in->HasPos()) {
            r0=_matrix.front()._in->getPos();
            name0 = _matrix.front()._in->getName();
            id0 = _matrix.front()._in->getId();
        }
    }

    double M = 0;

    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        auto bead = iter->_in;
        _out->ParentBeads().push_back(bead->getId());
        M+=bead->getM();
        if(bead->HasPos()) {
            vec r = top->BCShortestConnection(r0, bead->getPos());
            if(abs(r) > max_dist) {
                cout << r0 << " " << bead->getPos() << endl;
                throw std::runtime_error("coarse-grained bead is bigger than half the box \n (atoms "
                        + name0 + " (id " + boost::lexical_cast<string>(id0+1) + ")" + ", " + bead->getName() + " (id " + boost::lexical_cast<string>(bead->getId()+1) + ")" +  +" , molecule "
                        + boost::lexical_cast<string>(bead->getMolecule()->getId()+1) + ")" );
            }
            cg += (*iter)._weight * (r+r0);
            bPos=true;
        }
        if(bead->HasVel()) {
            vel += (*iter)._weight * bead->getVel();
            bVel = true;
        }
        if(bead->HasF()) {
            f += (*iter)._force_weight * bead->getF();
            bF = true;
        }
    }
    _out->setM(M);
    if(bPos)
        _out->setPos(cg);
    if(bVel)
        _out->setVel(vel);
    if(bF)
        _out->setF(f);
}

/// \todo implement this function
void Map_Ellipsoid::Apply()
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), c(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    matrix m(0.);
     bool bPos, bVel, bF;
    bPos=bVel=bF=false;

    // the following is needed for pbc treatment
    auto top = _out->getParent();
    double max_dist = 0.5*top->ShortestBoxSize();
    vec r0 = vec(0,0,0);
    if(_matrix.size() > 0) {
        if(_matrix.front()._in->HasPos()) {            
            r0=_matrix.front()._in->getPos();
        }
    }

    int n;
    n = 0;
    _out->ParentBeads().clear();
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
       auto bead = iter->_in;
       _out->ParentBeads().push_back(bead->getId());
       if(bead->HasPos()) {
            vec r = top->BCShortestConnection(r0, bead->getPos());
            if(abs(r) > max_dist)
                throw std::runtime_error("coarse-grained bead is bigger than half the box");
            cg += (*iter)._weight * (r+r0);
            bPos=true;
        }
        if(bead->HasVel() == true) {
            vel += (*iter)._weight * bead->getVel();
            bVel = true;
        }
        if(bead->HasF()) {
            /// \todo fix me, right calculation should be F_i = m_cg / sum(w_i) * sum(w_i/m_i*F_i)
            //f += (*iter)._weight * _in->getBeadF((*iter)._in);
            f += (*iter)._force_weight * bead->getF();
            bF = true;
        }
        
        if((*iter)._weight>0 && bead->HasPos()) {
            c += bead->getPos();
            n++;
        }
    }
    
    if(bPos)
        _out->setPos(cg);
    if(bVel)
        _out->setVel(vel);
    if(bF)
        _out->setF(f);

    if(!_matrix[0]._in->HasPos()) {
        _out->setU(vec(1.0,0,0));
        _out->setV(vec(.0,1,0));
        _out->setW(vec(.0,0,1));
        return;
    }
    
    // calculate the tensor of gyration
    c=c/(double)n;    
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
      if((*iter)._weight == 0) continue;
      auto bead = iter->_in;
      vec v = bead->getPos() - c;
      //v = vec(1, 0.5, 0) * 0.*(drand48()-0.5)
      //    + vec(0.5, -1, 0) * (drand48()-0.5)
      //    + vec(0, 0, 1) * (drand48()-0.5);

      //Normalize the tensor with 1/number_of_atoms_per_bead
      m[0][0] += v.getX()*v.getX()/(double)_matrix.size();
      m[0][1] += v.getX()*v.getY()/(double)_matrix.size();
      m[0][2] += v.getX()*v.getZ()/(double)_matrix.size();
      m[1][1] += v.getY()*v.getY()/(double)_matrix.size();
      m[1][2] += v.getY()*v.getZ()/(double)_matrix.size();
      m[2][2] += v.getZ()*v.getZ()/(double)_matrix.size();

    }
    m[1][0] = m[0][1];
    m[2][0] = m[0][2];
    m[2][1] = m[1][2];
    
    // calculate the eigenvectors
    matrix::eigensystem_t es;
    m.SolveEigensystem(es);
    
    vec u = es.eigenvecs[0];
    vec v = _matrix[1]._in->getPos() - _matrix[0]._in->getPos();
    v.normalize();
    
    _out->setV(v);
    
    vec w = _matrix[2]._in->getPos() - _matrix[0]._in->getPos();
    w.normalize();
    
    if((v^w)*u < 0) u=vec(0.,0.,0.)-u;
    _out->setU(u);
    
    //write out w
    w=u^v;
    w.normalize();
    _out->setW(w);
    
}

}}
