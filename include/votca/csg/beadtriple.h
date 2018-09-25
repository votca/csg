/* 
 * Copyright 2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef BEADTRIPLE_H
#define	BEADTRIPLE_H

#include <tuple>
#include <memory>

namespace votca { namespace csg {

/**
   \brief A particle pair
 
   This class defines a particle pair. The future plan is, that the Pair class
   can be overloaded and Particle list creates these inherited pairs.
 
 */

class BeadTriple   
    : public std::tuple<std::shared_ptr<Bead>,std::shared_ptr<Bead>,std::shared_ptr<Bead>>
{
public:
    BeadTriple() {}
    BeadTriple(std::shared_ptr<Bead> bead1, std::shared_ptr<Bead> bead2, std::shared_ptr<Bead> bead3, vec r12, vec r13, vec r23)    
        : std::tuple<std::shared_ptr<Bead>, std::shared_ptr<Bead>, std::shared_ptr<Bead>>(bead1, bead2, bead3), _r12(r12), _r13(r13), _r23(r23), _dist12(abs(r12)), _dist13(abs(r13)), _dist23(abs(r23)) {}
        
    virtual ~BeadTriple() {}

    /// \brief the vector connecting two beads
    vec &r12() { return _r12; }
    vec &r13() { return _r13; }
    vec &r23() { return _r23; }
    /// \brief the distance of the beads
    double &dist12() { return _dist12; }
    double &dist13() { return _dist13; }
    double &dist23() { return _dist23; }    

protected:
    vec _r12;
    vec _r13;
    vec _r23;
    double _dist12;
    double _dist13;
    double _dist23;
};

}}

#endif	/* BEADTRIPLE_H */

