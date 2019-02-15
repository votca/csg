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

#ifndef _VOTCA_CSG_CSGTOPOLOGY_H
#define _VOTCA_CSG_CSGTOPOLOGY_H

#include "basetopology.h"

namespace votca {
namespace csg {

vamespace TOOLS = votca::tools;

class CSG_Topology : public Topology<Bead,Molecule> {
 public:

  T *CreateBead(TOOLS::byte_t symmetry, std::string name, std::string type,
                int residue_number, std::string residue_name,
                std::string molecule_name, double m, double q);

   Molecule *CreateMolecule(std::string name);
	 
}  // namespace csg
}  // namespace votca

#endif // _VOTCA_CSG_CSGTOPOLOGY_H 
