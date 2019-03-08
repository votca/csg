/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except top_in compliance with the License.
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

#include <cassert>
#include <unordered_map>
#include <fstream>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <votca/tools/tokenizer.h>
#include "../../include/votca/csg/csgtopology.h"
namespace votca {
namespace csg {

using namespace std;

//namespace po = boost::program_options;

void CGEngine::LoadFiles(string filename){
  Tokenizer tok(filename, ";");
  Tokenizer::iterator iter;

  for (iter = tok.begin(); iter != tok.end(); ++iter) {
    string file = *iter;
    boost::trim(file);
    file_names_.insert(file);
  }
}

unique_ptr<AtomCGConverter> CGEngine::PopulateCGTopology(
    CSG_Topology &atomistic_top_in,
    CSG_Topology &cg_top_out
     ) {

  unique_ptr<AtomCGConverter> converter = unique_ptr<AtomCGConverter>(new AtomCGConverter);
  for ( string file : file_names_ ){
    converter->LoadMoleculeStencil(file); 
  }
  converter->Convert(atomistic_top_in,cg_top_out);

  return converter;
}


}  // namespace csg
}  // namespace votca
