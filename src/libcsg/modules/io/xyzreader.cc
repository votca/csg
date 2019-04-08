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

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <votca/csg/xyzreader.h>
#include <votca/tools/elements.h>
#include <votca/tools/getline.h>

namespace votca {
namespace csg {

using namespace boost;
using namespace std;
using namespace votca::tools;
/*
bool XYZReader::ReadTopology(string file, CSG_Topology &top) {
  top.Cleanup();

  _file = file;
  _fl.open(file);
  if (!_fl.is_open()) {
    throw std::ios_base::failure("Error on open topology file: " + file);
  }

  ReadFrame<true, CSG_Topology>(top);

  _fl.close();

  return true;
}*/

bool XYZReader::Open(const string &file) {
  _file = file;
  _fl.open(file);
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on open trajectory file: " + file);

  _line = 0;
  return true;
}

void XYZReader::Close() { _fl.close(); }

// bool XYZReader::FirstFrame(CSG_Topology &top) { return NextFrame(top); }
/*
bool XYZReader::NextFrame(CSG_Topology &top) {
  bool success = ReadFrame<false, CSG_Topology>(top);
  return success;
}*/
}  // namespace csg
}  // namespace votca
