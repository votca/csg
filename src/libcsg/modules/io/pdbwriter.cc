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

#include <boost/format.hpp>
#include <stdio.h>
#include <string>
#include <votca/csg/pdbwriter.h>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;
void PDBWriter::Open(string file, bool bAppend) {
  if (bAppend) {
    _out.open(file, std::ios_base::app);
  } else {
    _out.open(file);
  }
}

void PDBWriter::WriteHeader(std::string header) {
  if (_out.is_open()) {
    if (header.size() < 10 || header.substr(0, 10) != "HEADER    ") {
      _out << "HEADER    ";
    }
    _out << header;
    if (header.back() != '\n') _out << "\n";
  } else {
    throw runtime_error("Cannot write header to pdb file, file is not open.");
  }
}

void PDBWriter::WriteBox(const Eigen::Matrix3d &box) {
  boost::format boxfrmt("CRYST1%1$9.3f%2$9.3f%3$9.3f%4$7.2f%5$7.2f%6$7.2f\n");
  double a = box.col(0).norm();
  double b = box.col(1).norm();
  double c = box.col(2).norm();
  double alpha =
      180 / tools::conv::Pi * std::acos(box.col(1).dot(box.col(2)) / (b * c));
  double beta =
      180 / tools::conv::Pi * std::acos(box.col(0).dot(box.col(2)) / (a * c));
  double gamma =
      180 / tools::conv::Pi * std::acos(box.col(0).dot(box.col(1)) / (a * b));
  _out << boxfrmt % a % b % c % alpha % beta % gamma;
}

void PDBWriter::Close() {
  if (_out.is_open()) {
    _out.close();
  }
}

void PDBWriter::Write(CSG_Topology *conf) {
  if (_out.is_open()) {
    _out << boost::format("MODEL     %1$4d\n") % (conf->getStep() + 1)
         << std::flush;
    WriteContainer<CSG_Topology>(conf, *conf);
    _out << "ENDMDL" << std::endl;
  } else {
    throw runtime_error("Cannot write topology to file, file is not open.");
  }
}

}  // namespace csg
}  // namespace votca
