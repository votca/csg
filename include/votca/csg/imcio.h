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

#ifndef VOTCA_CSG_IMCIO_H
#define VOTCA_CSG_IMCIO_H

#include <list>
#include <string>
#include <votca/tools/eigen.h>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;
void imcio_write_dS(const std::string &file, Eigen::VectorXd &r,
                    Eigen::VectorXd &dS, std::list<int> *list = NULL);
void imcio_write_matrix(const std::string &file, Eigen::MatrixXd &gmc,
                        std::list<int> *list = NULL);
void imcio_write_index(const std::string &file, std::vector<std::string> &names,
                       std::vector<TOOLS::RangeParser> &ranges);

void imcio_read_dS(const std::string &file, Eigen::VectorXd &r,
                   Eigen::VectorXd &dS);
void imcio_read_matrix(const std::string &file, Eigen::MatrixXd &gmc);
void imcio_read_index(const std::string &file, std::vector<std::string> &names,
                      std::vector<TOOLS::RangeParser> &ranges);

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_IMCIO_H
