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

#include "bondedstatistics.h"
#include "../../include/votca/csg/interaction.h"

using namespace votca::tools;

namespace votca {
namespace csg {

void BondedStatistics::BeginCG(CSG_Topology *top, CSG_Topology *top_atom) {
  const vector<unique_ptr<Interaction>> &interactions =
      top->BondedInteractions();
  vector<unique_ptr<Interaction>>::const_iterator ia;

  _bonded_values.clear();
  for (ia = interactions.begin(); ia != interactions.end(); ++ia) {
    _bonded_values.CreateArray((*ia)->getLabel());
  }
}

void BondedStatistics::EndCG() {}

void BondedStatistics::EvalConfiguration(CSG_Topology *conf,
                                         CSG_Topology *conv_atom) {
  const vector<unique_ptr<Interaction>> &interactions =
      conf->BondedInteractions();
  vector<unique_ptr<Interaction>>::const_iterator ia;

  DataCollection<double>::container::iterator is;
  for (ia = interactions.begin(), is = _bonded_values.begin();
       ia != interactions.end(); ++ia, ++is) {
    vector<int> bead_ids = (*ia)->getBeadIds();
    unordered_map<int, const Eigen::Vector3d *> bead_positions =
        conf->getBeadPositions(bead_ids);
    double value =
        (*ia)->EvaluateVar(*(conf->getBoundaryCondition()), bead_positions);
    cout << value << endl;
    (*is)->push_back(value);
  }
}

}  // namespace csg
}  // namespace votca
