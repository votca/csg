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
#include "../../include/votca/csg/interaction.h"
#include <stdexcept>

namespace votca {
namespace csg {

using namespace std;

string InteractionTypeToString(const InteractionType interaction_type) {
  if (interaction_type == InteractionType::unassigned) {
    return string("unassigned");
  } else if (interaction_type == InteractionType::bond) {
    return string("bond");
  } else if (interaction_type == InteractionType::angle) {
    return string("angle");
  } else if (interaction_type == InteractionType::dihedral) {
    return string("dihedral");
  }
  throw invalid_argument(
      "InteractionType is not recognized cannot convert to string.");
}

string Interaction::getLabel() const {
  stringstream label;
  label << "molecule id " << mol_id_;
  label << ":group name " << group_;
  label << ":group id " << group_id_;
  label << ":index " << index_;
  return label.str();
}

}  // namespace csg
}  // namespace votca
