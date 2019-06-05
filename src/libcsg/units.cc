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

#include "../../include/votca/csg/units.h"
namespace votca {
namespace csg {

const tools::DistanceUnit Units::distance_unit =
    tools::DistanceUnit::nanometers;
const tools::MassUnit Units::mass_unit = tools::MassUnit::atomic_mass_units;
const tools::TimeUnit Units::time_unit = tools::TimeUnit::picoseconds;
const tools::ChargeUnit Units::charge_unit = tools::ChargeUnit::e;
const tools::EnergyUnit Units::energy_unit =
    tools::EnergyUnit::kilojoules_per_mole;
const tools::VelocityUnit Units::velocity_unit =
    tools::VelocityUnit::nanometers_per_picosecond;
const tools::ForceUnit Units::force_unit =
    tools::ForceUnit::kilojoules_per_mole_nanometer;

}  // namespace csg
}  // namespace votca