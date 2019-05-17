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

#pragma once
#ifndef VOTCA_CSG_CGENGINE_H
#define VOTCA_CSG_CGENGINE_H

#include "atomcgconverter.h"
#include "cgobserver.h"
#include "csgtopology.h"
#include <boost/program_options.hpp>
#include <list>
#include <map>
#include <votca/tools/datacollection.h>

#include "cgengine.h"
#include "molecule.h"
#include "nematicorder.h"
#include "topologyreader.h"
#include "trajectoryreader.h"
#include "trajectorywriter.h"
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

class CGEngine {
 public:
  /**
   * @brief Loads .xml files containing coarse graining stencils and mapping
   * information
   *
   * @param[in] filenames - a string of with .xml files, can be a single file or
   * multiple file names separated by ';'
   */
  void LoadFiles(std::string filenames);

  /**
   * @brief Takes the atomisitic topology and populates the coarse grained
   * topology with a coarse grained represenation
   *
   * @param[in] atomistic_top_in
   * @param[in,out] cg_top
   *
   * @return return an atom to cg converter which can be used to update the
   * positions vectors and forces of the coarse grained representation
   */
  std::unique_ptr<AtomCGConverter> PopulateCGTopology(
      CSG_Topology& atomistic_top_in, CSG_Topology& cg_top);

  /**
   * \brief Adds molecules that are to be ignored during the mapping process
   * \param molecule_type glob molecule_type for molecule molecule_type
   */
  void AddIgnore(std::string molecule_type) {
    ignores_.push_back(molecule_type);
  }

  /**
   * \brief checks whether molecule is ignored
   * \param ident identifyier of molecule
   * \return true if is ignored
   */
  bool IsIgnored(std::string molecule_type) const;

 private:
  std::unordered_set<std::string> file_names_;

  std::vector<std::string> ignores_;
};

inline bool CGEngine::IsIgnored(std::string molecule_type) const {
  for (const std::string& ignore : ignores_) {
    if (tools::wildcmp(ignore.c_str(), molecule_type.c_str())) return true;
  }
  return false;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGENGINE_H
