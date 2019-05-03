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

#ifndef VOTCA_CSG_GMXTOPOLOGYREADER_H
#define VOTCA_CSG_GMXTOPOLOGYREADER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "../topologyreader.h"
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/structureparameters.h>
#include <votca/tools/types.h>

#include <gromacs/fileio/tpxio.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/topology/atoms.h>
#include <gromacs/topology/topology.h>

// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

/**
    \brief reader for gromacs topology files

    This class encapsulates the gromacs reading functions and provides an
   interface to fill a topolgy class

*/
template <class Topology_T>
class GMXTopologyReader : public TopologyReader {
 public:
  GMXTopologyReader() {}

  /// read a topology file
  bool ReadTopology(std::string file, boost::any top);

 private:
};

template <class Topology_T>
bool GMXTopologyReader<Topology_T>::ReadTopology(std::string file,
                                                 boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gmx topology reader, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);

  gmx_mtop_t mtop;

  tools::Elements elements;

  int natoms;
  // cleanup topology to store new data
  top.Cleanup();

  t_inputrec ir;
  ::matrix gbox;

  (void)read_tpx((char *)file.c_str(), &ir, gbox, &natoms, NULL, NULL, &mtop);

  size_t ifirstatom = 0;

#if GROMACS_VERSION >= 20190000
  size_t nmolblock = mtop.molblock.size();
#else
  size_t nmolblock = mtop.nmolblock;
#endif

  for (size_t iblock = 0; iblock < nmolblock; ++iblock) {
    gmx_moltype_t *mol = &(mtop.moltype[mtop.molblock[iblock].type]);

    std::string molname = *(mol->name);

    // This is to check that you are not adding another residue with the same id
    // as one that was previously added
    // int res_offset = top->getMaxResidueId() + 1;

    t_atoms *atoms = &(mol->atoms);

    for (int imol = 0; imol < mtop.molblock[iblock].nmol; ++imol) {
      tools::StructureParameters params;
      params.set(tools::StructureParameter::MoleculeId, top.MoleculeCount());
      params.set(tools::StructureParameter::MoleculeType, molname);
      typename Topology_T::molecule_t *mi = top.CreateMolecule(params);

#if GROMACS_VERSION >= 20190000
      size_t natoms_mol = mtop.moltype[mtop.molblock[iblock].type].atoms.nr;
#else
      size_t natoms_mol = mtop.molblock[iblock].natoms_mol;
#endif
      // read the atoms
      for (size_t iatom = 0; iatom < natoms_mol; iatom++) {
        t_atom *a = &(atoms->atom[iatom]);
        std::string residue_name = *(atoms->resinfo[iatom].name);

        std::string bead_type = *(atoms->atomtype[iatom]);

        std::string element = tools::topology_constants::unassigned_element;
        if (elements.isEleShort(bead_type)) {
          element = bead_type;
        }
        std::string name_all_caps =
            boost::to_upper_copy<std::string>(bead_type);
        if (elements.isEleFull(name_all_caps)) {
          element = elements.getEleShort(name_all_caps);
        }

        tools::byte_t symmetry = 1;
        tools::StructureParameters params;
        params.set(tools::StructureParameter::Symmetry, symmetry);
        params.set(tools::StructureParameter::Mass, a->m);
        params.set(tools::StructureParameter::Charge, a->q);
        params.set(tools::StructureParameter::Element, element);
        params.set(tools::StructureParameter::BeadId, a->atomnumber);
        params.set(tools::StructureParameter::BeadType, bead_type);
        params.set(tools::StructureParameter::ResidueId, a->resind);
        params.set(tools::StructureParameter::ResidueType, residue_name);
        params.set(tools::StructureParameter::MoleculeId, mi->getId());
        typename Topology_T::bead_t *bead = top.CreateBead(params);
        mi->AddBead(bead);
      }

      // add exclusions
      for (size_t iatom = 0; iatom < natoms_mol; iatom++) {
        // read exclusions
        t_blocka *excl = &(mol->excls);
        // insert exclusions
        std::list<typename Topology_T::bead_t *> excl_list;
        for (int k = excl->index[iatom]; k < excl->index[iatom + 1]; k++) {
          excl_list.push_back(top.getBead(excl->a[k] + ifirstatom));
        }
        top.InsertExclusion(top.getBead(iatom + ifirstatom), excl_list);
      }
      ifirstatom += natoms_mol;
    }
  }

  Eigen::Matrix3d m;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m(i, j) = gbox[j][i];
    }
  }
  top.setBox(m);

  return true;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GMXTOPOLOGYREADER_H
