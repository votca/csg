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
#ifndef VOTCA_CSG_GMXTOPOLOGYREADER_PRIV_H
#define VOTCA_CSG_GMXTOPOLOGYREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double GMXTopologyReader<Topology_T>::formatDistance_(const double &distance) {
  return converter_.convert(distance_unit, Topology_T::distance_unit) *
         distance;
}

template <class Topology_T>
double GMXTopologyReader<Topology_T>::formatMass_(const double &mass) {
  return converter_.convert(mass_unit, Topology_T::mass_unit) * mass;
}

template <class Topology_T>
double GMXTopologyReader<Topology_T>::formatCharge_(const double &charge) {
  return converter_.convert(charge_unit, Topology_T::charge_unit) * charge;
}

template <class Topology_T>
std::string GMXTopologyReader<Topology_T>::formatElement_(
    const std::string &bead_type) {
  std::string element = tools::topology_constants::unassigned_element;
  if (elements_.isEleShort(bead_type)) {
    element = bead_type;
  }
  std::string name_all_caps = boost::to_upper_copy<std::string>(bead_type);
  if (elements_.isEleFull(name_all_caps)) {
    element = elements_.getEleShort(name_all_caps);
  }
  return element;
}

template <class Topology_T>
bool GMXTopologyReader<Topology_T>::ReadTopology(const std::string &file,
                                                 boost::any top_any) {

  if (typeid(Topology_T *) != top_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using gmx topology reader, incorrect "
        "topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(top_any);

  gmx_mtop_t mtop;

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
      params.set(tools::StructureParameter::MoleculeId,
                 static_cast<int>(top.MoleculeCount()));
      params.set(tools::StructureParameter::MoleculeType, molname);
      typename Topology_T::container_t &mi = top.CreateMolecule(params);

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

        std::string element = formatElement_(bead_type);
        tools::byte_t symmetry = 1;
        tools::StructureParameters params;
        params.set(tools::StructureParameter::Symmetry, symmetry);
        params.set(tools::StructureParameter::CSG_Mass, formatMass_(a->m));
        params.set(tools::StructureParameter::CSG_Charge, formatCharge_(a->q));
        params.set(tools::StructureParameter::Element, element);
        params.set(tools::StructureParameter::BeadId, a->atomnumber);
        params.set(tools::StructureParameter::BeadType, bead_type);
        params.set(tools::StructureParameter::ResidueId, a->resind);
        params.set(tools::StructureParameter::ResidueType, residue_name);
        params.set(tools::StructureParameter::MoleculeId, mi.getId());
        typename Topology_T::bead_t &bead = top.CreateBead(params);
        mi.AddBead(&bead);
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
      m(i, j) = formatDistance_(gbox[j][i]);
    }
  }
  top.setBox(m);

  return true;
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GMXTOPOLOGYREADER_PRIV_H
