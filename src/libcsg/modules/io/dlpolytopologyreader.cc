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

#include "../../../../include/votca/csg/csgtopology.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>  // IWYU pragma: keep
#include <boost/filesystem/convenience.hpp>      // IWYU pragma: keep
#include <boost/filesystem/path.hpp>             // IWYU pragma: keep
#include <fstream>
#include <stddef.h>
#include <votca/csg/molecule.h>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/types.h>

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "dlpolytopologyreader.h"

namespace votca {
namespace csg {
class Bead;
}  // namespace csg
}  // namespace votca

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

string DLPOLYTopologyReader::_NextKeyline(ifstream &fs, const char *wspace)
// function to find and read the next line starting with a keyword/directive
// (skipping comments starting with "#" or ";") NOTE: the line is returned
// case-preserved, not to alter molecule/atom names (consider if --no-map is
// used)
{
  string line;
  size_t i_nws = 0;

  do {
    getline(fs, line);

    if (fs.eof())
      throw std::runtime_error("Error: unexpected end of dlpoly file '" +
                               _fname + "'");

    i_nws = line.find_first_not_of(wspace);
  } while (line.substr(i_nws, 1) == "#" || line.substr(i_nws, 1) == ";");

  return line.substr(i_nws, line.size() - i_nws);
}

string DLPOLYTopologyReader::_NextKeyInt(ifstream &fs, const char *wspace,
                                         const string &word, int &ival)
// function to read the next line containing only a given keyword and an integer
// value after it (only skipping comments!) NOTE: this function must only be
// called when the next directive line has to contain the given keyword and an
// integer value
{
  stringstream sl(_NextKeyline(fs, wspace));
  string line, sval;

  sl >> line;  // allow user not to bother about the case
  boost::to_upper(line);

  if (line.substr(0, word.size()) != word)
    throw std::runtime_error("Error: unexpected line from dlpoly file '" +
                             _fname + "', expected '" + word + "' but got '" +
                             line + "'");

  sl >> sval;

  size_t i_num = sval.find_first_of(
      "0123456789");  // assume integer number straight after the only keyword

  if (i_num > 0)
    throw std::runtime_error("Error: missing integer number in directive '" +
                             line + "' in topology file '" + _fname + "'");

  ival = boost::lexical_cast<int>(sval);

  return sl.str();
}

bool DLPOLYTopologyReader::_isKeyInt(const string &line, const char *wspace,
                                     const string &word, int &ival)
// function to check if the given (last read) directive line starts with a given
// keyword and has an integer value at the end NOTE: comments are allowed beyond
// the integer (anything after the first integer is ignored)
{
  // split directives consisting of a few words: the sought keyword must be the
  // first one!
  Tokenizer tok(line, wspace);
  vector<string> fields;
  tok.ToVector(fields);

  ival = 0;

  if (fields.size() < 2) return false;

  boost::to_upper(fields[0]);

  if (fields[0].substr(0, word.size()) != word)
    throw std::runtime_error("Error: unexpected directive from dlpoly file '" +
                             _fname + "', expected keyword '" + word +
                             "' but got '" + fields[0] + "'");

  size_t i_num = string::npos;

  unsigned int i = 1;
  do {  // find integer number in the field with the lowest index (closest to
        // the keyword)
    i_num = fields[i++].find_first_of("0123456789");
  } while (i_num > 0 && i < fields.size());

  if (i_num > 0) return false;

  ival = boost::lexical_cast<int>(fields[i - 1]);

  return true;
}

bool DLPOLYTopologyReader::ReadTopology(string file, CSG_Topology &top) {
  const char *WhiteSpace = " \t";

  int matoms = 0;
  int natoms = 0;

  std::ifstream fl;
  boost::filesystem::path filepath(file.c_str());

  string line;

  if (boost::filesystem::basename(filepath).size() == 0) {
    if (filepath.parent_path().string().size() == 0) {
      _fname = "FIELD";  // DL_POLY uses fixed file names in current/working
                         // directory
    } else {
      _fname = filepath.parent_path().string() + "/FIELD";
    }
  } else {
    _fname = file;
  }

  fl.open(_fname.c_str());

  if (!fl.is_open()) {
    throw std::runtime_error("Error on opening dlpoly file '" + _fname + "'");
  } else {

    line = _NextKeyline(fl, WhiteSpace);  // read title line and skip it
    line = _NextKeyline(fl, WhiteSpace);  // read next directive line
    boost::to_upper(line);

    if (line.substr(0, 4) == "UNIT") {      // skip 'unit' line
      line = _NextKeyline(fl, WhiteSpace);  // read next directive line
      boost::to_upper(line);
    }

    if (line.substr(0, 4) == "NEUT") {  // skip 'neutral groups' line (DL_POLY
                                        // Classic FIELD format)
      line = _NextKeyline(fl, WhiteSpace);  // look for next directive line
      boost::to_upper(line);
    }

    int nmol_types;

    if (!_isKeyInt(line, WhiteSpace, "MOLEC", nmol_types))
      throw std::runtime_error("Error: missing integer number in directive '" +
                               line + "' in topology file '" + _fname + "'");

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : '" << line << "' - "
         << nmol_types << endl;
#endif

    string mol_name;

    for (int nmol_type = 0; nmol_type < nmol_types; nmol_type++) {

      mol_name = _NextKeyline(fl, WhiteSpace);
      Molecule *mi = top.CreateMolecule(top.MoleculeCount(), mol_name);

      int nreplica = 1;
      line = _NextKeyInt(fl, WhiteSpace, "NUMMOL", nreplica);

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << mol_name
           << "' - '" << line << "' - " << nreplica << endl;
#endif

      line = _NextKeyInt(fl, WhiteSpace, "ATOMS", natoms);

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << line << "' - "
           << natoms << endl;
#endif

      // read molecule
      int id_map[natoms];
      for (int i = 0; i < natoms;) {  // i is altered in repeater loop
        stringstream sl(_NextKeyline(fl, WhiteSpace));

#ifdef DEBUG
        cout << "Read atom specs from dlpoly topology : '" << sl.str() << "'"
             << endl;
#endif
        string beadtype;
        sl >> beadtype;

        Elements elements;
        string element;
        if (elements.isEleShort(beadtype)) {
          element = beadtype;
        }
        string full_name = boost::to_upper_copy<string>(beadtype);
        if (elements.isEleFull(full_name)) {
          element = full_name;
        }

        double mass;
        sl >> mass;
        double charge;
        sl >> charge;

        line = " ";
        sl >> line;  // rest of the atom line

        Tokenizer tok(line, WhiteSpace);
        vector<string> fields;
        tok.ToVector(fields);

#ifdef DEBUG
        cout << "Rest atom specs from dlpoly topology : '" << line << "'"
             << endl;
#endif

        int repeater = 1;
        if (fields.size() > 1) repeater = boost::lexical_cast<int>(fields[0]);

        for (int j = 0; j < repeater; j++) {

          // string beadname = beadtype + "#" + boost::lexical_cast<string>(i +
          // 1);
          byte_t symmetry = 1;

          Bead *bead =
              top.CreateBead(symmetry, beadtype, top.BeadCount(), mi->getId(),
                             topology_constants::unassigned_residue_id,
                             topology_constants::unassigned_residue_type,
                             element, mass, charge);

          mi->AddBead(bead);
          id_map[i] = bead->getId();
          i++;
#ifdef DEBUG
          cout << "Atom identification in maps '" << bead->getLabel() << "'"
               << endl;
#endif
        }
        matoms += repeater;
      }

      while (line != "FINISH") {
        stringstream nl(_NextKeyline(fl, WhiteSpace));
        nl >> line;
#ifdef DEBUG
        cout << "Read unit type# from dlpoly topology : '" << nl.str() << "'"
             << endl;
#endif

        boost::to_upper(line);
        line = line.substr(0, 6);
        if ((line == "BONDS") || (line == "ANGLES") || (line == "DIHEDR")) {
          string interaction_group = line;
          int count;
          nl >> count;
          for (int interaction_id = 0; interaction_id < count;
               interaction_id++) {

            stringstream sl(_NextKeyline(fl, WhiteSpace));
#ifdef DEBUG
            cout << "Read unit specs from dlpoly topology : '" << sl.str()
                 << "'" << endl;
#endif
            sl >> line;  // internal dlpoly bond/angle/dihedral function types
                         // are merely skipped (ignored)
            int ids[4];
            Interaction *ic = NULL;
            sl >> ids[0];
            sl >> ids[1];
            if (interaction_group == "BONDS") {
              // ic = new IBond(id_map[ids[0] - 1],
              //               id_map[ids[1] - 1]);  // -1 due to fortran vs c
              int bead_id1 = id_map[ids[0] - 1];
              int bead_id2 = id_map[ids[1] - 1];
              ic = top.CreateInteraction(
                  Interaction::interaction_type::bond, interaction_group,
                  interaction_id, mi->getId(), vector<int>{bead_id1, bead_id2});
            } else if (interaction_group == "ANGLES") {
              sl >> ids[2];
              // ic = new IAngle(id_map[ids[0] - 1], id_map[ids[1] - 1],
              //                id_map[ids[2] - 1]);  // -1 due to fortran vs c
              int bead_id1 = id_map[ids[0] - 1];
              int bead_id2 = id_map[ids[1] - 1];
              int bead_id3 = id_map[ids[2] - 1];
              ic = top.CreateInteraction(
                  Interaction::interaction_type::angle, interaction_group,
                  interaction_id, mi->getId(),
                  vector<int>{bead_id1, bead_id2, bead_id3});
            } else if (interaction_group.substr(0, 6) == "DIHEDR") {
              interaction_group = "DIHEDRALS";
              sl >> ids[2];
              sl >> ids[3];
              // ic = new IDihedral(id_map[ids[0] - 1], id_map[ids[1] - 1],
              //                   id_map[ids[2] - 1],
              //                 id_map[ids[3] - 1]);  // -1 due to fortran vs c
              int bead_id1 = id_map[ids[0] - 1];
              int bead_id2 = id_map[ids[1] - 1];
              int bead_id3 = id_map[ids[2] - 1];
              int bead_id4 = id_map[ids[3] - 1];
              ic = top.CreateInteraction(
                  Interaction::interaction_type::dihedral, interaction_group,
                  interaction_id, mi->getId(),
                  vector<int>{bead_id1, bead_id2, bead_id3, bead_id4});
            } else {
              throw std::runtime_error(
                  "Error: interaction_group should be BONDS, ANGLES or "
                  "DIHEDRALS");
            }
            // could one use bond/angle/dihedral function types for 1:1 mapping?
            // (CG map overwrites ic->Group anyway)
            // ic->setGroup(line);
            // ic->setGroup(interaction_group);
            // ic->setIndex(interaction_id);
            // ic->setMoleculeId(mi->getId());
            mi->AddInteraction(ic);
          }
        }
      }

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << line
           << "' - done with '" << mol_name << "'" << endl;
#endif

      // replicate molecule
      for (int replica = 1; replica < nreplica; replica++) {
        Molecule *mi_replica =
            top.CreateMolecule(top.MoleculeCount(), mol_name);
        vector<int> bead_ids = mi->getBeadIds();
        for (const int &bead_id : bead_ids) {
          Bead *bead = mi->getBead(bead_id);
          byte_t symmetry = 1;
          Bead *bead_replica = top.CreateBead(
              symmetry, bead->getType(), top.BeadCount(), mi_replica->getId(),
              bead->getResidueId(), bead->getResidueType(), bead->getElement(),
              bead->getMass(), bead->getQ());
          mi_replica->AddBead(bead_replica);
        }
        matoms += mi->BeadCount();
        vector<Interaction *> ics = mi->Interactions();
        for (vector<Interaction *>::iterator ic = ics.begin(); ic != ics.end();
             ++ic) {
          Interaction *ic_replica = NULL;
          int offset =
              mi_replica->getBead(0)->getId() - mi->getBead(0)->getId();
          if ((*ic)->BeadCount() == 2) {
            // ic_replica = new IBond((*ic)->getBeadId(0) + offset,
            //                       (*ic)->getBeadId(1) + offset);
            int bead_id1 = (*ic)->getBeadId(0) + offset;
            int bead_id2 = (*ic)->getBeadId(1) + offset;
            ic_replica = top.CreateInteraction(
                Interaction::interaction_type::bond, (*ic)->getGroup(),
                (*ic)->getIndex(), mi_replica->getId(),
                vector<int>{bead_id1, bead_id2});
          } else if ((*ic)->BeadCount() == 3) {
            // ic_replica = new IAngle((*ic)->getBeadId(0) + offset,
            //                        (*ic)->getBeadId(1) + offset,
            //                        (*ic)->getBeadId(2) + offset);
            int bead_id1 = (*ic)->getBeadId(0) + offset;
            int bead_id2 = (*ic)->getBeadId(1) + offset;
            int bead_id3 = (*ic)->getBeadId(2) + offset;
            ic_replica = top.CreateInteraction(
                Interaction::interaction_type::angle, (*ic)->getGroup(),
                (*ic)->getIndex(), mi_replica->getId(),
                vector<int>{bead_id1, bead_id2, bead_id3});
          } else if ((*ic)->BeadCount() == 4) {
            // ic_replica = new IDihedral(
            //    (*ic)->getBeadId(0) + offset, (*ic)->getBeadId(1) + offset,
            //    (*ic)->getBeadId(2) + offset, (*ic)->getBeadId(3) + offset);
            int bead_id1 = (*ic)->getBeadId(0) + offset;
            int bead_id2 = (*ic)->getBeadId(1) + offset;
            int bead_id3 = (*ic)->getBeadId(2) + offset;
            int bead_id4 = (*ic)->getBeadId(3) + offset;
            ic_replica = top.CreateInteraction(
                Interaction::interaction_type::dihedral, (*ic)->getGroup(),
                (*ic)->getIndex(), mi_replica->getId(),
                vector<int>{bead_id1, bead_id2, bead_id3, bead_id4});
          } else {
            throw std::runtime_error("Error: BeadCount not equal 2, 3 or 4");
          }
          // ic_replica->setGroup((*ic)->getGroup());
          // ic_replica->setIndex((*ic)->getIndex());
          // ic_replica->setMoleculeId(mi_replica->getId());
          mi_replica->AddInteraction(ic_replica);
        }
      }
    }
    top.RebuildExclusions();
  }

#ifdef DEBUG
  getline(fl, line);  // is "close" found?
  if (line == "close") {
    cout << "Read from dlpoly file '" << _fname << "' : '" << line
         << "' - done with topology" << endl;
  } else {
    cout << "Read from dlpoly file '" << _fname
         << "' : 'EOF' - done with topology (directive 'close' not read!)"
         << endl;
  }
#endif

  // we don't need the rest
  fl.close();

  return true;
}

}  // namespace csg
}  // namespace votca
