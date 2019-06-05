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
#ifndef VOTCA_CSG_DLPOLYTRAJECTORYREADER_PRIV_H
#define VOTCA_CSG_DLPOLYTRAJECTORYREADER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double DLPOLYTrajectoryReader<Topology_T>::formatDistance_(
    const double &distance) const {
  return converter_.convert(distance_unit, Topology_T::units::distance_unit) *
         distance;
}

template <class Topology_T>
Eigen::Vector3d DLPOLYTrajectoryReader<Topology_T>::formatDistance_(
    const Eigen::Vector3d &distance) const {
  return converter_.convert(distance_unit, Topology_T::units::distance_unit) *
         distance;
}

template <class Topology_T>
Eigen::Vector3d DLPOLYTrajectoryReader<Topology_T>::formatVelocity_(
    const Eigen::Vector3d &velocity) const {
  return converter_.convert(velocity_unit, Topology_T::units::velocity_unit) *
         velocity;
}

template <class Topology_T>
Eigen::Vector3d DLPOLYTrajectoryReader<Topology_T>::formatForce_(
    const Eigen::Vector3d &force) const {
  return converter_.convert(force_unit, Topology_T::units::force_unit) * force;
}

template <class Topology_T>
double DLPOLYTrajectoryReader<Topology_T>::formatTime_(
    const double &time) const {
  return converter_.convert(time_unit, Topology_T::units::time_unit) * time;
}

template <class Topology_T>
bool DLPOLYTrajectoryReader<Topology_T>::Open(const std::string &file)
// open the original dlpoly configuration or trajectory file
// NOTE: allowed file naming - <name>.dlpc or <name>.dlph (convention:
// ".dlpc"="CONFIG", ".dlph"="HISTORY")
{
  boost::filesystem::path filepath(file.c_str());
  std::string inp_name = "HISTORY";

  if (boost::filesystem::extension(filepath).size() == 0) {

    throw std::ios_base::failure(
        "Error on opening dlpoly file '" + file +
        "' - extension is expected, use .dlph or .dlpc");

  } else if (boost::filesystem::extension(filepath) == ".dlpc") {

    _isConfig = true;
    inp_name = "CONFIG";

  } else if (boost::filesystem::extension(filepath) == ".dlph") {

    _isConfig = false;

  } else {
    throw std::ios_base::failure("Error on opening dlpoly file '" + file +
                                 "' - wrong extension, use .dlph or .dlpc");
  }

  if (boost::filesystem::basename(filepath).size() == 0) {
    if (filepath.parent_path().string().size() == 0) {
      _fname = inp_name;
    } else {
      _fname = filepath.parent_path().string() + "/" + inp_name;
    }
  } else {
    _fname = file;
  }

  _fl.open(_fname.c_str());
  if (!_fl.is_open())
    throw std::ios_base::failure("Error on opening dlpoly file '" + _fname +
                                 "'");
  return true;
}

template <class Topology_T>
void DLPOLYTrajectoryReader<Topology_T>::Close() {
  _fl.close();
}

template <class Topology_T>
bool DLPOLYTrajectoryReader<Topology_T>::FirstFrame(boost::any &conf) {
  _first_frame = true;
  bool res = NextFrame(conf);
  _first_frame = false;
  return res;
}

template <class Topology_T>
bool DLPOLYTrajectoryReader<Topology_T>::NextFrame(boost::any &conf_any) {

  if (typeid(Topology_T) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using dlpoly trajectory reader, incorrect "
        "topology type provided.");
  }

  Topology_T &conf = boost::any_cast<Topology_T>(conf_any);
  static bool hasVs = false;
  static bool hasFs = false;
  static int mavecs =
      0;  // number of 3d vectors per atom = keytrj in DL_POLY manuals
  static int mpbct = 0;   // cell PBC type = imcon in DL_POLY manuals
  static int matoms = 0;  // number of atoms/beads in a frame
  // const double scale = tools::conv::ang2nm;

  static int nerrt = 0;

  std::string line;

  BoundaryCondition::eBoxtype pbc_type = BoundaryCondition::typeAuto;

  if (_first_frame) {

    getline(_fl, line);  // title

#ifdef DEBUG
    std::cout << "Read from dlpoly file '" << _fname << "' : '" << line
              << "' - header" << std::endl;
#endif

    getline(_fl, line);  // 2nd header line

#ifdef DEBUG
    std::cout << "Read from dlpoly file '" << _fname << "' : '" << line
              << "' - directives line" << std::endl;
#endif

    tools::Tokenizer tok(line, " \t");
    std::vector<std::string> fields;
    tok.ToVector(fields);

    if (fields.size() < 3)
      throw std::runtime_error("Error: too few directive switches (<3) in '" +
                               _fname + "' header (check its 2-nd line)");

    mavecs = boost::lexical_cast<int>(fields[0]);
    mpbct = boost::lexical_cast<int>(fields[1]);
    matoms = boost::lexical_cast<int>(fields[2]);

    hasVs = (mavecs > 0);  // 1 or 2 => in DL_POLY frame velocity std::vector
                           // follows coords for each atom/bead
    hasFs = (mavecs > 1);  // 2      => in DL_POLY frame force std::vector
                           // follows velocities for each atom/bead

#ifdef DEBUG
    if (hasVs != conf.HasVel() || hasFs != conf.HasForce()) {
      std::cout << "WARNING: N of atom vectors (keytrj) in '" << _fname
                << "' header differs from that read with topology" << std::endl;
    }
#endif

    conf.SetHasVel(hasVs);
    conf.SetHasForce(hasFs);

#ifdef DEBUG
    std::cout << "Read from dlpoly file '" << _fname << "' : keytrj - "
              << mavecs << ", hasV - " << conf.HasVel() << ", hasF - "
              << conf.HasForce() << std::endl;
#endif

    if (static_cast<size_t>(matoms) != conf.BeadCount())
      throw std::runtime_error("Number of atoms/beads in '" + _fname +
                               "' header differs from that read with topology");

    if (mpbct == 0) {
      pbc_type = BoundaryCondition::typeOpen;
    } else if (mpbct == 1 || mpbct == 2) {
      pbc_type = BoundaryCondition::typeOrthorhombic;
    } else if (mpbct == 3) {
      pbc_type = BoundaryCondition::typeTriclinic;
    }

#ifdef DEBUG
    std::cout << "Read from dlpoly file '" << _fname
              << "' : pbc_type (imcon) - '" << pbc_type << "'" << std::endl;

    if (pbc_type != conf.getBoxType())
      std::cout << "WARNING: PBC type in dlpoly file '" << _fname
                << "' header differs from that read with topology" << std::endl;
// throw std::runtime_error("Error: Boundary conditions in '"+_fname+"'
// header differs from that read with topology");
#endif
  } else if (_isConfig) {

    return false;
  }
  // read normal frame

  if (!_isConfig) {
    getline(_fl, line);  // timestep line - only present in HISTORY, and not
                         // in CONFIG
#ifdef DEBUG
    std::cout << "Read from dlpoly file '" << _fname << "' : '" << line << "'"
              << std::endl;
#endif
  }

  if (!_fl.eof()) {
    double dtime, stime;
    int nstep;
    int natoms;
    int navecs;
    int npbct;

    if (_isConfig) {
      // use the above read specs from the header, and skip the data missing
      // in CONFIG

      natoms = matoms;
      navecs = mavecs;
      npbct = mpbct;

      conf.SetHasVel(hasVs);
      conf.SetHasForce(hasFs);

#ifdef DEBUG
      std::cout << "Read from CONFIG: traj_key - " << navecs << ", hasV - "
                << conf.HasVel() << ", hasF - " << conf.HasForce() << std::endl;
#endif

    } else {

      tools::Tokenizer tok(line, " \t");
      std::vector<std::string> fields;
      tok.ToVector(fields);

      if (fields.size() < 6)
        throw std::runtime_error(
            "Error: too few directive switches (<6) in 'timestep' record");

      nstep = boost::lexical_cast<int>(fields[1]);
      natoms = boost::lexical_cast<int>(fields[2]);
      navecs = boost::lexical_cast<int>(fields[3]);
      npbct = boost::lexical_cast<int>(fields[4]);
      dtime = stod(fields[5]);  // normally it is the 5-th column in
                                // 'timestep' line
      stime = stod(fields[fields.size() - 1]);  // normally it is the last
                                                // column in 'timestep' line

#ifdef DEBUG
      std::cout << "Read from dlpoly file '" << _fname
                << "' : natoms = " << natoms << ", levcfg = " << fields[3];
      std::cout << ", dt = " << fields[5] << ", time = " << stime << std::endl;
#endif

      if (static_cast<size_t>(natoms) != conf.BeadCount())
        throw std::runtime_error(
            "Error: N of atoms/beads in '" + _fname +
            "' header differs from that found in topology");
      if (natoms != matoms)
        throw std::runtime_error(
            "Error: N of atoms/beads in '" + _fname +
            "' header differs from that found in the frame");
      if (navecs != mavecs)
        throw std::runtime_error(
            "Error: N of atom vectors (keytrj) in '" + _fname +
            "' header differs from that found in the frame");
      if (npbct != mpbct)
        throw std::runtime_error(
            "Error: boundary conditions (imcon) in '" + _fname +
            "' header differs from that found in the frame");

      // total time - calculated as product due to differences between DL_POLY
      // versions in HISTORY formats
      conf.setTime(nstep * formatTime_(dtime));
      conf.setStep(nstep);

      if (std::abs(stime - conf.getTime()) > 1.e-8) {
        nerrt++;
        if (nerrt < 11) {
          std::cout << "Check: nstep = " << nstep << ", dt = " << dtime
                    << ", time = " << stime << " (correct?)" << std::endl;
          // std::cout << "Check: nstep = " << nstep << ", dt = " << dtime << ",
          // time = " << conf.getTime() << " (correct?)" << std::endl;
        } else if (nerrt == 11) {
          std::cout
              << "Check: timestep - more than 10 mismatches in total time "
                 "found..."
              << std::endl;
        }
      }

      if (npbct == 0) {
        pbc_type = BoundaryCondition::typeOpen;
      } else if (npbct == 1 || npbct == 2) {
        pbc_type = BoundaryCondition::typeOrthorhombic;
      } else if (npbct == 3) {
        pbc_type = BoundaryCondition::typeTriclinic;
      }
    }

    Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; i++) {  // read 3 box/cell lines

      getline(_fl, line);

#ifdef DEBUG
      std::cout << "Read from dlpoly file '" << _fname << "' : '" << line
                << "' - box vector # " << i + 1 << std::endl;
#endif

      if (_fl.eof())
        throw std::runtime_error("Error: unexpected EOF in dlpoly file '" +
                                 _fname + "', when reading box vector" +
                                 boost::lexical_cast<std::string>(i));

      tools::Tokenizer tok(line, " \t");
      std::vector<double> fields;
      tok.ConvertToVector<double>(fields);
      // Angs -> nm
      box.col(i) = Eigen::Vector3d(formatDistance_(fields[0]),
                                   formatDistance_(fields[1]),
                                   formatDistance_(fields[2]));
    }

    conf.setBox(box, pbc_type);

    for (int i = 0; i < natoms; i++) {

      {
        getline(_fl, line);  // atom header line

#ifdef DEBUG
        std::cout << "Read from dlpoly file '" << _fname << "' : '" << line
                  << "'" << std::endl;
#endif

        if (_fl.eof())
          throw std::runtime_error("Error: unexpected EOF in dlpoly file '" +
                                   _fname + "', when reading atom/bead # " +
                                   boost::lexical_cast<std::string>(i + 1));

        std::vector<std::string> fields;
        tools::Tokenizer tok(line, " \t");
        tok.ToVector(fields);
        int id = boost::lexical_cast<int>(fields[1]);
        if (i + 1 != id)
          throw std::runtime_error(
              "Error: unexpected atom/bead index in dlpoly file '" + _fname +
              "' : expected " + boost::lexical_cast<std::string>(i + 1) +
              " but got " + boost::lexical_cast<std::string>(id));
      }

      typename Topology_T::bead_t &b = conf.getBead(i);
      Eigen::Matrix3d atom_vecs = Eigen::Matrix3d::Zero();
      for (int j = 0; j < std::min(navecs, 2) + 1; j++) {

        getline(_fl, line);  // read atom positions

#ifdef DEBUG
        std::cout << "Read from dlpoly file '" << _fname << "' : '" << line
                  << "'" << std::endl;
#endif

        if (_fl.eof())
          throw std::runtime_error(
              "Error: unexpected EOF in dlpoly file '" + _fname +
              "', when reading atom/bead vector # " +
              boost::lexical_cast<std::string>(j) + " of atom " +
              boost::lexical_cast<std::string>(i + 1));

        std::vector<double> fields;
        tools::Tokenizer tok(line, " \t");
        tok.ConvertToVector<double>(fields);
        // Angs -> nm
        atom_vecs.col(j) = Eigen::Vector3d(fields[0], fields[1], fields[2]);
      }

      b.setPos(formatDistance_(atom_vecs.col(0)));
#ifdef DEBUG
      std::cout << "Crds from dlpoly file '" << _fname
                << "' : " << atom_vecs.col(0) << std::endl;
#endif
      if (navecs > 0) {
        b.setVel(formatVelocity_(atom_vecs.col(1)));
#ifdef DEBUG
        std::cout << "Vels from dlpoly file '" << _fname
                  << "' : " << atom_vecs.col(1) << std::endl;
#endif
        if (navecs > 1) {
          b.setF(formatForce_(atom_vecs.col(2)));
#ifdef DEBUG
          std::cout << "Frcs from dlpoly file '" << _fname
                    << "' : " << atom_vecs.col(2) << std::endl;
#endif
        }
      }
    }
  }

  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_DLPOLYTRAJECTORYREADER_PRIV_H