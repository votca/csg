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

#ifndef SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
#define SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_

#include "../topologyreader.h"
#include "../trajectoryreader.h"
#include "hdf5.h"
#include <boost/any.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <votca/tools/matrix.h>

namespace votca {  // NOLINT
namespace csg {
/**
    \brief class for reading H5MD trajectory.

    This class implements the H5MD trajectory reading function. The format of
   the H5MD file is defined in Pierre de Buyl, Peter H. Colberg, Felix Höfling,
   H5MD: A structured, efficient, and portable file format for molecular data,
   http://dx.doi.org/10.1016/j.cpc.2014.01.018 The current reference is
   available here: http://nongnu.org/h5md/
*/
template <class Bead_T, class Molecule_T, class Topology_T>
class H5MDTrajectoryReader : public TrajectoryReader {
 public:
  H5MDTrajectoryReader();
  ~H5MDTrajectoryReader();

  /// Opens original trajectory file.
  bool Open(const std::string &file);

  /// Initialize data structures.
  void Initialize(Topology_T &top);

  /// Reads in the first frame.
  bool FirstFrame(boost::any conf);  // NOLINT

  /// Reads in the next frame.
  bool NextFrame(boost::any conf);  // NOLINT

  /// Closes original trajectory file.
  void Close();

 private:
  enum DatasetState { NONE, STATIC, TIMEDEPENDENT };

  /// Reads dataset that contains vectors.
  template <typename T1>
  T1 *ReadVectorData(hid_t ds, hid_t ds_data_type, int row) {
    hsize_t offset[3];
    offset[0] = row;
    offset[1] = 0;
    offset[2] = 0;
    hsize_t chunk_rows[3];
    chunk_rows[0] = 1;
    chunk_rows[1] = N_particles_;
    chunk_rows[2] = vec_components_;
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, NULL, chunk_rows, NULL);
    hid_t mspace1 = H5Screate_simple(vec_components_, chunk_rows, NULL);
    T1 *data_out = new T1[N_particles_ * vec_components_];
    herr_t status =
        H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status < 0) {
      throw std::runtime_error("Error ReadVectorData: " +
                               boost::lexical_cast<std::string>(status));
    } else
      return data_out;
  }

  /// Reads dataset with scalar values.
  template <typename T1>
  T1 *ReadScalarData(hid_t ds, hid_t ds_data_type, int row) {
    hsize_t offset[2];
    offset[0] = row;
    offset[1] = 0;
    hsize_t ch_rows[2];
    ch_rows[0] = 1;
    ch_rows[1] = N_particles_;
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, NULL, ch_rows, NULL);
    hid_t mspace1 = H5Screate_simple(2, ch_rows, NULL);
    T1 *data_out = new T1[N_particles_];
    herr_t status =
        H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status < 0) {
      throw std::runtime_error("Error ReadScalarData: " +
                               boost::lexical_cast<std::string>(status));
    } else {
      return data_out;
    }
  }

  template <typename T1>
  void ReadStaticData(hid_t ds, hid_t ds_data_type, T1 &outbuf) {
    herr_t status =
        H5Dread(ds, ds_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, outbuf);
    if (status < 0) {
      H5Eprint(H5E_DEFAULT, stderr);
    }
  }

  void CheckError(hid_t hid, std::string error_message) {
    if (hid < 0) {
      // H5Eprint(H5E_DEFAULT, stderr);
      throw std::runtime_error(error_message);
    }
  }

  bool GroupExists(hid_t file_id, std::string path) {
    H5G_stat_t info;
    herr_t status = H5Gget_objinfo(file_id, path.c_str(), 0, &info);
    if (status < 0) return false;
    return info.type == H5G_GROUP;
  }

  hid_t file_id_;
  hid_t ds_atom_position_;
  hid_t ds_atom_force_;
  hid_t ds_atom_velocity_;
  hid_t ds_atom_id_;
  hid_t ds_edges_group_;

  hid_t particle_group_;
  hid_t atom_position_group_;
  hid_t atom_force_group_;
  hid_t atom_velocity_group_;
  hid_t atom_id_group_;
  hid_t edges_group_;

  int rank_;

  std::string fname_;
  bool first_frame_;

  // Flags about datasets.
  DatasetState has_velocity_;
  DatasetState has_force_;
  DatasetState has_id_group_;
  DatasetState has_box_;

  bool file_opened_;

  // Current frame indicator.
  int idx_frame_;
  int max_idx_frame_;

  // Number of particles. This is static among time.
  int N_particles_;
  //
  int vec_components_;

  // Box matrix.
  Eigen::Matrix3d m;
};

template <class Bead_T, class Molecule_T, class Topology_T>
H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::H5MDTrajectoryReader() {
  has_velocity_ = H5MDTrajectoryReader::NONE;
  has_force_ = H5MDTrajectoryReader::NONE;
  has_id_group_ = H5MDTrajectoryReader::NONE;
  has_box_ = H5MDTrajectoryReader::NONE;
}

template <class Bead_T, class Molecule_T, class Topology_T>
H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::~H5MDTrajectoryReader() {
  if (file_opened_) {
    H5Fclose(file_id_);
    file_opened_ = false;
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
bool H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::Open(
    const std::string &file) {
  // Checks if we deal with hdf5 file.
  if (!H5Fis_hdf5(file.c_str())) {
    std::cout << file << " is not recognise as HDF5 file format" << std::endl;
    return false;
  }
  file_id_ = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  file_opened_ = true;

  // Check the version of the file.
  hid_t g_h5md = H5Gopen(file_id_, "h5md", H5P_DEFAULT);
  CheckError(g_h5md, "Unable to open /h5md group.");

  hid_t at_version = H5Aopen(g_h5md, "version", H5P_DEFAULT);
  CheckError(at_version, "Unable to read version attribute.");
  int version[2];
  H5Aread(at_version, H5Aget_type(at_version), &version);
  if (version[0] != 1 || (version[0] == 1 && version[1] > 0)) {
    std::cout << "Major version " << version[0] << std::endl;
    std::cout << "Minor version " << version[1] << std::endl;
    throw std::ios_base::failure("Wrong version of H5MD file.");
  }

  // Clean up.
  H5Aclose(at_version);
  H5Gclose(g_h5md);

  // Checks if the particles group exists and what is the number of members.
  particle_group_ = H5Gopen(file_id_, "particles", H5P_DEFAULT);
  CheckError(particle_group_, "Unable to open /particles group.");
  hsize_t num_obj = 0;
  H5Gget_num_objs(particle_group_, &num_obj);
  if (num_obj == 0) {
    throw std::ios_base::failure("The particles group is empty.");
  }

  first_frame_ = true;

  // Handle errors by internal check up.
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  return true;
}

template <class Bead_T, class Molecule_T, class Topology_T>
void H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::Close() {
  if (file_opened_) {
    H5Fclose(file_id_);
    file_opened_ = false;
  }
}

template <class Bead_T, class Molecule_T, class Topology_T>
void H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::Initialize(
    Topology_T &top) {
  std::string *particle_group_name_ = new std::string(top.getParticleGroup());
  if (*particle_group_name_ == "")
    throw std::ios_base::failure(
        "Missing particle group in topology. Please set `h5md_particle_group` "
        "tag with `name` attribute set to the particle group.");
  std::string *position_group_name =
      new std::string(*particle_group_name_ + "/position");
  atom_position_group_ =
      H5Gopen(particle_group_, position_group_name->c_str(), H5P_DEFAULT);
  CheckError(atom_position_group_,
             "Unable to open " + *position_group_name + " group");

  idx_frame_ = -1;
  ds_atom_position_ = H5Dopen(atom_position_group_, "value", H5P_DEFAULT);
  CheckError(ds_atom_position_,
             "Unable to open " + *position_group_name + "/value dataset");

  // Reads the box information.
  std::string *box_gr_name = new std::string(*particle_group_name_ + "/box");
  hid_t g_box = H5Gopen(particle_group_, box_gr_name->c_str(), H5P_DEFAULT);
  CheckError(g_box, "Unable to open " + *box_gr_name + " group");
  hid_t at_box_dimension = H5Aopen(g_box, "dimension", H5P_DEFAULT);
  CheckError(at_box_dimension, "Unable to open dimension attribute.");
  int dimension;
  H5Aread(at_box_dimension, H5Aget_type(at_box_dimension), &dimension);
  if (dimension != 3) {
    throw std::ios_base::failure("Wrong dimension " +
                                 boost::lexical_cast<std::string>(dimension));
  }
  // TODO: check if boundary is periodic.
  std::string *box_edges_name =
      new std::string(*particle_group_name_ + "/box/edges");
  if (GroupExists(particle_group_, *box_edges_name)) {
    g_box = H5Gopen(particle_group_, box_gr_name->c_str(), H5P_DEFAULT);
    edges_group_ = H5Gopen(g_box, "edges", H5P_DEFAULT);
    ds_edges_group_ = H5Dopen(edges_group_, "value", H5P_DEFAULT);
    std::cout << "H5MD: has /box/edges" << std::endl;
    has_box_ = H5MDTrajectoryReader::TIMEDEPENDENT;
  } else {
    std::cout << "H5MD: static box" << std::endl;
    hid_t ds_edges = H5Dopen(g_box, "edges", H5P_DEFAULT);
    CheckError(ds_edges, "Unable to open /box/edges");
    double *box = new double[3];
    ReadStaticData(ds_edges, H5T_NATIVE_DOUBLE, box);
    std::cout << "H5MD: Found box " << box[0] << " x " << box[1] << " x "
              << box[2] << std::endl;
    // Sets box size.
    m = Eigen::Matrix3d::Zero();
    m(0, 0) = box[0];
    m(1, 1) = box[1];
    m(2, 2) = box[2];
    top.setBox(m);
    has_box_ = H5MDTrajectoryReader::STATIC;
  }
  H5Gclose(g_box);

  // Gets the force group.
  std::string *force_group_name =
      new std::string(*particle_group_name_ + "/force");
  if (GroupExists(particle_group_, *force_group_name)) {
    atom_force_group_ =
        H5Gopen(particle_group_, force_group_name->c_str(), H5P_DEFAULT);
    ds_atom_force_ = H5Dopen(atom_force_group_, "value", H5P_DEFAULT);
    has_force_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    std::cout << "H5MD: has /force" << std::endl;
  } else {
    has_force_ = H5MDTrajectoryReader::NONE;
  }

  // Gets the velocity group.
  std::string *velocity_group_name =
      new std::string(*particle_group_name_ + "/velocity");
  if (GroupExists(particle_group_, *velocity_group_name)) {
    atom_velocity_group_ =
        H5Gopen(particle_group_, velocity_group_name->c_str(), H5P_DEFAULT);
    ds_atom_velocity_ = H5Dopen(atom_velocity_group_, "value", H5P_DEFAULT);
    has_velocity_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    std::cout << "H5MD: has /velocity" << std::endl;
  } else {
    has_velocity_ = H5MDTrajectoryReader::NONE;
  }

  // Gets the id group so that the atom id is taken from this group.
  std::string *id_group_name = new std::string(*particle_group_name_ + "/id");
  if (GroupExists(particle_group_, *id_group_name)) {
    atom_id_group_ =
        H5Gopen(particle_group_, id_group_name->c_str(), H5P_DEFAULT);
    ds_atom_id_ = H5Dopen(atom_id_group_, "value", H5P_DEFAULT);
    has_id_group_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    std::cout << "H5MD: has /id group" << std::endl;
  } else {
    has_id_group_ = H5MDTrajectoryReader::NONE;
  }

  // Gets number of particles and dimensions.
  hid_t fs_atom_position_ = H5Dget_space(ds_atom_position_);
  CheckError(fs_atom_position_, "Unable to open atom position space.");
  hsize_t dims[3];
  rank_ = H5Sget_simple_extent_dims(fs_atom_position_, dims, NULL);
  N_particles_ = dims[1];
  vec_components_ = dims[2];
  max_idx_frame_ = dims[0] - 1;

  // TODO: reads mass, charge and particle type.

  if (has_id_group_ == H5MDTrajectoryReader::NONE && top.BeadCount() > 0 &&
      static_cast<size_t>(N_particles_) != top.BeadCount()) {
    std::cout << "Warning: The number of beads (" << N_particles_ << ")";
    std::cout
        << " in the trajectory is different than defined in the topology ("
        << top.BeadCount() << ")" << std::endl;
    std::cout << "The number of beads from topology will be used!" << std::endl;
    N_particles_ = top.BeadCount();
  }

  delete id_group_name;
  delete velocity_group_name;
  delete force_group_name;
  delete box_edges_name;
  delete particle_group_name_;
}

template <class Bead_T, class Molecule_T, class Topology_T>
bool H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::FirstFrame(
    boost::any top_any) {  // NOLINT const
                           // reference

  if (first_frame_) {

    if (typeid(Topology_T *) != top_any.type()) {
      throw std::runtime_error(
          "Error Cannot read topology using h3md trajectory reader first "
          "frame, incorrect topology type provided.");
    }
    Topology_T &top = *boost::any_cast<Topology_T *>(top_any);
    first_frame_ = false;
    Initialize(top);
  }
  NextFrame(top_any);
  return true;
}

/// Reading the data.
template <class Bead_T, class Molecule_T, class Topology_T>
bool H5MDTrajectoryReader<Bead_T, Molecule_T, Topology_T>::NextFrame(
    boost::any conf_any) {  // NOLINT const
                            // reference
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using h5mdtrajectory reader next frame, "
        "incorrect topology type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  // Reads the position row.
  idx_frame_++;
  if (idx_frame_ > max_idx_frame_) return false;

  // Set volume of box because top on workers somehow does not have this
  // information.
  top.setBox(m);

  std::cout << '\r' << "Reading frame: " << idx_frame_;
  std::cout.flush();
  double *positions;
  double *forces = NULL;
  double *velocities = NULL;
  int *ids = NULL;

  try {
    positions = ReadVectorData<double>(ds_atom_position_, H5T_NATIVE_DOUBLE,
                                       idx_frame_);
  } catch (const std::runtime_error &e) {
    return false;
  }

  if (has_velocity_ != H5MDTrajectoryReader::NONE) {
    velocities = ReadVectorData<double>(ds_atom_velocity_, H5T_NATIVE_DOUBLE,
                                        idx_frame_);
  }

  if (has_force_ != H5MDTrajectoryReader::NONE) {
    forces =
        ReadVectorData<double>(ds_atom_force_, H5T_NATIVE_DOUBLE, idx_frame_);
  }

  if (has_id_group_ != H5MDTrajectoryReader::NONE) {
    ids = ReadScalarData<int>(ds_atom_id_, H5T_NATIVE_INT, idx_frame_);
  }

  // Process atoms.
  for (int at_idx = 0; at_idx < N_particles_; at_idx++) {
    double x, y, z;
    int array_index = at_idx * vec_components_;
    x = positions[array_index];
    y = positions[array_index + 1];
    z = positions[array_index + 2];
    // Set atom id, or it is an index of a row in dataset or from id dataset.
    int atom_id = at_idx;
    if (has_id_group_ != H5MDTrajectoryReader::NONE) {
      if (ids[at_idx] == -1)  // ignore values where id == -1
        continue;
      atom_id = ids[at_idx] - 1;
    }

    // Topology has to be defined in the xml file or in other
    // topology files. The h5md only stores the trajectory data.
    Bead_T *b = top.getBead(atom_id);
    if (b == NULL)
      throw std::runtime_error("Bead not found: " +
                               boost::lexical_cast<std::string>(atom_id));

    b->setPos(Eigen::Vector3d(x, y, z));
    if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double vx, vy, vz;
      vx = velocities[array_index];
      vy = velocities[array_index + 1];
      vz = velocities[array_index + 2];
      b->setVel(Eigen::Vector3d(vx, vy, vz));
    }

    if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double fx, fy, fz;
      fx = forces[array_index];
      fy = forces[array_index + 1];
      fz = forces[array_index + 2];
      b->setF(Eigen::Vector3d(fx, fy, fz));
    }
  }

  // Clean up pointers.
  delete[] positions;
  if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] forces;
  if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] velocities;
  if (has_id_group_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] ids;

  return true;
}

}  // namespace csg
}  // namespace votca
#endif  // SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
