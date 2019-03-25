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

#ifndef VOTCA_CSG_FMATCH_H
#define VOTCA_CSG_FMATCH_H

#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/tools/cubicspline.h>
#include <votca/tools/property.h>

namespace votca {
namespace csg {
namespace TOOLS = votca::tools;
/**
    \brief Implements force matching algorithm using cubic spline basis set
 *
 *  Force matching method to obtain a coarse-grained force field is implemented
 *  using cubic spline basis set. Block averaging over trajectory blocks
 *  is used for calculating CG forces and their errors.
 *
 * \todo force matching needs a big cleanup!
 **/

class CGForceMatching : public CsgApplication {
 public:
  std::string ProgramName() { return "csg_fmatch"; }
  void HelpText(std::ostream &out) {
    out << "Perform force matching (also called multiscale coarse-graining)";
  }

  bool DoTrajectory() { return true; }
  bool DoMapping() { return true; }

  void Initialize(void);
  bool EvaluateOptions();

  /// \brief called before the first frame
  void BeginEvaluate(CSG_Topology *top, CSG_Topology *top_atom);
  /// \brief called after the last frame
  void EndEvaluate();
  /// \brief called for each frame which is mapped
  void EvalConfiguration(CSG_Topology *conf, CSG_Topology *conf_atom = 0);
  /// \brief load options from the input file
  void LoadOptions(const std::string &file);

 protected:
  /// \brief structure, which contains tools::CubicSpline object with related
  /// parameters
  struct SplineInfo {
    /// \brief constructor
    SplineInfo(int index, bool bonded_, int matr_pos_,
               tools::Property *options);
    /// \brief number of spline functions
    int num_splinefun;
    /// \brief number of spline grid points
    int num_gridpoints;
    /// \brief number of grid points for output
    int num_outgrid;
    /// \brief interaction index
    int splineIndex;
    /// \brief true for bonded interactions, false for non-bonded
    bool bonded;
    /// \brief true if tabulated forces are periodic (e.g. for dihedral
    /// interactions)
    bool periodic;
    /// \brief true if non-bonded interaction is threebody interaction
    bool threebody;
    /// \brief additional variables for treating cutoff of
    double a;
    double sigma;
    double gamma;
    /// \brief tools::CubicSpline object
    tools::CubicSpline Spline;
    /// \brief position in the _A matrix (first coloumn which is occupied with
    /// this particular spline)
    int matr_pos;
    /// \brief dx for output. Calculated in the code
    double dx_out;
    /// \brief only for non-bonded interactions (seems like it is not used?)
    std::pair<int, int> beadTypes;

    /// \brief Result of 1 block calculation for f
    Eigen::VectorXd block_res_f;
    /// \brief Result of 1 block calculation for f''
    Eigen::VectorXd block_res_f2;
    /// \brief Final result: average over all blocks
    Eigen::VectorXd result;
    /// \brief accuracy of the final result
    Eigen::VectorXd error;
    /// \brief sum of all block_res (used to calculate error)
    Eigen::VectorXd resSum;
    /// \brief sum of all squares of block_res (used to calculate error)
    Eigen::VectorXd resSum2;

    // only needed for 3body nonbonded interactions as here force and potential
    // are calculated simultaneously
    /// \brief Final result of derivatives: average over all blocks
    Eigen::VectorXd resultDer;
    /// \brief accuracy of the final result
    Eigen::VectorXd errorDer;
    /// \brief sum of all block_res
    Eigen::VectorXd resSumDer;
    /// \brief sum of all squares of block_res (used to calculate error)
    Eigen::VectorXd resSumDer2;

    /// \brief Spline Name
    std::string splineName;
    /// \brief for non-bonded interactions: types of beads involved (type3 only
    /// used if threebody interaction)
    std::string type1, type2, type3;  //

    /// \brief pointer to tools::Property object to hande input options
    tools::Property *_options;
  };
  /// \brief tools::Property object to hande input options
  tools::Property _options;
  /// \brief std::list of bonded interactions
  std::list<tools::Property *> _bonded;
  /// \brief std::list of non-bonded interactions
  std::list<tools::Property *> _nonbonded;

  typedef std::vector<SplineInfo> SplineContainer;
  /// \brief std::vector of SplineInfo * for all interactions
  SplineContainer _splines;

  /// \brief matrix used to store force matching equations
  Eigen::MatrixXd _A;
  /// \brief std::vector used to store reference forces on CG beads (from
  /// atomistic simulations)
  Eigen::VectorXd _b;
  /// \brief Solution of matrix equation _A * _x = _b : CG force-field
  /// parameters
  Eigen::VectorXd _x;  //
  /// \brief Additional matrix to handle constrained least squares fit
  /// contains constraints, which allow to get a real (smooth) spline (see VOTCA
  /// paper)
  Eigen::MatrixXd _B_constr;

  /// \brief Counter for trajectory frames
  int _frame_counter;
  /// \brief Number of CG beads
  int _nbeads;

  /// \brief accuracy for evaluating the difference in bead positions
  double _dist;

  /// \brief Flag: true for constrained least squares, false for simple least
  /// squares
  bool _constr_least_sq;
  /// \brief used in EvalConf to distinguish constrained and simple least
  /// squares
  int _least_sq_offset;
  /// \brief Number of frames used in one block for block averaging
  int _nframes;
  /// \brief Current number of blocks
  int _nblocks;

  /// \brief Counters for lines and coloumns in _B_constr
  int _line_cntr, _col_cntr;

  bool _has_existing_forces;

  /// \brief Solves FM equations for one block and stores the results for
  /// further processing
  void FmatchAccumulateData();
  /// \brief Assigns smoothing conditions to matrices _A and _B_constr
  void FmatchAssignSmoothCondsToMatrix(Eigen::MatrixXd &Matrix);
  /// \brief For each trajectory frame writes equations for bonded interactions
  /// to matrix _A
  void EvalBonded(CSG_Topology *conf, SplineInfo *sinfo);
  /// \brief For each trajectory frame writes equations for non-bonded
  /// interactions to matrix _A
  void EvalNonbonded(CSG_Topology *conf, SplineInfo *sinfo);
  /// \brief For each trajectory frame writes equations for non-bonded threebody
  /// interactions to matrix _A
  void EvalNonbonded_Threebody(CSG_Topology *conf, SplineInfo *sinfo);
  /// \brief Write results to output files
  void WriteOutFiles();

  void OpenForcesTrajectory();

  CSG_Topology _top_force;
  TrajectoryReader *_trjreader_force;
};
}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_FMATCH_H
