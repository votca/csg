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

#include <votca/csg/beadmotif.h>

#include <votca/tools/graphalgorithm.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/**********************
 * Internal Functions *
 **********************/
bool BeadMotif::junctionExist_() {
  if (!junctionsUpToDate_) {
    junctions_ = graph_->getJunctions();
    junctionsUpToDate_ = true;
  }
  return junctions_.size() != 0;
}

bool BeadMotif::isSingle_() {
  if (BeadCount() == 1) return true;
  return false;
}

bool BeadMotif::isLine_() {
  if (junctionExist_()) return false;

  // Ensure that the degree of two of the vertices is 1
  // all other vertices must be 2
  int num_vertices_degree_1 = 0;

  auto vertices = reduced_graph_->getVertices();
  for (auto vertex : vertices) {
    int degree = reduced_graph_->getDegree(vertex);
    if (degree == 1) {
      ++num_vertices_degree_1;
    } else if (degree == 0) {
      return false;
    } else if (degree > 2) {
      return false;
    }
  }
  if (num_vertices_degree_1 != 2) return false;
  return true;
}

bool BeadMotif::isLoop_() {
  if (junctionExist_()) return false;

  // Ensure that the degree of every vertex is 2
  auto vertices = graph_->getVertices();
  for (auto vertex : vertices) {
    cout << "Checking vertex " << vertex << " has degree " << graph_->getDegree(vertex) << endl;
    if (graph_->getDegree(vertex) != 2) return false;
  }
  return true;
}

bool BeadMotif::isFusedRing_() {
  if (!junctionExist_()) return false;
  cout << "Junction exists" << endl;
  // Ensure that the degree of every vertex is 2 or greater
  auto vertices = graph_->getVertices();
  for (auto vertex : vertices) {
    if (graph_->getDegree(vertex) < 2) return false;
  }
  // If exploring from one branch of a junction lets you explore every
  // edge than it is a fused ring.
  junctionExist_();

  for( int junction : junctions_ ){
    vector<Edge> edges = reduced_graph_->getNeighEdges(junction);
    set<Edge>    all_edges_explored =
      exploreBranch(*reduced_graph_, junction, edges.at(0));
    cout << "Exploring junction " << junction << endl;
    for (const Edge& edge_next_to_junction : edges) {
      cout << "Edges next to junciton" << endl;
      cout << edge_next_to_junction << endl;
      if (!all_edges_explored.count(edge_next_to_junction)) {
        return false;
      }
    }
  }
  return true;
}
/***************************
 * Public Facing Functions *
 ***************************/

BeadMotif::motif_type BeadMotif::getType() { 
  if(!type_up_to_date_){
    CalculateType_();
  }
  return type_; 
}

void BeadMotif::CalculateType_() {
  if (BeadCount() == 0) {
    cout << "empty" << endl;
    type_ = motif_type::empty;
  } else if (isSingle_()) {
    cout << "Is single bead" << endl;
    type_ = motif_type::single_bead;
  } else if (!BeadStructure::isSingleStructure()) {
    cout << "isSingleStructure " << BeadStructure::isSingleStructure() << endl;
    cout << "Multiple structures" << endl;
    type_ = motif_type::multiple_structures;
  } else if (isLine_()) {
    cout << "Is line " << endl;
    type_ = motif_type::line;
  } else if (isLoop_()) {
    cout << "Is loop " << endl;
    type_ = motif_type::loop;
  } else if (isFusedRing_()) {
    cout << "Is fused ring " << endl;
    type_ = motif_type::fused_ring;
  } else {
    cout << "Is single structure " << endl;
    type_ = motif_type::single_structure;
  }
  type_up_to_date_ = true;
}

BaseBead* BeadMotif::getBead(int id) { return BeadStructure::getBead(id); }

void BeadMotif::AddBead(BaseBead* bead) {
  type_ = motif_type::undefined;
  BeadStructure::AddBead(bead);
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}

bool BeadMotif::isStructureEquivalent(BeadMotif& beadmotif) {
  return BeadStructure::isStructureEquivalent(beadmotif);
}

std::vector<BaseBead*> BeadMotif::getNeighBeads(int index) {
  return getNeighBeads(index);
}

void BeadMotif::InitializeGraph_(){
  BeadStructure::InitializeGraph_();
  cout << "BeadMotif calling initialize graph" << endl;
  reduced_graph_ = unique_ptr<ReducedGraph>( new ReducedGraph(move(reduceGraph(*graph_))));
  cout << *reduced_graph_ << endl;
}

int BeadMotif::BeadCount() { return BeadStructure::BeadCount(); }

void BeadMotif::ConnectBeads(int bead1_id, int bead2_id) {
  cout << "Calling Connect Beads " << endl;
  BeadStructure::ConnectBeads(bead1_id, bead2_id);
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}
}  // namespace csg
}  // namespace votca