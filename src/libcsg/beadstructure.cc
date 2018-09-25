/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <iostream>
#include <unordered_map>
#include <votca/csg/beadstructure.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>

using namespace votca::csg;
using namespace votca::tools;
using namespace std;

/**********************
 * Internal Functions *
 **********************/

shared_ptr<GraphNode> BaseBeadToGraphNode(shared_ptr<BaseBead> basebead) {
  unordered_map<string, double> attributes1;
  unordered_map<string, string> attributes2;

  attributes1["Mass"] = basebead->getMass();
  attributes2["Name"] = basebead->getName();

  /// Add graphnodes
  GraphNode graphnode;
  graphnode.setDouble(attributes1);
  graphnode.setStr(attributes2);

  return make_shared<GraphNode>(graphnode);
}

void BeadStructure::InitializeGraph_() {
  if (!graphUpToDate) {
    vector<Edge> connections_vector;
    for (auto edge : connections_) {
      connections_vector.push_back(edge);
    }

    unordered_map<int, GraphNode> graphnodes;
    for (auto id_bead_ptr_pair : beads_) {
      graphnodes_[id_bead_ptr_pair.first] =
          BaseBeadToGraphNode(id_bead_ptr_pair.second);
      graphnodes[id_bead_ptr_pair.first] =
          *(graphnodes_[id_bead_ptr_pair.first]);
    }
    graph_ = make_shared<Graph>(Graph(connections_vector, graphnodes));
    graphUpToDate = true;
  }
}

void BeadStructure::CalculateStructure_() {

  InitializeGraph_();
  if (!structureIdUpToDate) {
    findStructureId<GraphDistVisitor>(*graph_);
    structureIdUpToDate = true;
  }
}

/***************************
 * Public Facing Functions *
 ***************************/

void BeadStructure::AddBead(shared_ptr<BaseBead> bead) {
  if (beads_.count(bead->getId())) {
    string err = "Cannot add bead with Id ";
    err += to_string(bead->getId());
    err += " because each bead must have a unique Id and a bead with that Id ";
    err += "already exists within the beadstructure";
    throw invalid_argument(err);
  }
  auto numberOfBeads = beads_.size();
  beads_[bead->getId()] = bead;
  if (numberOfBeads != beads_.size()) {
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

void BeadStructure::ConnectBeads(int bead1_id, int bead2_id) {
  if (!(beads_.count(bead1_id)) || !(beads_.count(bead2_id))) {
    string err = "Cannot connect beads in bead structure that do not exist";
    throw invalid_argument(err);
  }
  if (bead1_id == bead2_id) {
    string err = "Beads cannot be self-connected";
    throw invalid_argument(err);
  }
  auto numberOfConnections = connections_.size();
  connections_.insert(Edge(bead1_id, bead2_id));
  if (numberOfConnections != connections_.size()) {
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

bool BeadStructure::isSingleMolecule() {

  InitializeGraph_();
  auto vertices = graph_->getVertices();
  if (vertices.size() == 0)
    return false;

  // Choose first vertex that is actually in the graph as the starting vertex
  Graph_BF_Visitor gv_breadth_first;
  gv_breadth_first.setStartingVertex(vertices.at(0));

  if (!singleNetwork(*graph_, gv_breadth_first))
    return false;
  if (beads_.size() == 0)
    return false;
  if (vertices.size() != beads_.size())
    return false;
  return true;
}

bool BeadStructure::isStructureEquivalent(BeadStructure &beadstructure) {
  if (!structureIdUpToDate)
    CalculateStructure_();
  if (!beadstructure.structureIdUpToDate)
    beadstructure.CalculateStructure_();
  return *graph_ == *beadstructure.graph_;
}

vector<shared_ptr<BaseBead>> BeadStructure::getNeighBeads(int index) {
  if (!graphUpToDate)
    InitializeGraph_();
  auto neighbor_ids = graph_->getNeighVertices(index);
  vector<shared_ptr<BaseBead>> neighbeads;
  for (auto node_id : neighbor_ids) {
    neighbeads.push_back(beads_[node_id]);
  }
  return neighbeads;
}

shared_ptr<BaseBead> BeadStructure::getBead(int id) {
  if(id<0){
    string err = "bead with negative id " + to_string(id);
    throw invalid_argument(err);
  }
  if(!beads_.count(id)){
    string err = "bead with id: " + to_string(id) + " is not found.";
    throw invalid_argument(err);
  }
  return beads_[id];
}

vector<int> BeadStructure::getIdsOfBeadsWithName(const string &name){ 
  vector<int> ids;
  auto iterator = beads_.begin();
  while(iterator!=beads_.end()){
    if(iterator->second->getName() == name){
      ids.push_back(iterator->second->getId());
    }
    ++iterator;
  }
  return ids;
}

vector<int> BeadStructure::getBeadIds(){
  vector<int> ids;
  auto iterator = beads_.begin();
  while(iterator!=beads_.end()){
    ids.push_back(iterator->second->getId());
    ++iterator;
  }
  return ids;
}

string BeadStructure::getBeadName(int id){
  if(id<0){
    string err = "bead with negative id " + to_string(id);
    throw invalid_argument(err);
  }
  if(!beads_.count(id)){
    string err = "bead with id: " + to_string(id) + " is not found.";
    throw invalid_argument(err);
  }
  return beads_[id]->getName();
}

vector<shared_ptr<BeadStructure>> BeadStructure::breakIntoMolecules() {
  vector<shared_ptr<BeadStructure>> structures;
  if (!graphUpToDate)
    InitializeGraph_();
  auto sub_graphs = decoupleIsolatedSubGraphs(*(this->graph_));
  for (auto sub_graph_it = sub_graphs.begin(); sub_graph_it != sub_graphs.end();
       ++sub_graph_it) {
    auto sub_graph_edges = (*sub_graph_it)->getEdges();
    auto sub_graph_vertices = (*sub_graph_it)->getVertices();
    BeadStructure beadstructure;
    for (auto vertex : sub_graph_vertices)
      beadstructure.AddBead(beads_[vertex]);
    for (auto edge : sub_graph_edges)
      beadstructure.ConnectBeads(edge.getV1(), edge.getV2());
    structures.push_back(make_shared<BeadStructure>(beadstructure));
  }
  return structures;
}
