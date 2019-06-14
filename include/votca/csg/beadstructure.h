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
#ifndef _VOTCA_CSG_BEADSTRUCTURE_H
#define _VOTCA_CSG_BEADSTRUCTURE_H

#include <cassert>
#include <iostream>
#include <unordered_map>

#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>

namespace votca {
namespace csg {

/**
 * \brief Designed to determine if the structure beads passed in
 *
 * Essentially it will have the functionality to determine if the stored beads
 * make up a single molecule. It can also break the stored beads up into
 * molecules. It can compare two bead structures and determine if they are the
 * same structure. The Ids of a bead will not determine if a structure is
 * equivalent or not. Each bead must have a unique Id.
 *
 * E.g. If a beadstructure is composed of the following:
 *
 * BeadStructure 1
 * Type:"A" Id:1 ---- Type:"B" Id:2
 *
 * BeadStucture 2
 * Type:"A" Id:3 ---- Type:"B" Id:4
 *
 * Here BeadStructure 1 and 2 will compare equal
 *
 **/
template <class T>
class BeadStructure {
 public:
  ~BeadStructure() {}

  typedef T *bead_t;
  typedef typename std::unordered_map<int, T *>::iterator iterator;
  typename std::unordered_map<int, T *>::iterator begin() {
    return beads_.begin();
  }
  typename std::unordered_map<int, T *>::iterator end() { return beads_.end(); }
  typename std::unordered_map<int, T *>::iterator begin() const {
    return beads_.begin();
  }
  typename std::unordered_map<int, T *>::iterator end() const {
    return beads_.end();
  }
  /**
   * \brief Determine if the bead structure consists of a single connected
   * structure
   *
   * This function will determine if all beads in the structure are connected
   * somehow to everyother bead. The connection does not have to be direct
   *
   * @return - returns a boolean true if it is a single Structure
   **/
  bool isSingleStructure();

  /**
   * \brief returns the number of beads in the bead structure
   **/
  size_t BeadCount() const noexcept { return beads_.size(); }

  /**
   * \brief add a bead to the bead structure
   *
   * The same bead cannot be added twice.
   **/
  void AddBead(T &bead);

  /**
   * \brief Get the bead with the specified id
   **/
  T &getBead(int id);

  std::vector<T *> getBeads() const;
  /**
   * @brief Grabs the const version of the bead pointer
   *
   * @param[in] id grabs the bead with the provided `id`
   *
   * @return const pointer to the bead
   */
  const T &getBead(int id) const;
  /**
   * \brief Create a connection between two beads in the structure
   *
   * A bead cannot be connected to itself. It also may not be connected to a
   * bead that has not yet been added to the structure.
   **/
  void ConnectBeads(const int &bead1_id, const int &bead2_id);

  /**
   * \brief Return a vector of all the beads neighboring the bead_id
   **/
  std::vector<T *> getNeighBeads(int bead_id);

  /**
   * @brief Returns a vector containing all the bead ids in the structure
   *
   * @return vector of ints
   */
  std::vector<int> getBeadIds() const;

  /**
   * @brief Returns a copy of the internal graph
   *
   * @return graph
   */
  tools::Graph getGraph();

  /**
   * \brief Compare the topology of two bead structures
   *
   * This function looks at how the beads are arranged within the bead structure
   * and determines if the topology is the same.
   *
   * @param[in] - beadstructure to compare with
   * @return - if the same returns true else false
   *
   **/
  bool isStructureEquivalent(BeadStructure<T> &beadstructure);

  /// Determine if a bead exists in the structure
  bool BeadExist(int bead_id) const { return beads_.count(bead_id); }

 protected:
  void InitializeGraph_();
  void CalculateStructure_();
  tools::GraphNode BaseBeadToGraphNode_(T *basebead);

  bool structureIdUpToDate = false;
  bool graphUpToDate = false;
  bool single_structureUpToDate_ = false;
  bool single_structure_ = false;
  std::string structure_id_ = "";
  tools::Graph graph_;
  std::set<tools::Edge> connections_;
  std::unordered_map<int, T *> beads_;
  std::unordered_map<int, tools::GraphNode> graphnodes_;
};

/**********************
 * Internal Functions *
 **********************/

template <class T>
void BeadStructure<T>::InitializeGraph_() {
  if (!graphUpToDate) {
    std::vector<tools::Edge> connections_vector;
    for (const tools::Edge &edge : connections_) {
      connections_vector.push_back(edge);
    }

    for (std::pair<const int, T *> &id_bead_ptr_pair : beads_) {
      graphnodes_[id_bead_ptr_pair.first] =
          BaseBeadToGraphNode_(id_bead_ptr_pair.second);
    }
    graph_ = tools::Graph(connections_vector, graphnodes_);
    graphUpToDate = true;
  }
}

template <class T>
void BeadStructure<T>::CalculateStructure_() {

  InitializeGraph_();
  if (!structureIdUpToDate) {
    structure_id_ = tools::findStructureId<tools::GraphDistVisitor>(graph_);
    structureIdUpToDate = true;
  }
}

template <class T>
tools::GraphNode BeadStructure<T>::BaseBeadToGraphNode_(T *basebead) {
  std::unordered_map<std::string, double> attributes1;
  std::unordered_map<std::string, std::string> attributes2;

  attributes1["Mass"] = basebead->getMass();
  attributes2["Type"] = basebead->getType();

  /// Add graphnodes
  tools::GraphNode graphnode;
  graphnode.setDouble(attributes1);
  graphnode.setStr(attributes2);

  return graphnode;
}

/***************************
 * Public Facing Functions *
 ***************************/

template <class T>
void BeadStructure<T>::AddBead(T &bead) {
  if (beads_.count(bead.getId())) {
    std::string err = "Cannot add bead with Id ";
    err += std::to_string(bead.getId());
    err += " because each bead must have a unique Id and a bead with that Id ";
    err += "already exists within the beadstructure";
    throw std::invalid_argument(err);
  }
  size_t numberOfBeads = beads_.size();
  beads_[bead.getId()] = &bead;

  if (numberOfBeads != beads_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

template <class T>
void BeadStructure<T>::ConnectBeads(const int &bead1_id, const int &bead2_id) {
  if (!(beads_.count(bead1_id)) || !(beads_.count(bead2_id))) {
    std::string err =
        "Cannot connect beads in bead structure that do not exist";
    throw std::invalid_argument(err);
  }
  if (bead1_id == bead2_id) {
    std::string err = "Beads cannot be self-connected";
    throw std::invalid_argument(err);
  }
  size_t numberOfConnections = connections_.size();
  connections_.insert(tools::Edge(bead1_id, bead2_id));
  if (numberOfConnections != connections_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

template <class T>
tools::Graph BeadStructure<T>::getGraph() {
  InitializeGraph_();
  return graph_;
}

template <class T>
bool BeadStructure<T>::isSingleStructure() {

  InitializeGraph_();
  if (single_structureUpToDate_ == false) {
    std::vector<int> vertices = graph_.getVertices();
    if (vertices.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    // Choose first vertex that is actually in the graph as the starting vertex
    tools::Graph_BF_Visitor gv_breadth_first;
    gv_breadth_first.setStartingVertex(vertices.at(0));
    if (!singleNetwork(graph_, gv_breadth_first)) {
      single_structure_ = false;
      return single_structure_;
    }
    if (beads_.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    if (vertices.size() != beads_.size()) {
      single_structure_ = false;
      return single_structure_;
    }
    single_structure_ = true;
    single_structureUpToDate_ = true;
  }
  return single_structure_;
}

template <class T>
bool BeadStructure<T>::isStructureEquivalent(BeadStructure<T> &beadstructure) {
  if (!structureIdUpToDate) {
    CalculateStructure_();
  }
  if (!beadstructure.structureIdUpToDate) {
    beadstructure.CalculateStructure_();
  }
  return structure_id_.compare(beadstructure.structure_id_) == 0;
}

template <class T>
std::vector<T *> BeadStructure<T>::getNeighBeads(int bead_id) {
  if (!graphUpToDate) {
    InitializeGraph_();
  }
  std::vector<int> neighbor_ids = graph_.getNeighVertices(bead_id);
  std::vector<T *> neighbeads;
  for (int &node_id : neighbor_ids) {
    neighbeads.push_back(beads_[node_id]);
  }
  return neighbeads;
}

template <class T>
T &BeadStructure<T>::getBead(int bead_id) {
  assert(beads_.count(bead_id));
  return *beads_[bead_id];
}

template <class T>
std::vector<T *> BeadStructure<T>::getBeads() const {
  std::vector<T *> beads;
  for (std::pair<int, T *> bead_pr : beads_) {
    beads.push_back(bead_pr.second);
  }
  return beads;
}

template <class T>
const T &BeadStructure<T>::getBead(int bead_id) const {
  assert(beads_.count(bead_id));
  return *beads_.at(bead_id);
}

template <class T>
std::vector<int> BeadStructure<T>::getBeadIds() const {
  std::vector<int> bead_ids;
  for (const std::pair<const int, T *> &id_bead_ptr_pair : beads_) {
    bead_ids.push_back(id_bead_ptr_pair.first);
  }
  return bead_ids;
}

}  // namespace csg
}  // namespace votca

#endif
