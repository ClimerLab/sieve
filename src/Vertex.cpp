#include "Vertex.h"
#include "assert.h"
#include <stdio.h>

Vertex::Vertex() :  vertex_weight(0.0) {}

Vertex::~Vertex() {}

std::size_t Vertex::getDegree() const {
  return edges.size();
}

double Vertex::getVertexWeight() const {
  return vertex_weight;
}

double Vertex::getEdgeWeight(const std::size_t index) const {
  assert(index < edges.size());

  return edges[index].getWeight();
}

std::size_t Vertex::getTargetId(const std::size_t index) const {
  assert(index < edges.size());

  return edges[index].getTargetId();
}

void Vertex::setEdgeWeight(const std::size_t index, const double weight) {
  assert(index < edges.size());
  edges[index].setWeight(weight);
}

bool Vertex::addEdge(const std::size_t target_id, double weight) {
  bool edge_added = false;
  Edge new_edge(target_id, weight);

  if(!isAdjacent(target_id)) {
    edges.push_back(new_edge);
    vertex_weight += weight;
    edge_added = true;
  }        

  return edge_added;
}

bool Vertex::removeEdgeById(const std::size_t target_id) {
  bool edge_removed = false;

  for(std::size_t i = 0; i < getDegree(); ++i) {
    if(edges[i].getTargetId() == target_id) {
      edges.erase(edges.begin() + i);
      edge_removed = true;
      break;
    }
  }

  return edge_removed;
}

bool Vertex::isAdjacent(const std::size_t target_id) const {
  bool id_found = false;
  for(std::size_t i = 0; i < edges.size(); ++i) {
    if(edges[i].getTargetId() == target_id) {
      id_found = true;
      break;            
    }
  }

  return id_found;
}

void Vertex::swapEdges(const std::size_t index1, const std::size_t index2) {
  assert(index1 < edges.size());
  assert(index2 < edges.size());

  std::size_t tmp_id = edges[index1].getTargetId();
  double tmp_weight = edges[index1].getWeight();

  edges[index1].setTargetId(edges[index2].getTargetId());
  edges[index1].setWeight(edges[index2].getWeight());

  edges[index2].setTargetId(tmp_id);
  edges[index2].setWeight(tmp_weight);
}

std::vector<std::size_t> Vertex::getIds() const {
  std::vector<std::size_t> ids(edges.size());

  for(std::size_t i = 0; i < edges.size(); ++i)
    ids[i] = edges[i].getTargetId();

  return ids;
}