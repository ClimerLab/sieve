#include "Edge.h"

Edge::Edge() :  target_id(0),
                weight(0.0) {}

Edge::Edge(const std::size_t _target_id,
           const double _weight) :  target_id(_target_id),
                                    weight(_weight) {}

Edge::~Edge() {}

void Edge::setTargetId(const std::size_t _target_id) {
  target_id = _target_id;
}

void Edge::setWeight(const double _weight) {
  weight = _weight;
}

std::size_t Edge::getTargetId() const {
  return target_id;
}

double Edge::getWeight() const {
  return weight;
}