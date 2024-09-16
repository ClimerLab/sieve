#ifndef VERTEX_H
#define VERTEX_H

#include "Edge.h"
#include <vector>

class Vertex {
  private:
    std::vector<Edge> edges;
    double vertex_weight;

  public:
    Vertex();
    ~Vertex();

    std::size_t getDegree() const;
    double getVertexWeight() const;
    double getEdgeWeight(const std::size_t) const;
    std::size_t getTargetId(const std::size_t) const;
    void setEdgeWeight(const std::size_t, const double);

    bool addEdge(const std::size_t, double);
    bool removeEdgeById(const std::size_t);
    bool isAdjacent(const std::size_t) const;

    void swapEdges(const std::size_t, const std::size_t);
    
    std::vector<std::size_t> getIds() const;
};

#endif
