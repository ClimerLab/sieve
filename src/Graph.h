#ifndef GRAPH_H
#define GRAPH_H

#include <cstdlib>
#include <vector>
#include <string>
#include <memory>

const bool graph_directed = false;

class Vertex;

class Graph {
  private:
    std::vector<Vertex> vertices;
    std::size_t num_edges;

    std::size_t total_num_nodes;
    std::size_t total_num_edges;
    bool weighted;

    void readGML(const std::string &);
	  void readGraph(const std::string &);

  public:
    Graph();    
    ~Graph();

    std::size_t getNumNodes() const;
    std::size_t getNumEdges() const;
    std::size_t getTotalNumNodes() const;
    std::size_t getTotalNumEdges() const;
    Vertex* getPtrToVertex(const std::size_t);

    void setTotalNumNodes(const std::size_t);
    void setTotalNumEdges(const std::size_t);

    void readFile(const std::string &);

    bool addEdge(const std::size_t, const std::size_t, const double);
    bool removeEdge(const std::size_t, const std::size_t);
    std::size_t getDegree(const std::size_t) const;

    std::vector<std::size_t> getMinNeighbors(const std::size_t, const std::size_t);
    std::vector<std::size_t> getMinSepFromNeighbors(const std::size_t, const std::size_t);
    bool isAdjacent(const std::size_t, const std::size_t) const;

    std::size_t seperateGraphIntoComponents(const std::string &, const std::string &);
    void saveAsClusterAssignment(const std::string &);
};

#endif
