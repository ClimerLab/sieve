#ifndef EDGE_H
#define EDGE_H

#include <cstdlib>

class Edge {
  private:
  std::size_t target_id;
  double weight;

  public:
    Edge();
    Edge(const std::size_t, const double);
    ~Edge();

    void setTargetId(const std::size_t);
    void setWeight(const double);

    std::size_t getTargetId() const;
    double getWeight() const;
};

#endif
