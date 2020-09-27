#ifndef TRUSS_EDGE_H
#define TRUSS_EDGE_H

class Edge {
 public:
  int s, e;
  int support;
  int id;
  bool alive;
  int pos;
  bool deleted;
  int trussClass;
  int originSupport;
  Edge(int edgeId, int i, int j) {
    id = edgeId;
    s = i;
    e = j;
    alive = true;
    deleted = false;
    trussClass = -1;
  }
};

class RawEdge {
 public:
  int s, e, support;
  RawEdge(int s, int e, int support) {
    this->s = s;
    this->e = e;
    this->support = support;
  }
};

#endif  // TRUSS_EDGE_H
