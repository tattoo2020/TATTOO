#ifndef DAVINCIPATTERNGEN_PGGRAPH_H
#define DAVINCIPATTERNGEN_PGGRAPH_H

#include <tuple>
#include <vector>

enum EdgeType {
  RealTruss = 0,  // in Gin
  Line = 1,
  Circle = 2,
  SimpleStar = 3,
  CombinedStar = 4,
  Individual = 5,
  Default = -1
};

struct NeighbVE {
  int neighbId;
  int edgeId;
  NeighbVE(int v, int e) {
    neighbId = v;
    edgeId = e;
  }
};

class PGEdge {
 public:
  int e, s;
  PGEdge(int i, int j) {
    s = i;
    e = j;
  }
};

#endif  // DAVINCIPATTERNGEN_PGGRAPH_H
