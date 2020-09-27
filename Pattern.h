#ifndef DAVINCIPATTERNGEN_PATTERN_H
#define DAVINCIPATTERNGEN_PATTERN_H

#include <unordered_set>
#include <vector>

#include "Edge.h"

using namespace std;
enum PatType {
  OneEdge = 0,
  TwoEdge = 1,
  Triangle = 2,
  Truss = 3,           // truss-generated patterns
  TrussExtend = 4,     // now unused
  TrussCombined1 = 5,  // k1-side-k2-main
  TrussCombined2 = 6,  // k1-side-k2-side
  InRemain = 7,        // outside truss patterns
  OutsideStar = 8,
  OutsideCombinedStar = 9,
  OutsideLine = 10,
  OutsideCircle = 11,
  OutsideIndividual = 12
};

class PatEdge {
 public:
  int s, e;
  PatEdge(int s, int e) {
    this->s = s;
    this->e = e;
  }
};

class StarState {
 public:
  vector<int> nodeId;
  vector<short> degree;
  int cursize;
  unordered_set<int> myv;
  StarState() {
    nodeId.clear();
    degree.clear();
    myv.clear();
    cursize = 0;
  }
};

class Pattern {
 public:
  static int curId;
  int score = 0;
  int scorePart[10];
  int type;
  int id;
  vector<PatEdge> edges;
  int nodesCnt;
  int trussClass;
  double scoreCoverage;
  double scoreCognition;
  double finalScore;
  pair<int, int> detilType;
  double dynamicScore;
  int historyMaxComE;
  int finalRank;
  vector<vector<pair<int, int>>> edgesList;  // s --> [(e, edge_id)]
  vector<short> degree;  // now only use for linear star search
  int outsideCoverage;

  Pattern(int n) {
    this->nodesCnt = n;
    edgesList.resize(n);
    for (int i = 0; i < 10; ++i) scorePart[i] = 0;
    this->id = curId++;
    scoreCoverage = 0;
    scoreCognition = 0;
    finalScore = 0;
    dynamicScore = 0x3f3f3f3f;
    historyMaxComE = 0;
    // degree.resize(n, 0);
    degree.clear();
  }

  void setType(int type) { this->type = type; }
  void addEdge(int s, int e) {
    edgesList[s].emplace_back(make_pair(e, edges.size()));
    edgesList[e].emplace_back(make_pair(s, edges.size()));
    edges.emplace_back(s, e);
  }
  int getEdgeId(int s, int e) const {
    for (auto p : edgesList[s])
      if (p.first == e) return p.second;
    return -1;
  }
  void setScore(int score) { this->score = score; }
  void setTrussClass(int trussClass) { this->trussClass = trussClass; }
  int getNodesCnt() const { return nodesCnt; }
  int getEdgeCnt() const { return edges.size(); }
  int getScore() { return score; }
  int getType() { return type; }
  int getTrussClass() { return trussClass; }
};

#endif  // DAVINCIPATTERNGEN_PATTERN_H
