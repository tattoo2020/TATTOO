#ifndef DAVINCIPATTERNGEN_GRAPH_H
#define DAVINCIPATTERNGEN_GRAPH_H

#include <bits/stdc++.h>

#include "Edge.h"
#include "Node.h"
#include "Truss.h"

#define MAX_SUP 1000000

#define FINAL_TRUSS 1
#define FINAL_REMAIN 2
#define PROCESS_FURTHER 3

using namespace std;

struct myHash {
  size_t operator()(pair<int, int> __val) const {
    return static_cast<size_t>(__val.first ^ __val.second);
  }
};

/**
 * Now unused, only for deep truss decomposition
 */
class GraphTask {
 public:
  GraphTask() {}
  GraphTask(int type, int k, bool withSupport, bool isDirect) {
    this->k = k;
    this->withSupport = withSupport;
    this->type = type;
    this->isDirect = isDirect;
  }
  GraphTask(int type, int k, bool withSupport, vector<RawEdge>& edges,
            bool isDirect) {
    this->k = k;
    this->withSupport = withSupport;
    this->type = type;
    for (int i = 0; i < edges.size(); i++) {
      this->edges.push_back(RawEdge(edges[i].s, edges[i].e, edges[i].support));
    }
    this->isDirect = isDirect;
  }
  void outputToFile(int id, string filename) {
    ofstream fout(filename, ios::app);
    for (int i = 0; i < edges.size(); i++) {
      fout << edges[i].s << "_" << id << " " << edges[i].e << "_" << id << endl;
    }
    fout.close();
  }
  string filename;
  int k;
  bool withSupport;
  /**
   * type == 1 final truss
   * type == 2 final remain
   * type == 3 processor further
   */
  int type;
  vector<RawEdge> edges;
  bool isDirect;
};

class RemainGraphInfo {
 public:
  int nodesCnt;
  int edgesCnt;
  RemainGraphInfo(int nodesCnt, int edgesCnt) {
    this->nodesCnt = nodesCnt;
    this->edgesCnt = edgesCnt;
  }
};

class FinalTrussInfo {
 public:
  int nodesCnt;
  int edgesCnt;
  int k;
  FinalTrussInfo(int k, int nodesCnt, int edgesCnt) {
    this->nodesCnt = nodesCnt;
    this->edgesCnt = edgesCnt;
    this->k = k;
  }
};

class Graph {
 public:
  int maxSup = 0;
  static int remainGraphId;
  static int trussGraphId;
  static int remainSubgraphId;
  static int finalTrussId;
  static int finalRemainId;
  static int inTrussEdges;
  static int remainEdges;
  unordered_map<int, Node> nodesMap;  // nodeId, NodeInfo, the nodeId may not be
                                      // continuous, so use map to store them
  vector<Edge> edges;
  vector<int> binStart;
  vector<int> binSize;
  vector<int> curLength;
  vector<int> sortedEdges;
  string outputPrefix;
  static vector<int> oriDelCnt;
  static vector<int> oriTrussEdgeCnt;
  static vector<int> finalTrussEdgeCnt;
  bool isDirect;
  int edgeCnt;
  Graph() {}
  Graph(string inputFileName, bool withSupport = false, bool isDirect = true) {
    if (withSupport == false) {
      loadGraphFromEdgesFile(inputFileName);
      calcSupport();
    } else {
      loadGraphFromEdgesFileWithSupport(inputFileName);
    }
    arrangeSortedEdges();
    this->isDirect = isDirect;
  }
  Graph(vector<RawEdge>& rawEdges, bool withSupport, bool isDirect) {
    lessThanTwo = 0;
    if (withSupport) {
      loadGraphFromEdgesWithSupport(rawEdges);
      arrangeSortedEdges();
      this->isDirect = isDirect;
    } else {
      loadGraphFromEdges(rawEdges);
      calcSupport();
      arrangeSortedEdges();
      this->isDirect = isDirect;
    }
  }

  Graph(vector<Edge>& graphEdges, bool isDirect) {
    loadGraphFromVectorWithSupport(graphEdges);
    this->isDirect = isDirect;
  }
  void loadGraphFromFile(string filename);
  void loadGraphFromEdgesFile(string filename);
  void loadGraphFromEdges(vector<RawEdge>& rawEdges);
  void loadGraphFromEdgesWithSupport(vector<RawEdge>& rawEdges);
  void loadGraphFromVectorWithSupport(vector<Edge>& edges);
  void loadGraphFromEdgesFileWithSupport(string filename);
  void calcSupport();
  void arrangeSortedEdges();
  void end();
  void readSupport(string outputPrefix);
  int getKTruss(int k);
  int lessThanTwo;
  int maxSize;
  void outputEdgeInfo(string outputPrefix);
  void outputSortedEdges();
  void test();
  void convertGraph(string inputFile, string outputFile);
  void trussDecomposition(queue<GraphTask>& graphTasks, int k, int maxK);
  void outputGraphToFile(string filename);
  vector<RemainGraphInfo> generateConnectedSubgraphs(string outputFile,
                                                     int& maxRemainGraphSize);
  void outputDenseRemain();
  void annoEdgeWithTrussClass(int maxK);
  vector<int> getSurroudings(int edgeId, int maxK);
  void outputEdgesWithInfo(string outputPrefix);

 private:
  int getInter(int s, int e);
  void decSupport(int s, int e, int k);
  int getEdgeId(int s, int e);
  void decSupportForEdge(int edgeId, int k);
  bool checkDec(int s, int e);
  void floodFill(queue<GraphTask>& graphTasks, int k, int maxK);
  void outputEdgesToFileWithSupport(string filename,
                                    unordered_set<int>& edgeIds);
  void outputEdgesToFile(string filename);
  void removeStarEdges();
  bool allEdgesRemoved();
};

#endif  // TRUSS_GRAPH_H