/**
 * the truss decomposition algorithm, may contain some functions for deep truss
 * decomposition(which is useless currently)
 */
#include "Graph.h"

using namespace std;

int Graph::remainGraphId = 0;
int Graph::trussGraphId = 0;
int Graph::remainSubgraphId = 0;
int Graph::finalTrussId = 0;
int Graph::finalRemainId = 0;
int Graph::inTrussEdges = 0;
int Graph::remainEdges = 0;
vector<int> Graph::oriDelCnt(16, 0);
vector<int> Graph::oriTrussEdgeCnt(16, 0);
vector<int> Graph::finalTrussEdgeCnt(16, 0);

/**
 * the input file only contain the edges
 * @param filename
 */
void Graph::loadGraphFromEdgesFile(string filename) {
  ifstream fin(filename.c_str());
  cout << "Input the graph..." << endl;
  int s, e;
  int edgeId = 0;
  while (fin >> s >> e) {
    if (s == e) continue;
    if (nodesMap.find(s) == nodesMap.end()) {
      Node node = Node(s);
      nodesMap.insert(make_pair(s, node));
    }
    if (nodesMap.find(e) == nodesMap.end()) {
      Node node = Node(e);
      nodesMap.insert(make_pair(e, node));
    }
    if (nodesMap[s].neighbours.find(e) != nodesMap[s].neighbours.end())
      continue;
    Edge edge = Edge(edgeId, s, e);
    edges.push_back(edge);
    nodesMap[s].neighbours.insert(make_pair(e, edgeId));
    nodesMap[e].neighbours.insert(make_pair(s, edgeId));
    edgeId++;
  }
  fin.close();
  sortedEdges.resize(edges.size(), 0);
  binStart.resize(MAX_SUP, 0);
  binSize.resize(MAX_SUP, 0);
  curLength.resize(MAX_SUP, 0);
  cout << "Nodes Cnt:" << nodesMap.size() << endl;
}

/**
 * from an array of edges
 * @param rawEdges
 */
void Graph::loadGraphFromEdges(vector<RawEdge> &rawEdges) {
  cout << "Input the graph..." << endl;
  int edgeId = 0;
  for (int i = 0; i < rawEdges.size(); i++) {
    int s = rawEdges[i].s;
    int e = rawEdges[i].e;
    if (s == e) continue;
    if (nodesMap.find(s) == nodesMap.end()) {
      Node node = Node(s);
      nodesMap.insert(make_pair(s, node));
    }
    if (nodesMap.find(e) == nodesMap.end()) {
      Node node = Node(e);
      nodesMap.insert(make_pair(e, node));
    }
    if (nodesMap[s].neighbours.find(e) != nodesMap[s].neighbours.end())
      continue;
    Edge edge = Edge(edgeId, s, e);
    edges.push_back(edge);
    nodesMap[s].neighbours.insert(make_pair(e, edgeId));
    nodesMap[e].neighbours.insert(make_pair(s, edgeId));
    edgeId++;
  }
  sortedEdges.resize(edges.size(), 0);
  binStart.resize(MAX_SUP, 0);
  binSize.resize(MAX_SUP, 0);
  curLength.resize(MAX_SUP, 0);
  cout << "Nodes Cnt:" << nodesMap.size() << endl;
  cout << "Edge Cnt:" << edges.size() << endl;
}

/**
 * the input file already containes the support of each edge which have been
 * processed
 * @param filename
 */
void Graph::loadGraphFromEdgesFileWithSupport(string filename) {
  ifstream fin(filename.c_str());
  cout << "Input the graph with support..." << endl;
  int s, e, support;
  int edgeId = 0;
  int maxSupport = 0;
  while (fin >> s >> e >> support) {
    if (s == e) continue;
    if (nodesMap.find(s) == nodesMap.end()) {
      Node node = Node(s);
      nodesMap.insert(make_pair(s, node));
    }
    if (nodesMap.find(e) == nodesMap.end()) {
      Node node = Node(e);
      nodesMap.insert(make_pair(e, node));
    }
    if (nodesMap[s].neighbours.find(e) != nodesMap[s].neighbours.end())
      continue;
    Edge edge = Edge(edgeId, s, e);
    edge.support = support;
    if (support > maxSupport) maxSupport = support;
    edges.push_back(edge);
    nodesMap[s].neighbours.insert(make_pair(e, edgeId));
    nodesMap[e].neighbours.insert(make_pair(s, edgeId));
    edgeId++;
  }
  fin.close();
  sortedEdges.resize(edges.size(), 0);
  binStart.resize(MAX_SUP, 0);
  binSize.resize(MAX_SUP, 0);
  curLength.resize(MAX_SUP, 0);
  cout << "Nodes Cnt:" << nodesMap.size() << endl;
}

/**
 * the array of edges have already contains the support
 * @param rawEdges
 */
void Graph::loadGraphFromEdgesWithSupport(vector<RawEdge> &rawEdges) {
  cout << "Input the graph with support..." << endl;
  int s, e, support;
  int edgeId = 0;
  int maxSupport = 0;
  for (int i = 0; i < rawEdges.size(); i++) {
    int s = rawEdges[i].s;
    int e = rawEdges[i].e;
    if (s == e) continue;
    int support = rawEdges[i].support;
    if (nodesMap.find(s) == nodesMap.end()) {
      Node node = Node(s);
      nodesMap.insert(make_pair(s, node));
    }
    if (nodesMap.find(e) == nodesMap.end()) {
      Node node = Node(e);
      nodesMap.insert(make_pair(e, node));
    }
    if (nodesMap[s].neighbours.find(e) != nodesMap[s].neighbours.end())
      continue;
    Edge edge = Edge(edgeId, s, e);
    edge.support = support;
    if (support > maxSupport) maxSupport = support;
    edges.push_back(edge);
    nodesMap[s].neighbours.insert(make_pair(e, edgeId));
    nodesMap[e].neighbours.insert(make_pair(s, edgeId));
    edgeId++;
  }
  sortedEdges.resize(edges.size(), 0);
  binStart.resize(MAX_SUP, 0);
  binSize.resize(MAX_SUP, 0);
  curLength.resize(MAX_SUP, 0);
  cout << "Nodes Cnt:" << nodesMap.size() << endl;
  cout << "Edges Cnt:" << edges.size() << endl;
}

void Graph::loadGraphFromVectorWithSupport(vector<Edge> &graphEdges) {
  cout << "Generate the graph with support..." << endl;
  int s, e, support;
  int edgeId = 0;
  int maxSupport = 0;
  for (int i = 0; i < graphEdges.size(); i++) {
    s = graphEdges[i].s;
    e = graphEdges[i].e;
    support = graphEdges[i].support;
    if (nodesMap.find(s) == nodesMap.end()) {
      Node node = Node(s);
      nodesMap.insert(make_pair(s, node));
    }
    if (nodesMap.find(e) == nodesMap.end()) {
      Node node = Node(e);
      nodesMap.insert(make_pair(e, node));
    }
    if (nodesMap[s].neighbours.find(e) != nodesMap[s].neighbours.end())
      continue;
    Edge edge = Edge(edgeId, s, e);
    edge.support = support;
    if (support > maxSupport) maxSupport = support;
    edges.push_back(edge);
    nodesMap[s].neighbours.insert(make_pair(e, edgeId));
    nodesMap[e].neighbours.insert(make_pair(s, edgeId));
    edgeId++;
  }
  sortedEdges.resize(edges.size(), 0);
  binStart.resize(MAX_SUP, 0);
  binSize.resize(MAX_SUP, 0);
  curLength.resize(MAX_SUP, 0);
  cout << "Nodes Cnt:" << nodesMap.size() << endl;
}

/**
 * get the neighbourhood of egde (s, e)
 */

int Graph::getInter(int s, int e) {
  if (nodesMap[s].neighbours.size() > nodesMap[e].neighbours.size()) {
    swap(s, e);
  }

  unordered_map<int, int> &nodeS = nodesMap[s].neighbours;
  unordered_map<int, int> &nodeE = nodesMap[e].neighbours;

  int ret = 0;
  for (unordered_map<int, int>::iterator itr = nodeS.begin();
       itr != nodeS.end(); itr++) {
    if (nodeE.find(itr->first) != nodeE.end()) {
      ret++;
    }
  }
  return ret;
}

/**
 * calc the support for all the edges
 */
void Graph::calcSupport() {
  cout << "calcSupport" << endl;
  int t0 = time(NULL);
  int i = 0;
  for (vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
    i++;
    int s = itr->s;
    int e = itr->e;
    itr->support = getInter(s, e);
    itr->originSupport = itr->support;
    if (itr->support > maxSup) maxSup = itr->support;
  }
  int t1 = time(NULL);
  cout << "Time cost: " << t1 - t0 << endl;
}

/**
 * sort the edges with their support
 */
void Graph::arrangeSortedEdges() {
  cout << "Arrange sorted edges..." << endl;
  for (int i = 0; i < edges.size(); i++) {
    binSize[edges[i].support]++;
  }

  binStart[0] = 0;
  for (int i = 1; i < binStart.size(); i++) {
    binStart[i] = binStart[i - 1] + binSize[i - 1];
  }

  for (vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
    int support = itr->support;
    curLength[support]++;
    sortedEdges[binStart[support] + curLength[support] - 1] = itr->id;
    itr->pos = binStart[support] + curLength[support] - 1;
  }
}

/**
 * support(edgeId) -= 1, rearrange its position in the array
 * @param edgeId
 * @param k
 */
void Graph::decSupportForEdge(int edgeId, int k) {
  if (!edges[edgeId].alive) return;
  if (edges[edgeId].support < k - 2) return;
  int oldSupport = edges[edgeId].support;
  int x = binStart[oldSupport];
  int pos = edges[edgeId].pos;
  int temp = sortedEdges[x];
  sortedEdges[x] = sortedEdges[pos];
  sortedEdges[pos] = temp;
  edges[sortedEdges[x]].pos = x;
  edges[sortedEdges[pos]].pos = pos;
  binStart[oldSupport]++;
  edges[edgeId].support--;
}

/**
 * remove an edge (s, e), delete the support of the connected edges which can
 * form a triangle with it
 * @param s
 * @param e
 * @param k
 */
void Graph::decSupport(int s, int e, int k) {
  if (nodesMap[s].neighbours.size() > nodesMap[e].neighbours.size()) {
    swap(s, e);
  }
  unordered_map<int, int> &neighS = nodesMap[s].neighbours;
  unordered_map<int, int> &neighE = nodesMap[e].neighbours;
  for (unordered_map<int, int>::iterator itr = neighS.begin();
       itr != neighS.end(); itr++) {
    if (neighE.find(itr->first) != neighE.end()) {
      int w = itr->first;
      int edgeId1 = neighS[w];
      int edgeId2 = neighE[w];
      // both of the two edges are alive, otherwise the triangle doesn't exists.
      if (edges[edgeId1].alive && edges[edgeId2].alive) {
        decSupportForEdge(edgeId1, k);
        decSupportForEdge(edgeId2, k);
      }
    }
  }
}

/**
 * get K truss
 * @param k
 * @return
 */
int Graph::getKTruss(int k) {
  unordered_set<pair<int, int>, myHash> removeEdges;
  int round = 0;
  cout << "---------------------------------------" << endl;
  cout << k << "-Truss" << endl;
  int tot = 0;
  while (true) {
    int cnt = 0;
    round++;
    cout << "Round:" << round << endl;
    for (int i = 0; i < edges.size(); i++) {
      Edge edge = edges[sortedEdges[i]];
      if (edge.support >= k - 2) continue;
      if (!edge.alive || edge.deleted) continue;
      edges[sortedEdges[i]].alive = false;
      // annotate the edge
      edges[sortedEdges[i]].trussClass = k - 1;
      decSupport(edge.s, edge.e, k);
      cnt++;
      if (isDirect) {
        // useless now, only for deep decomposition
        oriDelCnt[k] += 1;
      }
    }
    cout << "k:" << k << " Round:" << round << " Cnt:" << cnt << endl;
    if (cnt == 0) break;
    tot += cnt;
  }
  cout << "---------------------------------------" << endl;
  return tot;
}

/**
 * now useless, only for deep truss decomposition
 * @param graphTasks
 * @param k
 * @param maxK
 */
void Graph::floodFill(queue<GraphTask> &graphTasks, int k, int maxK) {
  int aliveCnt = 0;

  vector<GraphTask> ret;
  for (unordered_map<int, Node>::iterator itr = nodesMap.begin();
       itr != nodesMap.end(); itr++) {
    itr->second.visited = false;
    itr->second.color = -1;
  }
  int colorId = 0;
  for (unordered_map<int, Node>::iterator itr = nodesMap.begin();
       itr != nodesMap.end(); itr++) {
    Node &node = itr->second;
    if (!node.visited) {
      node.color = colorId++;
      queue<int> q = queue<int>();
      q.push(itr->first);
      node.visited = true;
      vector<int> compNodes;
      unordered_set<int> compEdges;
      compNodes.push_back(itr->first);
      while (!q.empty()) {
        int id = q.front();
        q.pop();
        unordered_map<int, int> &neighs = nodesMap[id].neighbours;
        int color = nodesMap[id].color;
        for (unordered_map<int, int>::iterator neigh = neighs.begin();
             neigh != neighs.end(); neigh++) {
          int j = neigh->first;
          int edgeId = neigh->second;
          Node &nodeJ = nodesMap[j];
          if (edges[edgeId].alive) compEdges.insert(edgeId);
          if (edges[edgeId].alive && !nodeJ.visited && nodeJ.color < 0) {
            nodeJ.visited = true;
            q.push(j);

            compNodes.push_back(j);
            nodeJ.color = color;
          }
        }
      }

      if (compEdges.size() == 3) {
        inTrussEdges += compEdges.size();
        aliveCnt += compEdges.size();
        if (isDirect && k == 3) {
          finalTrussEdgeCnt[3] += 3;
          oriDelCnt[4] += 3;
        }
      } else if (compEdges.size() > 3 && k == maxK) {
        GraphTask task = GraphTask(FINAL_TRUSS, k, false, isDirect);
        for (unordered_set<int>::iterator itr = compEdges.begin();
             itr != compEdges.end(); itr++) {
          task.edges.push_back(
              RawEdge(edges[*itr].s, edges[*itr].e, edges[*itr].support));
        }
        graphTasks.push(task);
        aliveCnt += compEdges.size();

      } else {
        if (compEdges.size() > 3) {
          this->trussGraphId++;
          cout << "Truss to process:" << compNodes.size() << endl;
          if (compEdges.size() >= 3) aliveCnt += compEdges.size();
          GraphTask task = GraphTask(PROCESS_FURTHER, k + 1, true, isDirect);
          for (unordered_set<int>::iterator itr = compEdges.begin();
               itr != compEdges.end(); itr++) {
            task.edges.push_back(
                RawEdge(edges[*itr].s, edges[*itr].e, edges[*itr].support));
          }
          graphTasks.push(task);
        }
      }
    }
  }

  // can not be decomposed more
  cout << "Alive Cnt:" << aliveCnt << endl;
  if (aliveCnt == 0 && k > 3) {
    GraphTask task = GraphTask(FINAL_TRUSS, k - 1, false, isDirect);
    for (int i = 0; i < edges.size(); i++) {
      task.edges.push_back(RawEdge(edges[i].s, edges[i].e, 0));
    }
    graphTasks.push(task);
    return;
  }

  this->remainGraphId++;
  int newK;
  if (k == 3)
    newK = 2;
  else
    newK = 3;
  GraphTask task = GraphTask(PROCESS_FURTHER, newK, false, false);
  for (int i = 0; i < edges.size(); i++) {
    if (!edges[i].alive) {
      task.edges.push_back(RawEdge(edges[i].s, edges[i].e, 0));
    }
  }
  if (task.edges.size() > 1) {
    graphTasks.push(task);
  } else {
    remainEdges += task.edges.size();
  }
  return;
}

void Graph::outputDenseRemain() {
  this->finalTrussId++;
  string filename =
      "./" + outputPrefix + "/" + "ft" + to_string(this->finalTrussId) + ".txt";
  outputEdgesToFile(filename);
}

void Graph::trussDecomposition(queue<GraphTask> &graphTasks, int k, int maxK) {
  getKTruss(k);
  floodFill(graphTasks, k, maxK);
}

/**
 * used for annotate each edge with trussness (truss number)
 * @param maxK, the largest k
 */
void Graph::annoEdgeWithTrussClass(int maxK) {
  for (int i = 3; i <= maxK; i++) {
    getKTruss(i);
  }
  for (int i = 0; i < edges.size(); i++) {
    if (edges[i].trussClass == -1) {
      edges[i].trussClass = maxK;
    }
  }
}

vector<int> Graph::getSurroudings(int edgeId, int maxK) {
  int s = edges[edgeId].s;
  int e = edges[edgeId].e;
  if (nodesMap[s].neighbours.size() > nodesMap[e].neighbours.size()) {
    swap(s, e);
  }

  unordered_map<int, int> &nodeS = nodesMap[s].neighbours;
  unordered_map<int, int> &nodeE = nodesMap[e].neighbours;

  int curTrussClass = edges[edgeId].trussClass;
  vector<int> ret(maxK + 1, 0);

  for (unordered_map<int, int>::iterator itr = nodeS.begin();
       itr != nodeS.end(); itr++) {
    int w = itr->first;
    if (nodeE.find(itr->first) != nodeE.end()) {
      int edgeId1 = nodeS[w];
      int edgeId2 = nodeE[w];
      int minK = min(curTrussClass,
                     min(edges[edgeId1].trussClass, edges[edgeId2].trussClass));
      for (int i = 3; i <= minK; i++) {
        ret[i]++;
      }
    }
  }
  return ret;
}

void Graph::outputGraphToFile(string filename) {
  ofstream fout(filename.c_str(), ios::app);
  for (int i = 0; i < edges.size(); i++) {
    fout << edges[i].s << " " << edges[i].e << endl;
  }
  fout.close();
}

void Graph::outputEdgesToFile(string filename) {
  ofstream fout(filename.c_str());
  for (int i = 0; i < edges.size(); i++) {
    fout << edges[i].s << " " << edges[i].e << endl;
  }
  fout.close();
}

vector<RemainGraphInfo> Graph::generateConnectedSubgraphs(
    string outputPrefix, int &maxRemainGraphSize) {
  vector<RemainGraphInfo> ret;
  for (unordered_map<int, Node>::iterator itr = nodesMap.begin();
       itr != nodesMap.end(); itr++) {
    Node &node = itr->second;
    node.visited = false;
  }
  for (unordered_map<int, Node>::iterator itr = nodesMap.begin();
       itr != nodesMap.end(); itr++) {
    Node &node = itr->second;
    if (!node.visited && !node.deleted) {
      queue<int> q = queue<int>();
      q.push(itr->first);
      node.visited = true;
      vector<int> compNodes;
      unordered_set<int> compEdges;
      compNodes.push_back(itr->first);
      while (!q.empty()) {
        int id = q.front();
        q.pop();
        unordered_map<int, int> &neighs = nodesMap[id].neighbours;
        for (unordered_map<int, int>::iterator neigh = neighs.begin();
             neigh != neighs.end(); neigh++) {
          int j = neigh->first;
          int edgeId = neigh->second;
          Node &nodeJ = nodesMap[j];
          if (!edges[edgeId].deleted) compEdges.insert(edgeId);
          if (!nodeJ.visited && !nodeJ.deleted) {
            nodeJ.visited = true;
            q.push(j);
            compNodes.push_back(j);
          }
        }
      }
      if (compEdges.size() > 1) {
        this->remainSubgraphId++;
        ofstream fout("./" + outputPrefix + "/remain.txt", ios::app);
        if (true) {
          for (unordered_set<int>::iterator edge = compEdges.begin();
               edge != compEdges.end(); edge++) {
            if (!edges[*edge].deleted) {
              fout << edges[*edge].s << " " << edges[*edge].e << endl;
            }
          }
        }
        fout.close();
        ret.push_back(RemainGraphInfo(compNodes.size(), compEdges.size()));
        if (compNodes.size() > maxRemainGraphSize) {
          maxRemainGraphSize = compNodes.size();
          cout << "Max Remain Size:" << maxRemainGraphSize << endl;
        }
      }
    }
  }
  return ret;
}

void Graph::outputEdgeInfo(string outputPrefix) {
  ofstream fout(outputPrefix + "/EdgeInfo.graph");
  for (int i = 0; i < edges.size(); i++) {
    fout << edges[i].s << " " << edges[i].e << " " << edges[i].originSupport
         << " " << edges[i].trussClass << endl;
  }
  fout.close();
}

void Graph::end() {}

void Graph::test() {}
