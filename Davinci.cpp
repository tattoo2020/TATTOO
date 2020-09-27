#include <fstream>
#include <iostream>

#include "FastPGGraph.h"

using namespace std;

string datastore;
int num, minSize, maxSize;

void showUsage(char *name) {
  std::cerr
      << "Usage:\n\n"
      << "step1: extract patterns and caluculate scores\n"
      << "\t" << name << " <datastore name> -step1\n\n"
      << "step2: pick patterns by pattern number, min size and max size\n"
      << "\t" << name
      << " <datastore name> -step2 <pattern number> <min size> <max size>\n";
}

bool checkGraphFile(string filename) {
  ifstream fin(filename);
  int n, m, s, e;
  try {
    if (!fin.good()) throw 0;
    fin >> n >> m;
    if (fin.fail() || n <= 0 || m <= 0) throw 1;
    while (fin >> s >> e) {
      if (fin.fail()) throw 2;
      --m;
    }
    //        if (m != 0) throw 3;
  } catch (int e) {
    return false;
  }
  fin.close();
  return true;
}

void outputProgress(string filename, int p) {
  ofstream fout(filename);
  fout << p << endl;
  fout.close();
}
/*the main for GUI*/
int main1(int argc, char *argv[]) {
  if (argc < 2) {
    showUsage(argv[0]);
    exit(1);
  }

  srand((unsigned)time(NULL));
  datastore.assign(argv[1]);

  FastPGGraph::finalSet4alldataset.clear();
  FastPGGraph::finalSet4alldataset.resize(10);

  if (strcmp("-step1", argv[2]) == 0) {
    if (!checkGraphFile("data/" + datastore + "/graph")) {
      ofstream fout("data/" + datastore + "/error");
      fout << "Invalid graph file\n";
      fout.close();
      printf("Invalid graph file!!\n");
      exit(0);
    }
    FastPGGraph pgGraph = FastPGGraph("data/" + datastore + "/graph");
    outputProgress("data/" + datastore + "/progress", 20);
    pgGraph.caltrussNumberForEdges();
    outputProgress("data/" + datastore + "/progress", 40);
    pgGraph
        .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
    pgGraph.genCombinedPatterns();
    pgGraph.calcScoreForComPatterns();
    outputProgress("data/" + datastore + "/progress", 70);
    pgGraph.searchMultiLinearStars();
    pgGraph.searchLinesAndCirclesAndIndividuals();
    outputProgress("data/" + datastore + "/progress", 80);
    pgGraph.calCoverageAndCognition("");
    outputProgress("data/" + datastore + "/progress", 90);
    pgGraph.outputStep1("data/" + datastore);
    outputProgress("data/" + datastore + "/progress", 100);
  } else if (strcmp("-step2", argv[2]) == 0) {
    num = stoi(argv[3]);
    minSize = stoi(argv[4]);
    maxSize = stoi(argv[5]);
    FastPGGraph pgGraph = FastPGGraph();
    pgGraph.inputStep1("data/" + datastore, minSize, maxSize);
    pgGraph.decideFinalPatterns("data/" + datastore + "/", 0.5, num + 10, 0);
    pgGraph.outputStep2("data/" + datastore, num, minSize, maxSize);
  }
  return 0;
}