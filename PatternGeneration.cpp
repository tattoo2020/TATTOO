#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>

#include "FastPGGraph.h"
#include "Pattern.h"
using namespace std;

const int datasetNum = 10;

vector<string> graphFiles(10);
vector<string> outputPrefix_datanames(10);
// final patterns diversity
vector<vector<Pattern>> finalPatterns(10);
vector<vector<int>> finalPatternsIndex(10);
vector<graph_t> patternSet[32][32];
vector<int> patternCount[32][32];

void prejob4file() {
  graphFiles[0] = "data/loc-brightkite/Brightkite_edges.txt";
  outputPrefix_datanames[0] = "loc-brightkite";

  graphFiles[1] = "data/loc-gowalla/Gowalla_edges.txt";
  outputPrefix_datanames[1] = "loc-gowalla";

  graphFiles[2] = "data/amazon/com-amazon.graph";
  outputPrefix_datanames[2] = "amazon";

  graphFiles[3] = "data/com-youtube/com-youtube.ungraph.txt";
  outputPrefix_datanames[3] = "com-youtube";

  graphFiles[4] = "data/com-dblp/com-dblp.ungraph.txt";
  outputPrefix_datanames[4] = "com-dblp";

  graphFiles[5] = "data/lj/com-lj.ungraph.txt";
  outputPrefix_datanames[5] = "lj";

  graphFiles[6] = "data/skitter/as-skitter.graph";
  outputPrefix_datanames[6] = "skitter";

  graphFiles[7] = "data/roadnet-ca/roadNet-CA.txt";
  outputPrefix_datanames[7] = "roadnet-ca";
  graphFiles[8] = "data/roadnet-tx/roadNet-TX.txt";
  outputPrefix_datanames[8] = "roadnet-tx";
  graphFiles[9] = "data/roadnet-pa/roadNet-PA.txt";
  outputPrefix_datanames[9] = "roadnet-pa";

  FastPGGraph::finalSet4alldataset.clear();
  FastPGGraph::finalSet4alldataset.resize(10);
}

// defaut + truss in final
void test4trussonly() {
  // test for truss-only in final
  for (int j = 0; j < 10; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    //		pgGraph.genQuery(4, 30, 100);							//random gen 100*27
    //query based on dataset
    pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
    ofstream myout9(outputPrefix + outputPrefix_datanames[j] +
                    "_default+truss_results.txt");
    pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
    pgGraph.finalSet.clear();
    assert(30 == pgGraph.patterns.size());
    for (int i = 0; i < pgGraph.patterns.size(); ++i)
      if (pgGraph.patterns[i].type <= 6)  // default, truss-based
        pgGraph.finalSet.emplace_back(i);
    cout << "num=" << pgGraph.finalSet.size() << endl;
    pgGraph.paratest(myout9, 4, 30, 100, pgGraph.finalSet, 0);
    myout9.close();
  }
}

// default + all truss
void test4alltruss() {
  // test for truss-only in final
  for (int j = 0; j < 10; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    //		pgGraph.genQuery(4, 30, 100);							//random gen 100*27
    //query based on dataset
    pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
    ofstream myout9(outputPrefix + outputPrefix_datanames[j] +
                    "_default+alltruss_results.txt");
    pgGraph.patterns.clear();
    pgGraph.genDefaultPatterns();
    pgGraph.trussEdgeCnt.resize(20, 0);
    pgGraph.genTrussPatterns();
    pgGraph.genCombinedPatterns();
    pgGraph.finalSet.clear();
    pgGraph.boostPatterns.clear();

    for (int i = 0; i < pgGraph.patterns.size(); ++i) {
      pgGraph.calCognition3Single(pgGraph.patterns[i]);
      pgGraph.finalSet.emplace_back(i);
      pgGraph.boostPatterns.emplace_back(
          pgGraph.pattern2boostGraph(pgGraph.patterns[i]));
    }
    cout << "num=" << pgGraph.finalSet.size() << endl;
    pgGraph.paratest(myout9, 4, 30, 100, pgGraph.finalSet, 0);
    myout9.close();
  }
}

void test4graphlet() {
  // test for Graphlet. of all 10 dataset
  for (int j = 0; j < 10; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    //		pgGraph.genQuery(4, 30, 100);							//random gen 100*27
    //query based on dataset
    pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
    ofstream myout9(outputPrefix + outputPrefix_datanames[j] +
                    "_graphlet_results9.txt");
    pgGraph.genGraphlet();
    pgGraph.graphletTest(9);  // 2 3 4-nodes
    pgGraph.paratest(myout9, 4, 30, 100, pgGraph.finalSet, 0);
    myout9.close();
    ofstream myout30(outputPrefix + outputPrefix_datanames[j] +
                     "_graphlet_results30.txt");
    pgGraph.graphletTest(30);  // 2 3 4 5-nodes
    pgGraph.paratest(myout30, 4, 30, 100, pgGraph.finalSet, 0);
    myout30.close();
    //		pgGraph.outputSummaryInfo(outputPrefix  + "_SummaryInfo.txt");
  }
}

void test4Random() {
  // test for Random. of all 10 dataset
  for (int j = 0; j < 10; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    // pgGraph.genQuery(4, 30, 3);							//random gen 100*27 query based
    // on dataset pgGraph.genDefaultPatterns(); for (int z = 4, i = 5, j = 0; z
    // < 30; ++z) {
    // pgGraph.patterns.emplace_back(pgGraph.queryPatterns[i][j]);
    //	++i;
    //	if (i > 15) {
    //		i = 7;
    //		++j;
    //	}
    //}
    // for (int z = 0; z < 30; ++z)
    //	pgGraph.outputPngPattern(pgGraph.patterns[z], outputPrefix +
    //"randomPatterns/"); pgGraph.randomTest(30);
    // pgGraph.saveFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
    pgGraph.loadFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
    pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
    // pgGraph.loadFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
    ofstream myout9(outputPrefix + outputPrefix_datanames[j] +
                    "_random_results9.txt");
    pgGraph.randomTest(9);
    pgGraph.paratest(myout9, 4, 30, 100, pgGraph.finalSet, 0);
    myout9.close();
    ofstream myout30(outputPrefix + outputPrefix_datanames[j] +
                     "_random_results30.txt");
    pgGraph.randomTest(30);
    pgGraph.paratest(myout30, 4, 30, 100, pgGraph.finalSet, 0);
    myout30.close();
    //		pgGraph.outputSummaryInfo(outputPrefix  + "_SummaryInfo.txt");
  }
}

void test4RandomDiv() {
  vector<int> othersFinalSet;
  ofstream myout("data/random_div30.txt");
  for (int j = 0; j < 10; ++j) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    pgGraph.loadFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
    double sumdiv1 = 0, sumdiv2 = 0;
    for (int i = 4; i < 30; ++i) {
      othersFinalSet.clear();
      for (int k = 0; k < 30; ++k)
        if (pgGraph.patterns[k].getEdgeCnt() < pgGraph.patterns[i].getEdgeCnt())
          othersFinalSet.emplace_back(k);
      int div1 = pgGraph.minStep4filledQuery(
          othersFinalSet, i, pgGraph.patterns[i].getEdgeCnt(), true);
      double div2 = div1 * 1.0 / pgGraph.patterns[i].getEdgeCnt();
      sumdiv1 += div1;
      sumdiv2 += div2;
    }
    myout << outputPrefix << "\n"
          << sumdiv1 / 26.0 << " " << sumdiv2 / 26.0 << endl;
  }
  myout.close();
}

void test4TattooDiv() {
  vector<int> othersFinalSet;
  ofstream myout("data/tattoo_div30.txt");
  for (int j = 0; j < 10; ++j) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph();  //* read the graph

    pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
    double sumdiv1 = 0, sumdiv2 = 0;
    for (int i = 4; i < 30; ++i) {
      othersFinalSet.clear();
      for (int k = 0; k < 30; ++k)
        if (pgGraph.patterns[k].getEdgeCnt() < pgGraph.patterns[i].getEdgeCnt())
          othersFinalSet.emplace_back(k);
      int div1 = pgGraph.minStep4filledQuery(
          othersFinalSet, i, pgGraph.patterns[i].getEdgeCnt(), true);
      double div2 = div1 * 1.0 / pgGraph.patterns[i].getEdgeCnt();
      sumdiv1 += div1;
      sumdiv2 += div2;
    }
    myout << outputPrefix << "\n"
          << sumdiv1 / 26.0 << " " << sumdiv2 / 26.0 << endl;
  }
  myout.close();
}

void testMaxPatternSize() {
  for (int j = 0; j <= 5; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph(graphFiles[j]);  //* read the graph
    pgGraph.genQuery(4, 30, 100);  // random gen 100*27 query based on dataset
    ofstream myout(outputPrefix + "testResult/" + outputPrefix_datanames[j] +
                   "_results.txt");
    for (int i = 9; i <= 18; i += 3) {
      pgGraph.MAX_PATTERN_SIZE = i;
      pgGraph.init_part1();
      pgGraph.caltrussNumberForEdges();
      pgGraph
          .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
      pgGraph.genCombinedPatterns();
      pgGraph.calcScoreForComPatterns();
      // pgGraph.searchStars();
      pgGraph.searchMultiLinearStars();
      pgGraph.searchLinesAndCirclesAndIndividuals();
      pgGraph.calCoverageAndCognition(outputPrefix);
      pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, 30, j);
      pgGraph.paratest(myout, 4, 30, 100, pgGraph.finalSet, 0);
      pgGraph.outputSummaryInfo(outputPrefix + to_string(i) +
                                "_SummaryInfo.txt");
    }
  }
}

void testStartDegree4star() {
  for (int j = 9; j >= 0; j--) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph(graphFiles[j]);  //* read the graph
    pgGraph.genQuery(4, 30, 100);  // random gen 100*27 query based on dataset
    ofstream myout(outputPrefix + "testResult/" + outputPrefix_datanames[j] +
                   "_results.txt");
    for (int i = 4; i <= 10; i++) {
      pgGraph.init_part1();
      pgGraph.STAR_START_DEGREE = i;
      pgGraph.caltrussNumberForEdges();
      pgGraph
          .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
      pgGraph.genCombinedPatterns();
      pgGraph.calcScoreForComPatterns();
      pgGraph.searchStars();
      pgGraph.searchLinesAndCirclesAndIndividuals();
      pgGraph.calCoverageAndCognition(outputPrefix);
      pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, 30, j);
      // cout << pgGraph.paratest(myout, 4, 30, 100, 0) << "\n\n";
    }
  }
}

void testTokKbyCov() {
  for (int j = 9; j >= 0; j--) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph(graphFiles[j]);  //* read the graph
    pgGraph.genQuery(4, 30, 100);  // random gen 100*27 query based on dataset
    ofstream myout(outputPrefix + "testResult/" + outputPrefix_datanames[j] +
                   "_results.txt");

    pgGraph.caltrussNumberForEdges();
    pgGraph
        .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
    pgGraph.genCombinedPatterns();
    pgGraph.calcScoreForComPatterns();
    pgGraph.searchStars();
    pgGraph.searchLinesAndCirclesAndIndividuals();
    pgGraph.calCoverageAndCognition(outputPrefix);

    pgGraph.predo4test();
    for (int i = 10; i <= 50; i += 10) {
      pgGraph.choosefinalpatternsbycov(i);
      // cout << pgGraph.paratest(myout, 4, 30, 100, 0) << "\n\n";
    }
  }
}

void testComPatterns4allDatasets() {
  for (int j = 9; j >= 0; j--) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph(graphFiles[j]);  //* read the graph
    ofstream myout(outputPrefix + "testResult/" + outputPrefix_datanames[j] +
                   "_results.txt");
    pgGraph.caltrussNumberForEdges();
    pgGraph
        .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
    pgGraph.genCombinedPatterns();
    pgGraph.calcScoreForComPatterns();
    pgGraph.searchStars();
    pgGraph.searchLinesAndCirclesAndIndividuals();
    pgGraph.calCoverageAndCognition(outputPrefix);
    pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, 50, j);
  }
  FastPGGraph::getcomset4alldataset();
}

void testFinalPatternSetSize() {
  for (int j = 0; j < 10; j++) {
    cout << outputPrefix_datanames[j] +
                "*************************************************\n\n\n\n";
    string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
    FastPGGraph pgGraph = FastPGGraph(graphFiles[j]);  //* read the graph
    pgGraph.genQuery(4, 30, 100);  // random gen 100*27 query based on dataset
    ofstream myout(outputPrefix + "testResult/" + outputPrefix_datanames[j] +
                   "_results.txt");

    pgGraph.caltrussNumberForEdges();
    pgGraph
        .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
    pgGraph.genCombinedPatterns();
    pgGraph.calcScoreForComPatterns();
    // pgGraph.searchStars();
    pgGraph.searchMultiLinearStars();
    pgGraph.searchLinesAndCirclesAndIndividuals();
    pgGraph.calCoverageAndCognition(outputPrefix);
    pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, 50,
                                j);  // gen 50 first

    pgGraph.prejob4finalPatternSetSize();
    pgGraph.Test4finalPatternSetSize(myout, 5);  // test top i
    for (int i = 10; i <= 40; i += 10) {
      pgGraph.Test4finalPatternSetSize(myout, i);  // test top i
    }
  }
}

/**
 * Count and print occurrences of patterns among all datasets.
 */
void testUniquePatterns() {
  cout << "starting unique patterns analysis...\n";
  for (int i = 0; i < datasetNum; ++i) {
    string outputPrefix = "data/" + outputPrefix_datanames[i] + "/";
    FastPGGraph pgGraph = FastPGGraph();
    pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
    FastPGGraph::finalSet4alldataset[i] = pgGraph.patterns;
  }
  for (int i = 0, j = 0; i < datasetNum; ++i)
    for (auto &p : FastPGGraph::finalSet4alldataset[i]) {
      auto g = FastPGGraph::pattern2boostGraph(p);
      bool found = false;
      int V = p.getNodesCnt();
      int E = p.getEdgeCnt();
      std::vector<graph_t::vertex_descriptor> f(V);
      for (int j = 0; j < patternSet[V][E].size(); ++j)
        if (isomorphism(
                g, patternSet[V][E][j],
                boost::isomorphism_map(boost::make_iterator_property_map(
                    f.begin(), boost::get(boost::vertex_index, g))))) {
          found = true;
          finalPatternsIndex[i].emplace_back(j);
          ++patternCount[V][E][j];
          break;
        }
      if (!found) {
        finalPatternsIndex[i].emplace_back(patternSet[V][E].size());
        patternSet[V][E].emplace_back(g);
        patternCount[V][E].emplace_back(1);
      }
    }
  cout << "printing unique patterns...\n";
  for (int i = 0; i < datasetNum; ++i)
    for (int j = 0; j < FastPGGraph::finalSet4alldataset[i].size(); ++j) {
      auto &p = FastPGGraph::finalSet4alldataset[i][j];
      int V = p.getNodesCnt();
      int E = p.getEdgeCnt();
      if (patternCount[V][E][finalPatternsIndex[i][j]] == 1)
        FastPGGraph::outputFinalPngPattern(
            p, "data/" + outputPrefix_datanames[i] + "/uniquePatterns/",
            "#" + to_string(j + 1));
      FastPGGraph::outputFinalPngPattern(
          p,
          "commonPatterns/" +
              to_string(patternCount[V][E][finalPatternsIndex[i][j]]) + "/",
          outputPrefix_datanames[i] + "=" + to_string(j + 1));
    }
}

/*datasetid maps in  function prejob4file(). */

void runSingleDataset(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph

  //	if (!pgGraph.LOADED)
  //	{
  //		pgGraph.genQuery(4, 30, 100);			//random gen 100*27
  //query based on dataset 		pgGraph.saveQueries(outputPrefix +
  //"random_query4-30-100.txt", 4, 30, 100);
  pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
  //	}
  //	else
  //	{
  //		pgGraph.loadQuerys_steps(outputPrefix + "query_steps.txt");
  //	}

  //	if (!pgGraph.LOADED)
  //	{
  pgGraph.caltrussNumberForEdges();
  pgGraph.genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
  pgGraph.genCombinedPatterns();
  pgGraph.outputGin4Gephi(outputPrefix + "gin.dl");
  pgGraph.outputGout4Gephi(outputPrefix + "gout.dl");
  pgGraph.calcScoreForComPatterns();
  pgGraph.searchMultiLinearStars();
  pgGraph.outputGr4Gephi(outputPrefix + "gr.dl");
  //		return;
  pgGraph.searchLinesAndCirclesAndIndividuals();
  pgGraph.calCoverageAndCognition(outputPrefix);
  pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns30/", 0.5, 30,
                              datasetid);
  //	}
  //	else
  //	{
  //		pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns.txt");
  //	}

  //		return;
  // if you want to test:
  ofstream myout(outputPrefix + outputPrefix_datanames[datasetid] +
                 "_results30.txt");

  //	pgGraph.testByPatternClass(myout, 4, 30, 100);
  cout << pgGraph.paratest(myout, 4, 30, 100, pgGraph.finalSet, 0) << "\n\n";
  // final patterns analysis
  pgGraph.outputSummaryInfo(outputPrefix + outputPrefix_datanames[datasetid] +
                            "_SummaryInfo30.txt");

  pgGraph.saveFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
  for (auto pid : pgGraph.finalSet)
    finalPatterns[datasetid].emplace_back(pgGraph.patterns[pid]);
}

void runSingleDataset9(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph

  //	if (!pgGraph.LOADED)
  //	{
  //		pgGraph.genQuery(4, 30, 100);			//random gen 100*27
  //query based on dataset 		pgGraph.saveQueries(outputPrefix +
  //"random_query4-30-100.txt", 4, 30, 100);
  pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
  //	}
  //	else
  //	{
  //		pgGraph.loadQuerys_steps(outputPrefix + "query_steps.txt");
  //	}

  //	if (!pgGraph.LOADED)
  //	{
  pgGraph.caltrussNumberForEdges();
  pgGraph.genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
  pgGraph.genCombinedPatterns();
  pgGraph.calcScoreForComPatterns();
  pgGraph.searchMultiLinearStars();
  //		return;
  pgGraph.searchLinesAndCirclesAndIndividuals();
  pgGraph.calCoverageAndCognition(outputPrefix);
  pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns9/", 0.5, 9,
                              datasetid);
  //	}
  //	else
  //	{
  //		pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns.txt");
  //	}

  //		return;
  // if you want to test:
  ofstream myout(outputPrefix + outputPrefix_datanames[datasetid] +
                 "_results9.txt");

  //	pgGraph.testByPatternClass(myout, 4, 30, 100);
  cout << pgGraph.paratest(myout, 4, 30, 100, pgGraph.finalSet, 0) << "\n\n";
  // final patterns analysis
  pgGraph.outputSummaryInfo(outputPrefix + outputPrefix_datanames[datasetid] +
                            "_SummaryInfo9.txt");

  pgGraph.saveFinalpatterns(outputPrefix + "finalPatterns9.txt", 9);
  for (auto pid : pgGraph.finalSet)
    finalPatterns[datasetid].emplace_back(pgGraph.patterns[pid]);
}

// generate final results and random results
void genSingleDatasetResult(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string dataPrefix = "data_13092020-cov2-cog7/";
  string outputPrefix = dataPrefix + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph

  pgGraph.caltrussNumberForEdges();
  pgGraph.genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
  pgGraph.genCombinedPatterns();
  pgGraph.calcScoreForComPatterns();
  pgGraph.searchMultiLinearStars();
  pgGraph.searchLinesAndCirclesAndIndividuals();
  pgGraph.calCoverageAndCognition(outputPrefix);
  int num = 30;
  pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, num,
                              datasetid);

  int cnt = 0;
  for (auto pid : pgGraph.finalSet)
    if (cnt < num && pgGraph.patterns[pid].getEdgeCnt() >= 4 &&
        pgGraph.patterns[pid].getEdgeCnt() <= 30) {  // size 4-30
      finalPatterns[datasetid].emplace_back(pgGraph.patterns[pid]);
      ++cnt;
    }
  /*vector<int> randomPatternsNumPerSize(16, 0);
  vector<Pattern> randomPatterns;
  pgGraph.queryPatterns.clear();
  pgGraph.queryPatterns.resize(16);
  while (cnt > 0) {
          for (int i = 4; i <= 15 && cnt > 0; ++i) {
                  ++randomPatternsNumPerSize[i];
                  --cnt;
          }
  }
  for (int i = 4; i <= 15; ++i)
          pgGraph.genQueryRandomly("", randomPatternsNumPerSize[i], i);
  for (int i = 4; i <= 15; ++i)
          for (auto q : pgGraph.queryPatterns[i])
                  randomPatterns.emplace_back(q);*/

  pgGraph.outputLocalPatterns4UI(finalPatterns[datasetid],
                                 outputPrefix + "finalPatterns4UI.txt");
  pgGraph.outputSummaryInfo(
      outputPrefix + outputPrefix_datanames[datasetid] + "_SummaryInfo30.txt",
      outputPrefix_datanames[datasetid], false);
  pgGraph.outputSummaryInfo(dataPrefix + "SummaryInfo30.txt",
                            outputPrefix_datanames[datasetid], true);
  pgGraph.saveFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
  // pgGraph.outputLocalPatterns4UI(randomPatterns, outputPrefix +
  // "randomPatterns4UI.txt"); FastPGGraph pgGraph = FastPGGraph();
  // pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns30.txt", 30);
  // pgGraph.outputLocalPatterns4UI(pgGraph.patterns, outputPrefix +
  // "finalPatterns4UI_new.txt");
}

void testCoverageLowerBound(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph

  if (true) {
    pgGraph.caltrussNumberForEdges();
    pgGraph
        .genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
    pgGraph.genCombinedPatterns();
    pgGraph.calcScoreForComPatterns();
    pgGraph.searchMultiLinearStars();
    pgGraph.searchLinesAndCirclesAndIndividuals();

    // before subgraph check and add frequency
    // outsideCoverage.resize(pgGraph.patterns.size(), 0);
    for (int i = 0; i < pgGraph.patterns.size(); ++i) {
      auto &p = pgGraph.patterns[i];
      if (p.type > PatType::OutsideCombinedStar)
        p.outsideCoverage = p.score * p.getEdgeCnt();
    }

    pgGraph.calCoverageAndCognition(outputPrefix);
    pgGraph.decideFinalPatterns(outputPrefix + "finalPatterns/", 0.5, 30,
                                datasetid);
    // pgGraph.decideFinalPatternsByCov(30);
    vector<Pattern> pats;
    for (int i : pgGraph.finalSet) pats.emplace_back(pgGraph.patterns[i]);
    pgGraph.outputLocalPatterns4UI(pats, outputPrefix +
                                             outputPrefix_datanames[datasetid] +
                                             "_cov_only_patterns.txt");
  } else {
    pgGraph.caltrussNumberForEdges();
    pgGraph.loadFinalpatterns(outputPrefix + "finalPatterns.txt", 30);
  }

  ofstream fout(outputPrefix + outputPrefix_datanames[datasetid] +
                "chord_out_cov.txt");
  int e = pgGraph.NUM_E, v = pgGraph.NUM_V;
  int x = 0, y = 0, cx = 0, cy = 0, minTruss = 1000, minStar = e;
  for (int i = 4; i < pgGraph.finalSet.size(); ++i) {
    int id = pgGraph.finalSet[i];
    cout << "id=" << id << " p_size=" << pgGraph.patterns.size() << endl;
    auto &p = pgGraph.patterns[id];
    if (p.type <= PatType::TrussCombined2) {
      if (p.type == PatType::Truss) {
        int k = (p.getEdgeCnt() + 3) / 2;
        if (k < minTruss) minTruss = k;
      }
      ++cx;
    } else {
      if (p.type <= PatType::OutsideCombinedStar) {
        if (p.type == PatType::OutsideStar && p.getEdgeCnt() < minStar)
          minStar = p.getEdgeCnt();
      } else
        y += p.outsideCoverage;
      ++cy;
    }
  }
  int gin = 0, gout = 0;
  for (int i = 0; i < e; ++i)
    if (pgGraph.trussClass[i] >= 3) {
      ++gin;  // G_in
      if (minTruss <= pgGraph.trussClass[i]) ++x;
    } else {
      if (pgGraph.GoutDegree[pgGraph.edges[i].s] >= minStar ||
          pgGraph.GoutDegree[pgGraph.edges[i].e] >= minStar)
        ++y;
      ++gout;
    }
  fout << gin << " " << gout << endl;
  fout << cx << " " << cy << endl;
  fout << x << " " << y << endl;
  fout.close();
}

void runChordAndTrussComparison(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph
  // pgGraph.loadQuerys_steps(outputPrefix + "query_steps.txt");

  ofstream myout(outputPrefix + outputPrefix_datanames[datasetid] +
                 "_chord_truss.txt");

  pgGraph.caltrussNumberForEdges();
  pgGraph.genQuery(4, 30, 100);  // random gen 100*27 query based on dataset
  // for (int i = 4; i <= 30; ++i)
  // 	for (int j = 0; j < 10; ++j)
  // 		pgGraph.outputPngPattern(pgGraph.queryPatterns[i][j], outputPrefix
  // + "QueryPatterns/");
  pgGraph.genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
  pgGraph.genCombinedPatterns();
  pgGraph.genDefaultPatterns();
  pgGraph.calCognition3("");
  vector<int> ids;
  for (int i = 0; i < pgGraph.patterns.size(); ++i) ids.emplace_back(i);
  for (auto p : pgGraph.patterns)
    pgGraph.outputPngPattern(p, outputPrefix + "kChordPatterns/");
  pgGraph.boostPatterns.resize(pgGraph.patterns.size());
  for (int j = 0; j < pgGraph.patterns.size(); ++j)
    pgGraph.boostPatterns[j] = pgGraph.pattern2boostGraph(pgGraph.patterns[j]);
  pgGraph.paratest(myout, 4, 30, 100, ids, true);

  pgGraph.patterns.clear();
  pgGraph.genDefaultPatterns();
  pgGraph.searchKTrussPatterns(15);
  pgGraph.calCognition3("");
  ids.clear();
  for (int i = 0; i < pgGraph.patterns.size(); ++i) ids.emplace_back(i);
  for (auto p : pgGraph.patterns)
    pgGraph.outputPngPattern(p, outputPrefix + "kTrussPatterns/");
  pgGraph.boostPatterns.resize(pgGraph.patterns.size());
  for (int j = 0; j < pgGraph.patterns.size(); ++j)
    pgGraph.boostPatterns[j] = pgGraph.pattern2boostGraph(pgGraph.patterns[j]);
  pgGraph.paratest(myout, 4, 30, 100, ids, true);

  myout.close();
}

void runUserCogComparison(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[datasetid]);  //* read the graph
  pgGraph.caltrussNumberForEdges();
  pgGraph.genTrussPatterns();  //* generate k(4--9) edge (k-2 triangle) pattern
  pgGraph.genCombinedPatterns();
  pgGraph.calcScoreForComPatterns();
  pgGraph.searchMultiLinearStars();
  pgGraph.searchLinesAndCirclesAndIndividuals();
  pgGraph.calCoverageAndCognition(outputPrefix);
  pgGraph.genQuery(4, 30, 1);
  sort(pgGraph.patterns.begin(), pgGraph.patterns.end(),
       [](Pattern a, Pattern b) -> bool {
         return a.scoreCognition < b.scoreCognition;
       });
  int size = pgGraph.patterns.size();
  int pick_num = 5;
  double prop[5] = {1.0 / 15, 1.0 / 8, 1.0 / 5, 2.0 / 5, 4.5 / 5};
  bool found_chord = false;
  for (int i = 0; i < pick_num; ++i) {
    int j = floor(size * prop[i]);
    if (!found_chord && pgGraph.patterns[j].type > PatType::TrussCombined2)
      for (int j1 = j - 10; j1 <= j + 10; ++j1) {
        if (j1 >= 0 && j1 < size &&
            pgGraph.patterns[j1].type <= PatType::TrussCombined2) {
          j = j1;
          break;
        }
      }
    cout << "j=" << j << endl;
    auto &p = pgGraph.patterns[j];
    if (p.type <= PatType::TrussCombined2) found_chord = true;
    auto sc = p.scoreCognition;
    pgGraph.calCognition4Single(p);
    auto sc_alter = p.scoreCognition;
    string yn = "";
    for (int q = 25; q <= 30; ++q)
      yn += (pgGraph.isSubGraph(p, pgGraph.queryPatterns[q][0]) ? "y" : "n");
    pgGraph.outputFinalPngPattern(
        p, "cogComparison/" + outputPrefix_datanames[datasetid] + "/",
        "pattern#" + to_string(i) + "=" + to_string(sc) + "=" +
            to_string(sc_alter) + "=" + yn);
  }
  // clique
  srand((unsigned)time(NULL));
  int k = rand() % 3 + 5;
  Pattern p(k);
  for (int i = 0; i + 1 < k; ++i)
    for (int j = i + 1; j < k; ++j)
      if (rand() % 100 < 90) p.addEdge(i, j);
  cout << "d=" << k * (k - 1) / 2 - p.getEdgeCnt() << endl;
  pgGraph.calCognition3Single(p);
  auto sc = p.scoreCognition;
  pgGraph.calCognition4Single(p);
  auto sc_alter = p.scoreCognition;
  string yn = "";
  for (int q = 25; q <= 30; ++q)
    yn += (pgGraph.isSubGraph(p, pgGraph.queryPatterns[q][0]) ? "y" : "n");
  pgGraph.outputFinalPngPattern(
      p, "cogComparison/" + outputPrefix_datanames[datasetid] + "/",
      "pattern#" + to_string(5) + "=" + to_string(sc) + "=" +
          to_string(sc_alter) + "=" + yn);
  for (int i = 25; i <= 30; ++i)
    pgGraph.outputFinalPngPattern(
        pgGraph.queryPatterns[i][0],
        "cogComparison/" + outputPrefix_datanames[datasetid] + "/",
        "query#" + to_string(pgGraph.queryPatterns[i][0].getEdgeCnt()));
}

void testChordAndTrussComparison() {
  for (int i = 5; i < 10; ++i) runChordAndTrussComparison(i);
}

void testUserCogComparison() {
  runUserCogComparison(0);
  runUserCogComparison(2);
  runUserCogComparison(7);
}

void genRandomPatterns4UI(int datasetid) {
  cout << outputPrefix_datanames[datasetid] +
              "*************************************************\n\n\n\n";
  string outputPrefix = "data/" + outputPrefix_datanames[datasetid] + "/";
  FastPGGraph pgGraph = FastPGGraph();  //* read the graph
  pgGraph.loadFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
  pgGraph.outputLocalPatterns4UI(pgGraph.patterns,
                                 outputPrefix + "randomPatterns4UI.txt");
}

void testExpt3Expt5_tattoo_graphlets() {
  string newDataFolder = "data_23082020-cov2-cog6";
  // string newDataFolder = "data_13092020-cov2-cog7";
  ofstream expt3_fout("data/cog6-expt3.txt");
  ofstream expt5_fout("data/cog6-expt5.txt");
  vector<int> finalSet;
  for (int j = 0; j < 10; j++)
    if (j != 0 && j != 4) {
      cout << outputPrefix_datanames[j] +
                  "*************************************************\n\n\n\n";
      string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
      string newOutputPrefix =
          newDataFolder + "/" + outputPrefix_datanames[j] + "/";
      FastPGGraph pgGraph = FastPGGraph();
      // pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30,
      // 100); assert(pgGraph.queryPatterns[4].size() == 100 &&
      // pgGraph.queryPatterns[30].size() == 100);
      pgGraph.loadQueries2(outputPrefix + "query2.txt",
                           pgGraph.query2Setting[7]);
      assert(pgGraph.query2Patterns.size() == pgGraph.query2Setting[7]);
      pgGraph.loadFinalpatterns(newOutputPrefix + "finalPatterns30.txt", 30);
      assert(pgGraph.patterns.size() == 30);
      assert(pgGraph.boostPatterns.size() == 30);

      expt3_fout << "==============================="
                 << outputPrefix_datanames[j]
                 << "===============================\n";
      // default
      expt3_fout << "======default======" << endl;
      finalSet.clear();
      for (int i = 0; i < 4; ++i) finalSet.emplace_back(i);
      // pgGraph.paratest(expt3_fout, 4, 30, 100, finalSet, 0);
      pgGraph.paratest2(expt3_fout, finalSet);
      // default + chord
      expt3_fout << "======default+chord======" << endl;
      finalSet.clear();
      for (int i = 0; i < pgGraph.patterns.size(); ++i)
        if (pgGraph.patterns[i].type <= PatType::TrussCombined2)
          finalSet.emplace_back(i);
      expt3_fout << "size=" << (int)finalSet.size() << "\n";
      // pgGraph.paratest(expt3_fout, 4, 30, 100, finalSet, 0);
      pgGraph.paratest2(expt3_fout, finalSet);
      /*
      // default + chord + non-chord (TATTOO)
      expt3_fout << "======default+chord+non-chord (TATTOO)======\n";
      finalSet.clear();
      for (int i = 0; i < 30; ++i)
              finalSet.emplace_back(i);
      pgGraph.paratest(expt3_fout, 4, 30, 100, finalSet, 0);*/

      expt5_fout << "==============================="
                 << outputPrefix_datanames[j]
                 << "===============================\n";
      // TATTOO
      expt5_fout << "======TATTOO======" << endl;
      finalSet.clear();
      for (int i = 0; i < 30; ++i) finalSet.emplace_back(i);
      // pgGraph.paratest(expt5_fout, 4, 30, 100, finalSet, 0);
      pgGraph.paratest2(expt5_fout, finalSet);
      // graphlets
      expt5_fout << "======graphlets======" << endl;
      pgGraph.patterns.clear();
      pgGraph.genGraphlet();
      assert(pgGraph.patterns.size() == 30);
      pgGraph.boostPatterns.clear();
      for (auto &pattern : pgGraph.patterns) {
        pgGraph.calCognition6Single(pattern);
        pgGraph.boostPatterns.emplace_back(pgGraph.pattern2boostGraph(pattern));
      }
      finalSet.clear();
      for (int i = 0; i < 30; ++i) finalSet.emplace_back(i);
      // pgGraph.paratest(expt5_fout, 4, 30, 100, finalSet, 0);
      pgGraph.paratest2(expt5_fout, finalSet);
    }
  expt3_fout.close();
  expt5_fout.close();
}

void testExpt4_random() {
  ofstream expt4_fout("data/cog6-expt4.txt" /*, ios_base::app*/);
  vector<int> finalSet;
  for (int j = 0; j < 10; j++)
    if (j != 0 && j != 4) {
      cout << outputPrefix_datanames[j] +
                  "*************************************************\n\n\n\n";
      string outputPrefix = "data/" + outputPrefix_datanames[j] + "/";
      FastPGGraph pgGraph = FastPGGraph();
      // pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30,
      // 100); assert(pgGraph.queryPatterns[4].size() == 100 &&
      // pgGraph.queryPatterns[30].size() == 100);
      pgGraph.loadQueries2(outputPrefix + "query2.txt",
                           pgGraph.query2Setting[7]);
      cout << "query2 size=" << pgGraph.query2Patterns.size() << endl;
      assert(pgGraph.query2Patterns.size() == pgGraph.query2Setting[7]);
      pgGraph.loadFinalpatterns(outputPrefix + "random_patterns30.txt", 30);
      assert(pgGraph.patterns.size() == 30);
      assert(pgGraph.boostPatterns.size() == 30);

      expt4_fout << "==============================="
                 << outputPrefix_datanames[j]
                 << "===============================\n";
      // random 5
      expt4_fout << "======random 5======" << endl;
      finalSet.clear();
      for (int i = 0; i < 5; ++i) finalSet.emplace_back(i);
      pgGraph.paratest2(expt4_fout, finalSet);
      // random 15
      expt4_fout << "======random 15======" << endl;
      finalSet.clear();
      for (int i = 0; i < 15; ++i) finalSet.emplace_back(i);
      pgGraph.paratest2(expt4_fout, finalSet);
      // random 15
      expt4_fout << "======random 30======" << endl;
      finalSet.clear();
      for (int i = 0; i < 30; ++i) finalSet.emplace_back(i);
      pgGraph.paratest2(expt4_fout, finalSet);
    }
  expt4_fout.close();
}

void testAmazonCatapult() {
  // string newDataFolder = "data_23082020-cov2-cog6";
  string newDataFolder = "data_13092020-cov2-cog7";
  ofstream expt5_fout("amazon_catapult_comparison.txt");
  vector<int> finalSet;

  string outputPrefix = "data/" + outputPrefix_datanames[2] + "/";
  string newOutputPrefix =
      newDataFolder + "/" + outputPrefix_datanames[2] + "/";
  FastPGGraph pgGraph = FastPGGraph();
  // pgGraph.loadQueries(outputPrefix + "random_query4-30-100.txt", 4, 30, 100);
  // assert(pgGraph.queryPatterns[4].size() == 100 &&
  // pgGraph.queryPatterns[30].size() == 100);
  pgGraph.loadQueries2(outputPrefix + "query2.txt", pgGraph.query2Setting[7]);
  assert(pgGraph.query2Patterns.size() == pgGraph.query2Setting[7]);
  pgGraph.loadFinalpatterns(newOutputPrefix + "finalPatterns30.txt", 30);
  assert(pgGraph.patterns.size() == 30);
  assert(pgGraph.boostPatterns.size() == 30);
  for (auto &pattern : pgGraph.patterns) pgGraph.calCognition6Single(pattern);

  expt5_fout << "===============================" << outputPrefix_datanames[2]
             << "===============================\n";
  // TATTOO
  expt5_fout << "======TATTOO======" << endl;
  finalSet.clear();
  for (int i = 0; i < 30; ++i) finalSet.emplace_back(i);
  // pgGraph.paratest(expt5_fout, 4, 30, 100, finalSet, 0, true);
  pgGraph.paratest2(expt5_fout, finalSet);
  // CATAPULT
  expt5_fout << "======CATAPULT======" << endl;
  pgGraph.patterns.clear();

  {
    ifstream fin("Amazon_Catapult_withDefaults.txt");
    char ch;
    int a, b, c;
    fin >> ch;
    while (!fin.eof()) {
      if (ch == 't') {
        fin >> ch;
        assert(ch == '#');
        fin >> a >> b;
        assert(pgGraph.patterns.size() == a);
        Pattern p = Pattern(b);
        for (int i = 0; i < p.nodesCnt; ++i) {
          fin >> ch;
          assert(ch == 'v');
          fin >> a >> b;
          assert(i == a && 0 == b);
        }
        while (!fin.eof()) {
          fin >> ch;
          if (ch == 't') break;
          assert(ch == 'e');
          fin >> a >> b >> c;
          assert(a < p.nodesCnt && b < p.nodesCnt && 0 == c);
          p.addEdge(a, b);
        }
        pgGraph.patterns.emplace_back(p);
      } else
        assert(false);
    }
    fin.close();
  }

  assert(pgGraph.patterns.size() == 28);
  pgGraph.boostPatterns.clear();
  for (auto &pattern : pgGraph.patterns) {
    pgGraph.calCognition6Single(pattern);
    pgGraph.boostPatterns.emplace_back(pgGraph.pattern2boostGraph(pattern));
  }
  finalSet.clear();
  for (int i = 0; i < 28; ++i) finalSet.emplace_back(i);
  // pgGraph.paratest(expt5_fout, 4, 30, 100, finalSet, 0, true);
  pgGraph.paratest2(expt5_fout, finalSet);

  expt5_fout.close();
}

void genQuery2(int dataset_id) {
  cout << "==============================="
       << outputPrefix_datanames[dataset_id]
       << "===============================\n";
  string outputPrefix = "data/" + outputPrefix_datanames[dataset_id] + "/";
  FastPGGraph pgGraph = FastPGGraph(graphFiles[dataset_id]);
  pgGraph.caltrussNumberForEdges();
  pgGraph.genTrussPatterns();
  pgGraph.genCombinedPatterns();
  pgGraph.calcScoreForComPatterns();
  pgGraph.searchMultiLinearStars();
  pgGraph.searchLinesAndCirclesAndIndividuals();
  /*printf("circleRecord size=%d\n", pgGraph.circleRecord.size());
  for (const auto& r : pgGraph.circleRecord)
          printf("%d,", r.size());
  printf("\n");
  return;*/
  pgGraph.genQueries2(outputPrefix + "random_query4-30-100.txt",
                      outputPrefix + "query2.txt");
  /*int group_pre_cnt = 0;
  for (int i = 0, group_id = 0; i < pgGraph.query2Patterns.size(); ++i) {
          while (i - group_pre_cnt == pgGraph.query2Setting[group_id]) {
                  group_pre_cnt = i;
                  ++group_id;
          }
          if (i % 50 == 0) printf("i=%d, group_id=%d\n", i, group_id);
          pgGraph.outputFinalPngPattern(pgGraph.query2Patterns[i], outputPrefix
  + "query2_png/", to_string(i) + "-" + pgGraph.query2SettingName[group_id]);
  }*/
}

int main(int argc, char *argv[]) {
  prejob4file();
  // for (int i = 2; i <= 3; ++i)
  //	genRandomPatterns4UI(i);

  /*
  you can choose one of your dataset to run,(dataset id)
    or
  when you need to test para, choose function you needed to run.
  */
  //	for (int i = 0; i < 10; ++i)
  //		runSingleDataset9(i);
  // genSingleDatasetResult(2); // amazon
  // genSingleDatasetResult(3); // youtube
  // genSingleDatasetResult(6); // skitter
  // genSingleDatasetResult(7); // roadnet-ca
  /*for (int i = 0; i < 10; ++i)
          if (i != 0 && i != 4)
                  genQuery2(i);*/
  // testAmazonCatapult();
  // for (int i = 0; i < 10; ++i)
  //	genSingleDatasetResult(i);
  testExpt3Expt5_tattoo_graphlets();
  testExpt4_random();
  //	test4trussonly();
  //	test4alltruss();
  //	testUniquePatterns();
  // testMaxPatternSize();
  //	test4TattooDiv();
  //	test4graphlet();
  // testComPatterns4allDatasets();
  // testFinalPatternSetSize();
  // testStartDegree4star();
  // testTokKbyCov();
  /*unsigned concurentThreadsSupported = boost::thread::hardware_concurrency();
  cout << concurentThreadsSupported << endl;     to see the kenel of cpu*/
  // testChordAndTrussComparison();
  // testUserCogComparison();
  system("pause");
  return 0;
}
