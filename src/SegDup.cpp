/*
 * SegDup.cpp
 *
 *  Created on: 2 Aug 2022
 *      Author: mac
 */

#include <sys/types.h>
#include <algorithm>
#include <chrono>
#include <cmath>	// for exp
#include <cstdlib>
#include <cstring> // for strcmp
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../utility/debugging.h"
#include "../utility/myrandom.h"
//#include "../utility/niceties.h"

#include "Contender.h"
#include "CophyMap.h"
#include "CophyMultiMap.h"
#include "DupMove.h"
#include "EventCount.h"
#include "EmptyMove.h"
#include "Node.h"
#include "NodeMap.h"
#include "SingleNodeMove.h"
#include "SingleVertexMove.h"
#include "Tree.h"

/**
 * Input to an instance includes whether each node is mapped to a D or S
 */

using namespace std;
using namespace segdup;

bool _debugging(true);
bool _cacheEventCounts(false);
bool _silent(false);
bool _outputProbabilities(false);
bool _saveFinal(false);
bool _saveTrace(false);
bool _saveSampledDistribution(false);
bool _verbose(false);
int outputInterval(100);
int nSteps(10000);
double Tinitial(10.0);
double Tfinal(0.0);
//double SATempSpread(1.0);
//double SATempDecay(1.0);
int nFinal(0);

extern std::mt19937 generator;
extern unsigned seed;

std::ofstream summaryfile;

double lossCost(defLossCost);						// XXX
double duplicationCost(defDuplicationCost);	// XXX MAGIC number!
									//
namespace segdup {

	double CSD(const EventCount& ec) {
		double cost(0.0);
		cost += ec.dups * duplicationCost;
		cost += ec.losses * lossCost;
		return cost;
	}

}

string hline("============================================================================\n");


void SelectNextConfiguration(CophyMultiMap& CMM, double T, vector<DupMove*>& moves, vector<double> probs) {
	// select which move type is to be used at random
	// for selected move type, Boltzmann sampling to get instance of move
	// perform the instance on CMM
	double total(0.0);
	for (auto p : probs)
		total += p;
	double d(dran(total));

	for (int i = 0; i < moves.size(); i++) {
		if (d < probs[i]) {
			//select this move
			moves[i]->apply(CMM, T);
			
			break;
		}
		d -= probs[i];
	}
}

void Algorithm2(CophyMultiMap& CMM, vector<DupMove*> moves, vector<double> probs, map<string, int>* sampledDistribution) {
	string bestMMap;
	EventCount bestEventCount;
	double bestCost(1e100);	// 10^100 should be enough!!
	string mapDescription;

	bool _debugging(false);
	DEBUG(cout << hline << "Algorithm1" << endl << hline);
	//DEBUG(cout << "Input multi-map:" << endl << CMM);
	string outputLine;
	char charBuffer[256];

	bool _doEarlyReconciliation(false);
	if (_doEarlyReconciliation) {
		CMM.doEarlyReconciliation();
	} else {
		CMM.doPageReconciliation();
	}
	CMM.toCompactString(bestMMap);
	DEBUG(cout << "Initial reconciliation complete:" << endl << CMM);

	CMM.putAllMoveableNodes();

	double T(0.0);
	ofstream ftrace;
	if (_saveTrace) {
		ftrace.open("segdup-trace.csv", std::ofstream::out);
		ftrace << "step,codivs,dups,losses,cost,temperature" << endl;
	}
	ostringstream bestPrettyMap;
	CMM.calcEventCount();
	bestEventCount = CMM.countEvents();
	bestCost = CSD(bestEventCount);

	DEBUG(cout << "number of movable nodes: " << CMM.getAllMoveableNodes().size() << endl);

	if (CMM.getAllMoveableNodes().size() != 0) {
		if (_verbose || !_silent) {
			cout << "i,#dups,#losses,cost\n";
		}
		for (int t(1); t <= nSteps + nFinal; ++t) {
			int nullMoves(0);
			EventCount ec;
			if (t <= nSteps)
				T = (Tinitial-Tfinal)*(1.0 - (1.0 * (t-1.0) / nSteps)) + Tfinal;
			else
				T = Tfinal;

//			if (t % 1000000 == 0)
//				cout << t << endl;
			//if (t == 2752238)

			SelectNextConfiguration(CMM, T, moves, probs);	// modifies CMM
													
			if (_saveSampledDistribution && (!_saveFinal || t > nSteps)) {
				if (t % outputInterval == 0) {
					CMM.toCompactString(mapDescription);
					EventCount totalEC = CMM.countEvents();
					mapDescription += "-D" + to_string(totalEC.dups) + "L" + to_string(totalEC.losses);
					(*sampledDistribution)[mapDescription] += 1;
				}
			}

			ec = CMM.countEvents();
			double cost = CSD(ec);
			if (cost < bestCost) {
				bestCost = cost;
				DEBUG(cout << "Best cost = " << bestCost << endl);
				bestEventCount = ec;
				DEBUG(if (bestCost < 0) {
					cout << "step " << t << ": best cost is NEGATIVE!" << endl;
					cout << "\tbest cost = " << bestCost << endl;
					cout << "\tbest event count = " << ec << endl;
					bestPrettyMap.str("");
					bestPrettyMap << CMM;
					exit(-1);
				});
				CMM.toCompactString(bestMMap);
				bestPrettyMap.str("");
				bestPrettyMap << CMM;
				DEBUG(
						cout << CMM << endl;
				);
			}
			if ((_verbose || !_silent) && (t % max(outputInterval,100) == 0)) {
				if (!_verbose) {
					for (unsigned int i(0); i < outputLine.size(); ++i) {
						cout << '\b';
					}
				}
//				outputLine.clear();
//				outputLine = to_string(t+1) + "," + to_string(ec.dups) + ",";
//				outputLine += to_string(ec.losses) + "," + to_string(CSD(ec));
				sprintf(charBuffer, "%d,%d,%d,%g", t, ec.dups, ec.losses, CSD(ec));
				outputLine = charBuffer;
				cout << outputLine;
				cout.flush();
				if (_verbose) {
					cout << endl;
				}
			}

			if (_saveTrace && (!_saveFinal || t > nSteps) ) {
				if (t % outputInterval == 0) {
	//					ftrace << sampleNumber << ",\"" << mapDescription << "\"," << to_string(CSD(ec)) << endl;
					ftrace << t << ',' << ec.codivs << ',' << ec.dups << ','
						<< ec.losses << ',' << to_string(CSD(ec)) << ','
						<< T
						<< endl;
				}
			}
	//		DEBUG(cout << nei.getLabel() << '\t' << nei.getScore() << endl);
		}
	}
	cout << endl;
	if (_saveSampledDistribution) {
		ofstream fout("segdup-samples.csv");
		fout << hline << "Sampled Distribution of Solutions:" << endl
				<< "Map-Events\tSamples" << endl;
		for (auto dis : (*sampledDistribution)) {
			fout << dis.first << "\t" << dis.second << endl;
		}
		fout.close();
	}
	cout << hline << "BEST Multiple CophyMap found by Algorithm 1:" << endl;
	cout << bestMMap << '\t' << endl;
	cout << bestEventCount << '\t' << bestCost << endl;
	cout << hline << endl;

	if (_saveTrace) {
		ftrace.close();
	}
	summaryfile << bestEventCount.codivs << ',' << bestEventCount.dups << ',' << bestEventCount.losses << ',' << bestCost << endl;
}

string segdupHelp("SegDup Help:\n"
		"\t>./segdup [options]\n"
		"\t? or -h\n\t\tto print this help message\n"
		"\t-S <newickformatspeciestree>\n"
		"\t\tNote that the species tree MUST be defined BEFORE the gene trees else the program will crash.\n"
		"\t\tAlso note that you MUST put trees in matched quotes if invoking from the command-line.\n"
		"\t-G <newickformatgenetree> <leafassociations>\n"
		"\t\tAssociation list MUST be a quoted string of space-separated pairs such as 'p:A q:B' to mean\n"
		"\t\tgene p is on species leaf A, and gene q is on species leaf B.\n"
		"\t-n <int>\n\t\tto supply the number of steps for Algorithm 1 (default: " + to_string(nSteps) + ")\n"
		"\t-Tinit <float>\n\t\tto supply the initial temperature (default: " + to_string(Tinitial) + ")\n"
		"\t-Tfinal <float>\n\t\tto supply the final temperature (default value " + to_string(Tfinal) + ")\n"
		"\t-nfinal <float>\n\t\tto supply the number of steps at the final temperature (default value " + to_string(nFinal) + ")\n"
		"\t--seed <int>\n\t\tto set the random number generator seed\n"
//		"\t--sat-spread <float>\n\t\tto set the Simulated Annealing \"spread\" parameter (default: " + to_string(SATempSpread) + ")\n"
//		"\t--sat-decay <float>\n\t\tto set the Simulated Annealing \"decay\" parameter (default: " + to_string(SATempDecay) + ")\n"
		"\t-d <float>\n\t\tto set the duplication event cost (default: " + to_string(duplicationCost) + ")\n"
		"\t-l <float>\n\t\tto set the loss event cost (default: " + to_string(lossCost) + ")\n"
		"\t-o (samples|trace|interval <int>|final)\n"
		"\t\tsamples to save the sampled distribution of maps (default value false)\n"
		"\t\ttrace to save a trace of the progress of the search (default value false)\n"
		"\t\tinterval to save only every few samples (default value 1)\n"
		"\t\tfinal to save only at the final temperature (default value false)\n"
		"\t--verbose\n\t\tto output lots of stuff (default: FALSE);\n"
		"\t--silent\n\t\tto output as little as possible (default: FALSE);\n"
	);

int main(int argn, char** argv) {
	bool _debugging(false);
	if (argn <= 1) {
		cout << segdupHelp << endl;
		return 0;
	}
	cout << "SegDup!" << endl;
	CophyMultiMap CMM;
	summaryfile.open("summary.csv", std::ios_base::app);
	summaryfile << "codivs,dups,losses,cost\n";
	vector<CophyMap*> M;
	uint numGeneTrees(0);
	Tree *S(nullptr);
	vector<Tree*> G;
	for (int i(1); i < argn; ++i) {
		DEBUG(cout << "Parsing argument " << i << " = " << argv[i] << endl);
		if (!strcmp(argv[i], "?") || !strcmp(argv[i], "-h")) {
			cout << segdupHelp;
			return 0;
		}
		if (!strcmp(argv[i], "-S")) {
			++i;
			string newick(argv[i]);
			S = new Tree('s', newick);
			S->setLabel("S");
			if (!_silent) {
				cout << "Input Species tree S:" << endl << (*S) << endl;
			}
		} else if (!strcmp(argv[i], "-G")) {
			++i;
			string newick(argv[i]);
			if (newick[0] != '(') {
				cerr << "WARNING: Expecting a non-trivial Newick-format tree here, but got this:\n\t"
						<< newick
						<< "\nSegDup is skipping this argument and the next, which should be a set of associations.\n";
				++i;	// skip this AND the next argument
				continue;	// ... and go on to the next argument.  This isn't a non-trivial tree.
			}
			++numGeneTrees;
			G.push_back(new Tree('g', newick));
			Tree *P = G[numGeneTrees-1];
			P->setLabel("G" + to_string(numGeneTrees));
			if (!_silent) {
				cout << "Input Gene tree " << P->getLabel() << ":" << endl << (*P) << endl;
			}
			++i;
			string assoc(argv[i]);
			NodeMap* A = new NodeMap(S, P, assoc);
			if (!_silent) {
				cout << "Input Associations:" << endl << (*A);
			}
			CophyMap* M = new CophyMap(*A);
			P->setInfo(M->getInfo());
			P->setShowInfo(true);
			CMM.addCophyMap(M);
		} else if (!strcmp(argv[i], "-n")) {
			++i;
			nSteps = atoi(argv[i]);
			if (!_silent) {
				cout << "Setting Number of steps to " << nSteps << endl;
			}
		} else if (!strcmp(argv[i], "-d")) {
			++i;
			duplicationCost = atof(argv[i]);
			if (!_silent) {
				cout << "Setting DuplicationCost to " << duplicationCost << endl;
			}
			CMM.setDuplicationCost(duplicationCost);
		} else if (!strcmp(argv[i], "-l")) {
			++i;
			lossCost = atof(argv[i]);
			CMM.setLossCost(lossCost);
			cout << "Setting LossCost to " << lossCost << endl;
		} else if (!strcmp(argv[i], "-Tinit")) {
			++i;
			Tinitial = atof(argv[i]);
			cout << "Setting Initial Temperature to " << Tinitial << endl;
		} else if (!strcmp(argv[i], "-Tfinal")) {
			++i;
			Tfinal = atof(argv[i]);
			cout << "Setting Final Temperature to " << Tfinal << endl;
		} else if (!strcmp(argv[i], "-nfinal")) {
			++i;
			nFinal = atoi(argv[i]);
			cout << "Setting Number of steps at final temperature to " << nFinal << endl;
		} else if (!strcmp(argv[i], "-o")) {
			++i;
//			if (!strcmp(argv[i], "probs")) {
//				_outputProbabilities = true;
//			} else
			if (!strcmp(argv[i], "samples")) {
				if (!_silent) {
					cout << "Setting Show_Samples to true" << endl;
				}
				_saveSampledDistribution = true;
			} else if (!strcmp(argv[i], "interval")) {
				cout << "Output every n steps" << endl;
				++i;
				outputInterval = atoi(argv[i]);
			} else if (!strcmp(argv[i], "trace")) {
				cout << "setting save a trace to true" << endl;
				_saveTrace = true;
			} else if (!strcmp(argv[i], "final")) {
				cout << "Setting save at final temperature to true" << endl;
				_saveFinal = true;
			}
//		} else if (!strcmp(argv[i], "--sat-spread")) {
//			++i;
//			SATempSpread = atof(argv[i]);
//			if (_verbose) {
//				cout << "Setting Simulated Annealing \"spread\" parameter to " << SATempSpread << endl;
//			}
//		} else if (!strcmp(argv[i], "--sat-decay")) {
//			++i;
//			SATempDecay = atof(argv[i]);
//			if (!_silent) {
//				cout << "Setting Simulated Annealing \"decay\" parameter to " << SATempDecay << endl;
//			}
		} else if (!strcmp(argv[i], "--verbose")) {
			_verbose = true;
		} else if (!strcmp(argv[i], "--silent")) {
			_silent = true;
		} else if (!strcmp(argv[i], "--seed")) {
			++i;
			seed = atoi(argv[i]);
			cout << "Setting random number seed to " << seed << " for testing / repeatability." << endl;
			generator.seed(seed);
		} else {
			cout << "I cannot understand this argument: \"" << argv[i] << "\", which is number " << i << " in the input." << endl
					<< "Rather than continue with a potentially erroneous input I am quitting." << endl;
			return 1;
		}
	}
	cout << hline;
	//CMM.doPageReconciliation();
	map<string, int> sampledDistribution;
//	Algorithm1(CMM, sampledDistribution);
	std::vector<DupMove*> moves;
	moves.push_back(new SingleNodeMove);
	//moves.push_back(new EmptyMove);
	moves.push_back(new SingleVertexMove);
	std::vector<double> probs;
	probs.push_back(0.5);
	//probs.push_back(0.5);
	probs.push_back(0.5);
	Algorithm2(CMM, moves, probs, &sampledDistribution);


//	bool _debugging(true);

//	testPageReconciliation2();
//	testPageReconciliation3();
//	testPerfectMatchReconciliation();
//	testAvailableHosts();
//	testDuplicationHeight();
//	ybcTestCases();

//	doAlgorithmTest2();
//	doContenderTest();
	//	cout << exp(-31.0) << endl;

//	set<Contender> C;
//	Contender a(0.1,"fred");
//	Contender b(0.2,"harry");
//	Contender c(0.1,"mary");
//	Contender d(0.3,"harry");
//	C.insert(a);
//	C.insert(b);
//	C.insert(c);
//	C.insert(d);
//	for (auto x : C) {
//		cout << x << endl;
//	}
	return 0;
}

