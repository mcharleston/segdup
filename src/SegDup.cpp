/*
 * SegDup.cpp
 *
 *  Created on: 2 Aug 2022
 *      Author: mac
 */

#include <chrono>
#include <cmath>	// for exp
#include <cstring> // for strcmp
#include <functional> // for transform
#include <map>
#include <random>
#include <sstream>
#include <fstream>
#include <stdio.h>

#include "Node.h"
#include "Tree.h"
#include "CophyMap.h"
#include "CophyMultiMap.h"
#include "Contender.h"

#include "../utility/appexception.h"
#include "../utility/debugging.h"

/**
 * Input to an instance includes whether each node is mapped to a D or S
 */

using namespace std;
using namespace segdup;

bool _debugging(true);
bool _silent(false);
bool _outputProbabilities(false);
bool _saveTrace(false);
bool _showSampledDistribution(false);
bool _verbose(false);
int nSteps(1000);
double Tinitial(0.1);
double SATempSpread(100);
double SATempDecay(4);

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

std::default_random_engine generator(seed);
std::uniform_real_distribution<float> runif(0.0, 1.0);
std::uniform_real_distribution<double> dunif(0.0, 1.0);

float fran() {
	return runif(generator);
}

double dran(double mult = 1.0) {
	return mult * dunif(generator);
}

uint plran(float l, float u, float r) {
	float y(fran());
	float ex(r+1.0);
	return static_cast<uint>( std::pow( std::pow(u,ex) - std::pow(l, ex)*y + std::pow(l, ex), 1.0/ex ) );
}

double lossCost(defLossCost);						// XXX
double duplicationCost(defDuplicationCost);	// XXX MAGIC number!

double CSD(const EventCount& ec) {
	double cost(0.0);
	cost += ec.dups * duplicationCost;
	cost += ec.losses * lossCost;
	return cost;
}

string hline("============================================================================\n");


void roughTest() {
	Tree T;
	Node* r = new Node("r");
	T.setRoot(r);
	Node *x = new Node("x");
	Node *a = new Node("a"), *b = new Node("b"), *c = new Node("c");
	r->bifurcate(x, c);
	x->bifurcate(a, b);
	cout << T;
	T.gatherVertices();
	map<string,Node*>::iterator it = T.getVertices().begin();
	for ( ; it != T.getVertices().end(); ++it) {
		cout << it->first << " -> " << *(it->second) << endl;
	}
	Tree G1("blue"), G2("red");
	Node* br = new Node("br");
	G1.setRoot(br);
	Node *bx = new Node("bx"), *by = new Node("by"), *ba1 = new Node("ba1"), *ba2 = new Node ("ba2"), *bb = new Node("bb"), *bc = new Node("bc");
	br->bifurcate(bx, bc);
	bx->bifurcate(by, bb);
	by->bifurcate(ba1, ba2);
	Node *w = new Node("w");
	Node *v = new Node("v");
	bb->setSibling(v);
	ba2->setSibling(w);
	cout << G1 << endl;
	G1.compressTraverseWrite(cout);

	Tree H("(a,(b,c));");
	H.getRoot()->setLabel("h0");
	Tree P("(p,(q,r))");
	std::map<std::string, Node*>& VH = H.getVertices();
//	std::map<std::string, Node*>& VP = P.getVertices();
	cout << "Finding a vertex by name: looking for a in V(H) (which should exist):" << endl;
	cout << VH["a"]->getLabel() << endl;
	cout << "Did it work?" << endl;
	Node* ha = VH["a"];
	if (ha == nullptr) {
		cout << "ha is null" << endl;
	}
	Node* hb = VH["b"];
	if (hb == nullptr) {
		cout << "hb is null" << endl;
	}
	Node* hc = VH["c"];
	if (hc == nullptr) {
		cout << "hc is null" << endl;
	}
	Node* hr = ha->getParent();
	if (hr == nullptr) {
		cout << "hr is null" << endl;
	}
	Node* v1 = VH["v1"];
	if (v1 == nullptr) {
		cout << "v1 is null" << endl;
	}
	cout << "Root of H:" << hr->getLabel() << endl;
	Node* lcabc = H.LCA(hb,hc);
	cout << "lcabc = " << lcabc->getLabel() << "; should be " << v1->getLabel() << endl;
	H.compressTraverseWrite(cout);
	P.compressTraverseWrite(cout);
	CophyMap CM(&H,&P);
	CM.setPhi("p", "a");
	CM.setPhi("q", "b");
	CM.setPhi("r", "c");
	cout << H.LCA(ha,hb)->getLabel() << endl;
	cout << H.LCA(ha,hc)->getLabel() << endl;
	cout << H.LCA(hc,hb)->getLabel() << endl;
	cout << H.LCA(ha,ha)->getLabel() << endl;
	CM.doPageReconciliation();

}

void test1() {
	//	Tree H("((a,b),c));");
	//	H.getRoot()->setLabel("h0");
	//	Tree P("(pa,(pb,pc));");
	//	P.getRoot()->setLabel("p0");
	////	cout << "Vertices of P:" << endl;
	////	for (auto pr : P.getVertices()) {
	////		cout << pr.first << " at " << pr.second << endl;
	////	}
	////	cout << H;
	////	cout << P;
	//	NodeMap Phi( {
	//				pair<Node*, Node*>(P["pa"], H["a"]),
	//				pair<Node*, Node*>(P["pb"], H["b"]),
	//				pair<Node*, Node*>(P["pc"], H["c"]),
	//				} );
	////	Phi.setPHAssociation(P["pa"], H["b"], 0);
	////	cout << "P[pa] = " << P["pa"] << endl;
	////	cout << "P.at(pa) = " << P.at("pa") << endl;
	////	cout << Phi;
	////	cout << "Attempting to set associations using string labels" << endl;
	////	cout << "&P = " << &P << endl;
	////	cout << "&H = " << &H << endl;
	////	Phi.setPHAssociation("pa", "a", 0);
	////	cout << "Finished setting associations using string labels" << endl;
	//	cout << Phi;
	//	CophyMap M(&H, &P, Phi);
	//	DEBUG(cout << "CophyMap M:" << endl << "Host =\n" << H << "Parasite = \n" << P);
	//	M.doPageReconciliation();
	//	cout << M.getPhi();
	//	P.setShowInfo(true);
	//	cout << P;
	//	EventCount E = M.countEvents();
	//	cout << "Event counts: " << E << endl;

}

void testPageReconciliation2() {
	Tree H2("((ha,hb),(hc,hd))");
	H2.getRoot()->setLabel("h0");
	Tree P2("((pa,(pb,pc)),pd)");
	P2.getRoot()->setLabel("p0");
	cout << "P2 = " << endl << P2;
	cout << "H2 = " << endl << H2;
	NodeMap Assocs(&H2, &P2, "pa:ha, pb:hb, pc:hc, pd:hd");
	CophyMap M2(Assocs);
	M2.doPageReconciliation();
	P2.setShowInfo(true);
	cout << "P2 = " << endl << P2;
	cout << "H2 = " << endl << H2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts: " << E2 << endl;

}

void testPageReconciliation3() {
	Tree H("((ha,hb),(hc,hd))");
	H.getRoot()->setLabel("h0");
	Tree P("((pa,(pb,pc,pe)),pd)");
	P.getRoot()->setLabel("p0");
	cout << "P = " << endl << P;
	cout << "H = " << endl << H;
	NodeMap Assocs(&H, &P, "pa:ha, pb:hb, pc:hc, pe:hc, pd:hd");
	CophyMap M(Assocs);
	M.doPageReconciliation();
	P.setShowInfo(true);
	cout << "P = " << endl << P;
	cout << "H = " << endl << H;
	EventCount E = M.countEvents();
	cout << "Event counts: " << E << endl;
}

//void testPerfectMatchReconciliation() {
//	Tree H("(((a,b),c),(d,e))");
//	H.getRoot()->setLabel("h0");
//	Tree P("(((a,b),c),(d,e))");
//	P.getRoot()->setLabel("p0");
//	cout << "P = " << endl << P;
//	cout << "H = " << endl << H;
//	NodeMap Assocs(&H, &P, "a:a, b:b, c:c, d:d, e:e");
//	CophyMap M(Assocs);
//	M.doPageReconciliation();
//	P.setShowInfo(true);
//	cout << "P = " << endl << P;
//	cout << "H = " << endl << H;
//	EventCount E = M.countEvents();
//	cout << "Event counts: " << E << endl;
//}
//
void testAvailableHosts() {
	Tree H("(((A,B),C),D)");
	H.setLabel("H");
	Tree G("(((a,b),d),c)");
	G.setLabel("G");
	G.setNodeLabel(G.LCA(G["a"], G["b"]), "p");
	G.setNodeLabel(G.LCA(G["d"], G["p"]), "q");
	G.setNodeLabel(G.LCA(G["c"], G["q"]), "r");
	H.setNodeLabel(H.LCA(H["A"], H["B"]), "W");
	H.setNodeLabel(H.LCA(H["W"], H["C"]), "X");
	H.setNodeLabel(H.LCA(H["X"], H["D"]), "Y");
	NodeMap phi(&H,&G,"a:A b:B c:C d:D");
	CophyMap M(phi);
	G.setShowInfo(true);
	M.doPageReconciliation();
	M.storeAssociationInfo();
	cout << M;
	M.calcAvailableNewHosts(G["a"]);
	cout << G;
	cout << M;
	cout << M.countEvents() << endl;
	std::set<pair<Node*, eventType>> avail = M.calcAvailableNewHosts(G["p"]);
	for (auto a : avail) {
		cout << "available host : " << a.first->getLabel() << ", event type " << a.second << endl;
	}
	M.moveToHost(G["p"], H["X"]);
	M.inferEvents();
	M.storeAssociationInfo();
	cout << M;
	cout << M.countEvents() << endl;
}

void testDuplicationHeight() {
	Tree S("(((A,B),(C,D)),E)");
	Tree G("((((a1,b1),((a2,b2),(a3,b3))),(((c1,d1),c2),d2)),e)");
//	cout << S;
//	cout << G;
	NodeMap A1(&S, &G, "a1:A a2:A a3:A b1:B b2:B b3:B c1:C c2:C d1:D d2:D e:E");
	CophyMap M(A1);
	M.doPageReconciliation();
	M.storeAssociationInfo();
	G.setShowInfo(true);
	cout << M << M.countEvents() << endl;

	Tree G2("(((((a1,b1),(a4,b4)),((a2,b2),(a3,b3))),((c1,d1),c2)),e)");
	NodeMap A2(&S, &G2, "a1:A a2:A a3:A a4:A b1:B b2:B b3:B b4:B c1:C c2:C d1:D e:E");
	CophyMap M2(A2);
	M2.doPageReconciliation();
	M2.storeAssociationInfo();
	G2.setShowInfo(true);
	cout << M2 << M2.countEvents() << endl;

	Tree G3("(((c1,d1),((c2,d2),((c3,d3),(c4,d4)))),e)");
	NodeMap A3(&S, &G3, "c1:C c2:C c3:C c4:C d1:D d2:D d3:D d4:D e:E");
	CophyMap M3(A3);
	M3.doPageReconciliation();
	M3.storeAssociationInfo();
	G3.setShowInfo(true);
	cout << M3 << M3.countEvents() << endl;

	CophyMultiMap CMM;
	CMM.addCophyMap(&M);
	CMM.addCophyMap(&M2);
	CMM.addCophyMap(&M3);
	int h = CMM.calcDuplicationHeight(S["v4"]);
	cout << "dh(v4) = " << h << endl;
	h = CMM.calcDuplicationHeight(S["v3"]);
	cout << "dh(v3) = " << h << endl;
	h = CMM.calcDuplicationHeight(S["v5"]);
	cout << "dh(v5) = " << h << endl;
	cout << "Counting events with segmental duplications" << endl << CMM.countEvents() << endl;
}

void Algorithm1(CophyMultiMap& CMM, map<string, int>& sampledDistribution) {
	bool _debugging(false);
	DEBUG(cout << hline << "Algorithm1" << endl << hline);
	DEBUG(cout << "Input multi-map:" << endl << CMM);
	CMM.doPageReconciliation();

	DEBUG(cout << "Initial reconciliation complete:" << endl << CMM);
//	DEBUG(for (auto mpr : CMM.getMaps()) {
//		CophyMap* M = mpr.second;
//		cout << M->getParasiteTree()->getLabel() << " event counts: " << M->countEvents() << endl;
//	});

//	DEBUG(cout << "total event counts: " << CMM.countEvents() << endl);
	CophyMultiMap nuMM(CMM);
	map<CophyMap*, set<Node*>> iV;	// internal vertices of each parasite / gene tree
	for (auto mpr : CMM.getMaps()) {
		CophyMap* M = mpr.second;string segdupHelp("SegDup Help:\n"
				"\t>./segdup -[hvDn]");
		Tree* G = M->getParasiteTree();
		G->putInternalVertices(iV[M]);
//		DEBUG(cout << "Tree " << G->getLabel() << " has internal vertices { ";
//			for (Node* n : iV[M]) {
//				cout << n->getLabel() << " ";
//			}
//			cout << "}" << endl;
//		);
	}

	double T(0.0);
//	double BoltzmannConst(1.3806503e-23);
	double fudgeFactor(100.0);
	std::set<Contender> neighbours;
	EventCount ecoriginal = CMM.countEvents();
	cout << "Initial" << endl;
	cout << ecoriginal.codivs << ',' << ecoriginal.dups << ',' << ecoriginal.losses << endl;

	string mapDescription;
	CMM.toCompactString(mapDescription);
	if (_verbose) {
		cout << "ORIGINAL MAP:" << endl << CMM << "Events\tScore\tMap\n" << ecoriginal
			<< '\t' << CSD(ecoriginal) << '\t' << mapDescription << endl;
	}
	string bestMMap;
	EventCount bestEventCount;
	double bestCost(1e100);	// 10^100 should be enough!!
	ofstream ftrace;
	int sampleNumber(0);
	if (_saveTrace) {
		ftrace.open("segdup-trace.csv", std::ofstream::out);
		ftrace << "i,T,dups,losses,csd\n";
	}
	ostringstream bestPrettyMap;

	EventCount ec, oldEC;
	int noChangeInHeights(0), heightsChanged(0);
	for (int t(0); t < nSteps; ++t) {
		DEBUG(cout << hline << "t = " << t << endl << hline << endl);
		int nullMoves(0);
		T = Tinitial*(1.0 - (1.0 * t / nSteps)); // function 0 / linear
//		T = Tinitial*(1.0 - (1.0 * t / nSteps))*(1.0 - (1.0 * t / nSteps)); // function 1 / linear2
//		double x = 1.0 * t / nSteps;
//		double tempParam(2.5);
//		T = Tinitial*(exp(-tempParam * x * x * x));	// function 2 / exp3	// XXX not very good on ybc-case3
//		T = Tinitial*(1.0/(1-sqrt(2)))*(1-sqrt(1.0 + 1.0/(1.0 +tempParam * x * x * x))); // function 3 / sqrts // XXX not very good on ybc-case3
//		T = Tinitial*1.5/(1.0 + 0.2 * sqrt(x*x*t+1.0));	// function 4 / sqrts	// XXX looks good on ybc-case3 and on Guido et al..
//		double paramA(4.4);
//		double paramB(400.0);
//		T = Tinitial*(1.0 / (1.0 + log(1 + pow(t/SATempSpread, SATempDecay)))); // function 5 / logpow
		double total(0.0);
		neighbours.clear();
		oldEC = CMM.countEvents();
		ec = oldEC;
		for (auto mpr : CMM.getMaps()) {
			CophyMap* M = mpr.second;
//			DEBUG(cout << hline << "ORIGINAL MAP:" << (*M));
			for (Node* p : iV[M]) {
				set<pair<Node*, eventType>> nextImages = M->calcAvailableNewHosts(p);
				DEBUG(
					cout << "Node " << p->getLabel() << " has possible images { ";
					for (auto a : nextImages) {
							cout << eventSymbol[a.second] << a.first->getLabel() << " ";
						}
					cout << "}" << endl
				);
				for (auto a : nextImages) {
//					DEBUG(cout << "***" << endl);
					if (a.first == M->getHost(p) && a.second == M->getEvent(p)) {
//						DEBUG(cout << "No change: " << p->getLabel() << " staying on " << a.first->getLabel() << endl);
//						DEBUG(cout << "Combined event counts: " << ecoriginal << "; cost = " << CSD(ecoriginal) << endl);
//						DEBUG(cout << "Probability of sampling proportional to: " << exp(-CSD(ecoriginal) / T) << endl);

//						ec = oldEC;//CMM.countEvents();	// XXX should just use the old one but if I do, it all breaks!
						double score = exp(-CSD(ec) / ( fudgeFactor * T));
						Contender noChange( score, p, a.first, a.second, M );

						noChange._noMove = true;
						noChange.setLabel("leaving " + p->getLabel() + " on [" + eventSymbol[a.second] + "]" + a.first->getLabel());

						neighbours.insert(noChange);
//						DEBUG(cout << "Adding a noChange move: " << noChange << endl);
						++nullMoves;
						continue;
					}
//					DEBUG(cout << "Testing moving " << p->getLabel() << " to host " << eventSymbol[a.second] << a.first->getLabel() << endl);
					Node* oldHost = M->getHost(p);
					eventType oldEvent = M->getEvent(p);
					int oldSourceDuplicationHeight = CMM.calcDuplicationHeight(oldHost);
					int oldTargetDuplicationHeight = CMM.calcDuplicationHeight(a.first);
					M->moveToHost(p, a.first, a.second);
					CMM.invMap[oldHost].erase(p);
					CMM.invMap[a.first].insert(p);
					// do something here about checking the heights...
					if ((oldSourceDuplicationHeight != CMM.calcDuplicationHeight(oldHost)) ||
							(oldTargetDuplicationHeight != CMM.calcDuplicationHeight(a.first))) {
//						cout << CMM;
//						cout << "p moving " << p->getLabel() << endl;
//						cout << "old host: " << oldHost->getLabel() << endl;
//						cout << "new host: " << a.first->getLabel() << endl;
//						cout << "old event: " << M->describeEvent(oldEvent) << endl;
//						cout << "new event: " << M->describeEvent(a.second) << endl;
//						cout << "oldSourceDuplicationHeight = " << oldSourceDuplicationHeight << endl;
//						cout << "oldTargetDuplicationHeight = " << oldTargetDuplicationHeight << endl;
//						cout << "new height for previous host = " << CMM.calcDuplicationHeight(oldHost) << endl;
//						cout << "new height for new host = " << CMM.calcDuplicationHeight(a.first) << endl;
						++heightsChanged;
//						return;
						ec = CMM.countEvents();
					} else {
						++noChangeInHeights;
						ec = oldEC;
						ec = CMM.countEvents();
//						cout << "old EC: " << oldEC << ";\t";
//						cout << "counted EC: " << ec << endl;
						CMM.toCompactString(mapDescription);
						CMM.storeEventCount(mapDescription, ec);
					}
					// XXX I just CANNOT work out why this isn't working. I'm going to have to leave it until I've had more sleep, and just run simulations the slow way.
					// XXX Fix after grant application goes in, perhaps.
//					ec = CMM.countEvents();
					/**
					 * SPEEDUP ideas:
					 * keep track of all tree duplication heights at this branch;
					 * or that and the max. If the current one is lowering in height past this max then we have to recalculate the height
					 * of this branch; otherwise we don't. Etc.
					 */
//					DEBUG(cout << "New event counts: " << ec << endl);
					Association ass(p, a.first, a.second);
					double score = exp(-1.0* CSD(ec) / (fudgeFactor * T));
					DEBUG(cout << "CSD(ec) = " << CSD(ec) << endl);
					DEBUG(cout << "Exponent = " << (-1.0 *CSD(ec) / (fudgeFactor * T)) << endl);
					DEBUG(cout << "Score = " << score << endl);
					Contender con( score, p, a.first, a.second, M );
					string label = "moving " + p->getLabel() + " to [" + eventSymbol[a.second] + "]" + a.first->getLabel();
					con.setLabel(label);
					neighbours.insert(con);
//					DEBUG(cout << "Complete Multi-Map:" << CMM << "Total multi-map event counts: " << ec << endl);
//					DEBUG(cout << "Probability of sampling proportional to: " << con.getScore() << endl);
					M->moveToHost(p, oldHost, oldEvent);
					CMM.invMap[oldHost].insert(p);
					CMM.invMap[a.first].erase(p);

				}
				total = 0.0;
				for (auto nei : neighbours) {
					total += nei.getScore();
//					DEBUG(cout << "total = " << total << endl);
				}
//				DEBUG(cout << "total score = " << total << endl);
//				DEBUG(cout << "Summary of sampling options, events & costs:" << endl);
//				for (auto nei : neighbours) {
//					DEBUG(cout << "\t" << nei.getLabel() //<< " " <<  nei.getParasite()->getLabel() << ':' << nei.getHost()->getLabel()
//							<< " has cost " << nei.getScore() << endl);
//				}
//				DEBUG(cout << '\t' << nullMoves << " null moves with cost " << CSD(ecoriginal) << endl);
//				 DO THE SAMPLING HERE
//				DEBUG(
//						cout << "Scores of neighbours: ";
//						for (auto nei : neighbours) {
//							cout << nei.getScore() << " ";
//						}
//						cout << endl
//						);
			}
//				DEBUG(cout << endl);
		}
		double r = dran(total);
		DEBUG(cout << "Total probability proportional to " << total << endl);
		if (neighbours.size() == 0) {
			throw new app_exception("No neighbours at all!!");
		}
		for (auto nei : neighbours) {
			DEBUG(cout << "This neighbour has score " << nei.getScore() << endl);
			if (r <= nei.getScore()) {
				if (nei._noMove) {
					DEBUG(cout << "Selected move: No change (probability = " << (nei.getScore()/total) << ")" << endl);
				} else {
					// sample this one
					DEBUG(cout << "Selected move: " << nei.getLabel() << " (probability = " << (nei.getScore()/total) << ")" << endl);
					nei.getMap()->moveToHost(nei.getParasite(), nei.getHost(), nei.getEvent());
					nei.getMap()->checkValidHostOrdering();
//							DEBUG(cout << nei.getLabel() << endl << *(nei.getMap()->getParasiteTree()));
//							DEBUG(cout << t << '\t' << nei.getLabel() << endl);
				}
				CMM.toCompactString(mapDescription);
				sampledDistribution[mapDescription] += 1;
				ec = CMM.getEventCount(mapDescription);
				double cost = CSD(ec);
				if (cost < bestCost) {
					bestCost = cost;
					bestEventCount = ec;
					bestMMap = mapDescription;
					bestPrettyMap.str("");
					bestPrettyMap << CMM;
				}
				DEBUG(cout << nei.getLabel() << '\t' << nei.getScore() << endl);
				if (_saveTrace) {
					++sampleNumber;
					ftrace << sampleNumber << ',';
//					ftrace << mapDescription << "\",";
					ftrace << T << ',';
					ftrace << ec.dups << ',';
					ftrace << ec.losses << ',';
					ftrace << to_string(CSD(ec)) << endl;
					DEBUG(cout << "\tSaving new sampled solution" << endl);
				}
				break;
			}
			DEBUG(cout << "r reducing from " << r);
			r -= nei.getScore();
			DEBUG(cout << " to " << r << endl);
		}
		oldEC = ec;
	}
	bool _showDistribution(false);
	if (_showDistribution) {
		cout << hline << "Sampled Distribution of Solutions:" << endl
				<< "Event Counts; Score";
		cout << "\t";
//		if (_outputProbabilities) {
//			cout << "\tProb";
//		}Input Species tree:
		cout << "\tnumSamples/" << nSteps << endl;
		for (auto dis : sampledDistribution) {
//			cout << "Looking for event count for this map: " << dis.first << endl;
			EventCount ec = CMM.getEventCount(dis.first);
			cout << ec << '\t' << dis.first << "\t" << dis.second << endl;
		}
		cout << hline << "FINAL Multiple CophyMap found by Algorithm 1:" << endl;
		EventCount ecFinal(CMM.countEvents());
		cout << CMM << ecFinal << endl << hline << endl;
		cout << "final CSD: " << CSD(ecFinal) << endl;
	}
//	cout << "FINAL Multiple CophyMap found by Algorithm 1:" << endl;
//	EventCount ecFinal(CMM.countEvents());
//	cout << CMM << ecFinal << endl << hline << endl;
//	cout << "final CSD: " << CSD(ecFinal) << endl;
//	if (_verbose) {
//		cout << hline << "BEST Multiple CophyMap found by Algorithm 1:" << endl;
		cout << "nCospec,nSegDup,nLoss" << endl;
		cout << bestEventCount.codivs << ',' << bestEventCount.dups << ',' << bestEventCount.losses << endl;
//		cout << bestEventCount << '\t' << bestCost << '\t' << bestMMap << '\t' << endl;
//		cout << bestPrettyMap.str();
//		cout << hline << endl;
//	}
	if (_saveTrace) {
		ftrace.close();
	}
//	cout << "Number of times heights don't change: " << noChangeInHeights << endl;
//	cout << "Number of times heights DO change: " << heightsChanged << endl;
}

void doTestCase1() {
	cout << "Test case 1 (trivial match)" << endl;
	Tree T1("(A,(B,C))");
	Tree S1("(a,(b,c))");
	T1.setLabel("T1");
	S1.setLabel("S1");
	S1.setShowInfo(true);
	NodeMap assocs1(&T1, &S1, "a:A b:B c:C");
	CophyMap M1(assocs1);
	M1.doPageReconciliation();
	M1.inferEvents();
	M1.storeAssociationInfo();
	cout << M1 << M1.countEvents() << endl << hline;
}
void doTestCase2() {
	cout << "Test case 2:" << endl;
	Tree S("(A,(B,C))");
	Tree G("(a1,(a2,((b,c1),c2)))");
	S.setLabel("S");
	G.setLabel("G1");
	cout << S.getLabel() << ":\n" << S << endl;
	cout << G.getLabel() << ":\n" << G << endl;
	NodeMap Assocs(&S, &G, "a1:A, a2:A, b:B, c1:C, c2:C");
	CophyMap M(Assocs);
	M.doPageReconciliation();
	M.inferEvents();
	M.storeAssociationInfo();
	G.setShowInfo(true);
	cout << "G:" << endl << G;
	EventCount E = M.countEvents();
	cout << "Event counts: " << E << endl << hline;
}

void doTestCase3a() {
	cout << "Test case 3a: LCA, FIRST gene tree" << endl;
	cout << "Testing moving a node and correctly calculating the events again";
	Tree S("(A,(B,C))");
	S.setLabel("S");
	Tree G1("(a,(b1,((b2,c1),c2)))");
	Tree G2("((a1,a2),((b,c1),c2)))");
	G1.setLabel("BLUE");
	G1.setShowInfo(true);
	G2.setLabel("GREEN");
	G2.setShowInfo(true);
	NodeMap Assocs1(&S, &G1, "a:A b1:B b2:B c1:C c2:C");
	CophyMap M1(Assocs1);
	M1.doPageReconciliation();
	M1.inferEvents();
	M1.storeAssociationInfo();
	cout << S;
	cout << G1 << M1.countEvents() << endl << hline;
	cout << "Now moving LCA(b2,c1) to a DUPLICATION on the host node above." << endl;
	Node* p = G1.LCA("b2", "c1");
	Node* h = S.LCA("B", "C");
	M1.moveToHost(p, h, duplication);
	M1.storeAssociationInfo();
	cout << G1 << M1.countEvents() << endl << hline;
}
void doTestCase3() {
	CophyMultiMap CMM;
	cout << "Test case 3: LCA, 2 gene trees" << endl;
	Tree S("(A,(B,C))");
	S.setLabel("S");
	Tree G1("(a,(b1,((b2,c1),c2)))");
	Tree G2("((a1,a2),((b,c1),c2)))");
	G1.setLabel("BLUE");
	G1.setShowInfo(true);
	G2.setLabel("GREEN");
	G2.setShowInfo(true);
	NodeMap Assocs1(&S, &G1, "a:A b1:B b2:B c1:C c2:C");
	CophyMap M1(Assocs1);
	M1.doPageReconciliation();
	M1.inferEvents();
	Node* p = G1.LCA("b2", "c1");
	Node* h = S.LCA("B", "C");
	M1.moveToHost(p, h, duplication);
	M1.storeAssociationInfo();
	cout << "M1 events (BLUE tree): " << M1.countEvents() << endl;
	NodeMap Assocs2(&S, &G2, "a1:A, a2:A, b:B, c1:C, c2:C");
	CophyMap M2(Assocs2);
	M2.doPageReconciliation();
	M2.inferEvents();
	M2.storeAssociationInfo();
	cout << "M2 events (GREEN tree): " << M2.countEvents() << endl;
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G1.getLabel() << endl << G1;
	cout << "Gene tree " << G2.getLabel() << endl << G2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts for " << G2.getLabel() << ":" << endl << E2 << endl;
	CMM.addCophyMap(&M1);
	CMM.addCophyMap(&M2);
	cout << " Test case 3 event counts: " << CMM.countEvents() << endl << hline;
}
void doTestCase4() {
	cout << "Test case 4: Not LCA, 1 gene tree" << endl;
	Tree S("(A,(B,C))");
	Tree G("(a,((b1,b2),((b3,c1),c2)))");
	G.setLabel("G");
	G.setShowInfo(true);
	NodeMap assocs(&S, &G, "a:A b1:B b2:B b3:B c1:C c2:C");
	CophyMap M(assocs);
	M.doPageReconciliation();
	Node* p = G.LCA("b1", "b2");
	Node* newHost = S.LCA("B", "C");
	M.moveToHost(p, newHost);
	M.inferEvents();
	// move first child
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G.getLabel() << endl << G;
	CophyMultiMap MM;
	MM.addCophyMap(&M);
	cout << "Test case 4 event counts:\n" << MM.countEvents() << endl << hline;
}
void doTestCase5() {
	cout << "Test case 5: Not LCA, 2 gene trees" << endl;
	Tree S("(A,(B,C))");
	S.setLabel("S");
	Tree G1("(a,(b1,b2))");
	Tree G2("(a,(c1,c2))");
	G1.setLabel("G1");
	G1.setShowInfo(true);
	G2.setLabel("G2");
	G2.setShowInfo(true);
	NodeMap Assocs1(&S, &G1, "a:A b1:B b2:B");
	CophyMap M1(Assocs1);
	M1.doPageReconciliation();
	Node* h = S.LCA("B", "C");
	Node* p = G1.LCA("b1", "b2");
	M1.moveToHost(p,h);
	M1.inferEvents();
	NodeMap Assocs2(&S, &G2, "a:A c1:C c2:C");
	CophyMap M2(Assocs2);
	M2.doPageReconciliation();
	p = G2.LCA(G2["c1"], G2["c2"]);
	M2.moveToHost(p,h);
	M2.inferEvents();
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G1.getLabel() << endl << G1;
	cout << "Gene tree " << G2.getLabel() << endl << G2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts for " << G2.getLabel() << ":" << endl << E2 << endl;
	CophyMultiMap MCM;
	MCM.addCophyMap(&M1);
	MCM.addCophyMap(&M2);
	cout << " Test case 5 event counts: " << MCM.countEvents() << endl << hline;
}
void doTestCase6() {
	cout << "Test case 6: Not LCA, 2 gene trees" << endl;
	Tree S("(A,(B,C))");
	S.setLabel("S");
	Tree G1("(a1,a2)");
	Tree G2("(c1,c2)");
	G1.setLabel("G1");
	G1.setShowInfo(true);
	G2.setLabel("G2");
	G2.setShowInfo(true);
	NodeMap Assocs1(&S, &G1, "a1:A a2:A");
	CophyMap M1(Assocs1);
	M1.doPageReconciliation();
	Node* h = S.LCA("A", "C");
	Node* p = G1.LCA("a1", "a2");
	M1.moveToHost(p,h);
	M1.inferEvents();
	NodeMap Assocs2(&S, &G2, "c1:C c2:C");
	CophyMap M2(Assocs2);
	M2.doPageReconciliation();
	p = G2.LCA("c1", "c2");
	M2.moveToHost(p,h);
	M2.inferEvents();
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G1.getLabel() << endl << G1;
	cout << "Gene tree " << G2.getLabel() << endl << G2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts for " << G2.getLabel() << ":" << endl << E2 << endl;
	CophyMultiMap MCM;
	MCM.addCophyMap(&M1);
	MCM.addCophyMap(&M2);
	cout << " Test case 6 event counts: " << MCM.countEvents() << endl << hline;
}
void doTestCase7() {
	cout << "Test case 7: a bit more complicated" << endl;
	Tree S("(A,(B,C))");
	S.setLabel("S");
	Tree G1("(a1,((a2,(b,c1)), c2))");
	Tree G2("(a,(c1,(c2,c3)))");
	G1.setLabel("G1");
	G1.setShowInfo(true);
	G2.setLabel("G2");
	G2.setShowInfo(true);
	cout << G1 << G2 << endl;
	NodeMap Assocs1(&S, &G1, "a1:A a2:A b:B c1:C c2:C");
	CophyMap M1(Assocs1);
	M1.doPageReconciliation();
	M1.inferEvents();
	NodeMap Assocs2(&S, &G2, "a:A c1:C c2:C c3:C");
	CophyMap M2(Assocs2);
	M2.doPageReconciliation();
	M2.inferEvents();
	CophyMultiMap MCM;
	MCM.addCophyMap(&M1);
	MCM.addCophyMap(&M2);
	cout << "G1 = " << endl << G1;
	cout << "G2 = " << endl << G2;
	cout << "After standard reconciliation, event counts: " << MCM.countEvents() << endl;

	Node* h = S.LCA("A", "C");
	Node* p = G1.LCA("a1", "a2");
	M1.moveToHost(p,h);
	M1.inferEvents();
	p = G2.LCA("a", "c2");
	M2.moveToHost(p, h, duplication);
	p = G2.LCA("c1", "c2");
	M2.moveToHost(p, h, duplication);
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G1.getLabel() << endl << G1;
	cout << "Gene tree " << G2.getLabel() << endl << G2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts for " << G2.getLabel() << ":" << endl << E2 << endl;
	cout << "Duplication height at root of tree, should be 2: " << MCM.calcDuplicationHeight(h) << endl;
	cout << "Duplication height at (B,C), should be 0: " << MCM.calcDuplicationHeight(S.LCA("A","B")) << endl;
	cout << " Test case 7 event counts: " << MCM.countEvents() << endl << hline;
}
void ybcTestCases() {
//	doTestCase1();
	doTestCase2();
	doTestCase3a();
//	doTestCase3();
	doTestCase4();
	doTestCase5();
	doTestCase6();
//	doTestCase7();
}

void doAlgorithmTest() {
	CophyMultiMap CMM;
	Tree S('s', "(A,(B,C))");
	S.setLabel("S");
	Tree G1('g', "(a1,((b1,c1),(c2,(c3,(c4,c5)))))");
	Tree G2('g', "(a2,(c6,(c7,(c8,(c9,b2)))))");
	G1.setLabel("G1");
	G2.setLabel("G2");
	cout << "S:\n" << S << "G1:\n" << G1 << "G2:\n" << G2 << endl;
	NodeMap Assocs1(&S, &G1, "a1:A b1:B c1:C c2:C c3:C c4:C c5:C");
	CophyMap M1(Assocs1);
	G1.setInfo(M1.getInfo());
	G1.setShowInfo(true);
	NodeMap Assocs2(&S, &G2, "a2:A b2:B c6:C c7:C c8:C c9:C");
	CophyMap M2(Assocs2);
	G2.setInfo(M2.getInfo());
	G2.setShowInfo(true);
	CMM.addCophyMap(&M1);
	CMM.addCophyMap(&M2);
	map<string, int> sampledDistribution;
	Algorithm1(CMM, sampledDistribution);
}

void doAlgorithmTest2() {
	CophyMultiMap CMM;
	Tree S('s', "(A,(B,C))");
	S.setLabel("S");
	Tree G1('g', "(a1,(a2,(a3,a4)))");
	Tree G2('g', "(b1,(b2,(b3,b4)))");
	Tree G3('g', "(c1,(c2,(c3,c4)))");
	G1.setLabel("G1");
	G2.setLabel("G2");
	G3.setLabel("G3");
	cout << "S:\n" << S << "G1:\n" << G1 << "G2:\n" << G2 << "G3:\n" << G3 << endl;
	NodeMap Assocs1(&S, &G1, "a1:A a2:A a3:A a4:A");
	CophyMap M1(Assocs1);
	G1.setInfo(M1.getInfo());
	G1.setShowInfo(true);
	NodeMap Assocs2(&S, &G2, "b1:B b2:B b3:B b4:B");
	CophyMap M2(Assocs2);
	G2.setInfo(M2.getInfo());
	G2.setShowInfo(true);

	NodeMap Assocs3(&S, &G3, "c1:C c2:C c3:C c4:C");
	CophyMap M3(Assocs3);
	G3.setInfo(M3.getInfo());
	G3.setShowInfo(true);
	CMM.addCophyMap(&M1);
	CMM.addCophyMap(&M2);
	CMM.addCophyMap(&M3);
	lossCost = 0.1;
	duplicationCost = 1.0;
	map<string, int> sampledDistribution;
	Algorithm1(CMM, sampledDistribution);
}

void doAlgorithmTest3() {
	CophyMultiMap CMM;
	cout << "Test case 1 (trivial match)" << endl;
	Tree T1("(A,(B,C))");
	Tree S1("(a,(b,c))");
	T1.setLabel("T1");
	S1.setLabel("S1");
	S1.setShowInfo(true);
	NodeMap assocs1(&T1, &S1, "a:A b:B c:C");
	CophyMap M1(assocs1);
	S1.setInfo(M1.getInfo());	// XXX TODO The code breaks (possibly on displaying the map) if info is not set.
	CMM.addCophyMap(&M1);
	cout << CMM;
	lossCost = 0.1;
	duplicationCost = 1.0;
	map<string, int> sampledDistribution;
	Algorithm1(CMM, sampledDistribution);
}

void doContenderTest() {
	/**
	 * the Contender object should be stored in decreasing order of score in containers: checking this!
	 */
	set<Contender> C;
	Node p("p"), h("h");
	CophyMap M;
	for (uint i(0); i < 100; ++i) {
		Contender c(dran(1.0), &p, &h, noevent, &M);
		C.insert(c); // should make a copy
	}
	for (auto c : C) {
		cout << c.getScore() << " ";
	}
	cout << endl;
}
/**
 * DONE count events on a single cophymap XXX test
 * DONE sum events across multiple cophymaps ("Reconciliation")
 * DONE move mapped nodes within their feasible range (cannot move something above the image of its ancestor, for example)
 * TODO move all nodes above the root
 * TODO do minimal cost cophymap (image of parents mapped to LCA of images of all children)
 */

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
		"\t-d <float>\n\t\tto set the duplication event cost (default: " + to_string(duplicationCost) + ")\n"
		"\t-l <float>\n\t\tto set the loss event cost (default: " + to_string(lossCost) + ")\n"
		"\t-o (samples)\n\t\tto show the sampled maps (default: FALSE).\n"
		"\t--sat-spread <float>\n\t\tto set the Simulated Annealing \"spread\" parameter (default: " + to_string(SATempSpread) + ")\n"
		"\t--sat-decay <float>\n\t\tto set the Simulated Annealing \"decay\" parameter (default: " + to_string(SATempDecay) + ")\n"
		"\t--verbose\n\t\tto output lots of stuff (default: FALSE);\n"
		"\t--silent\n\t\tto output as little as possible (default: FALSE);\n"
	);

int main(int argn, char** argv) {
	if (argn <= 1) {
		cout << segdupHelp << endl;
		return 0;
	}
	unsigned ranseed = std::chrono::system_clock::now().time_since_epoch().count();
	generator.seed(ranseed);

	CophyMultiMap CMM;
	vector<CophyMap*> M;
	uint numGeneTrees(0);
	Tree *S(nullptr);
	vector<Tree*> G;
	for (int i(1); i < argn; ++i) {
		if (!strcmp(argv[i], "?") || !strcmp(argv[i], "-h")) {
			cout << segdupHelp;
			return 0;
		}
		if (!strcmp(argv[i], "-S")) {
			++i;
			string newick(argv[i]);
			S = new Tree('s', newick);
			S->setLabel("S");
			if (_verbose) {
				cout << "Input Species tree:" << endl << (*S) << endl;
			}
		} else if (!strcmp(argv[i], "-G")) {
			++numGeneTrees;
			++i;
			string newick(argv[i]);
			G.push_back(new Tree('g', newick));
			Tree *P = G[numGeneTrees-1];
			P->setLabel("G" + to_string(numGeneTrees));
			if (_verbose) {
				cout << "Input Gene tree:" << endl << (*P) << endl;
			}
			++i;
			string assoc(argv[i]);
			NodeMap* A = new NodeMap(S, P, assoc);
			if (_verbose) {
				cout << "Input Associations:" << endl << (*A);
			}
			CophyMap* M = new CophyMap(*A);
			P->setInfo(M->getInfo());
			P->setShowInfo(true);
			CMM.addCophyMap(M);
		} else if (!strcmp(argv[i], "-n")) {
			++i;
			nSteps = atoi(argv[i]);
			if (_verbose) {
				cout << "Setting Number of steps to " << nSteps << endl;
			}
		} else if (!strcmp(argv[i], "-d")) {
			++i;
			duplicationCost = atof(argv[i]);
			if (_verbose) {
				cout << "Setting DuplicationCost to " << duplicationCost << endl;
			}
			CMM.setDuplicationCost(duplicationCost);
		} else if (!strcmp(argv[i], "-l")) {
			++i;
			lossCost = atof(argv[i]);
			CMM.setLossCost(lossCost);	// TODO Settle on either a global variable for this cost or just the instance variable!
			if (_verbose) {
				cout << "Setting LossCost to " << lossCost << endl;
			}
		} else if (!strcmp(argv[i], "-Tinit")) {
			++i;
			Tinitial = atof(argv[i]);
			if (_verbose) {
				cout << "Setting Initial Temperature to " << Tinitial << endl;
			}
		} else if (!strcmp(argv[i], "-o")) {
			++i;
			if (!strcmp(argv[i], "probs")) {
				_outputProbabilities = true;
			} else if (!strcmp(argv[i], "samples")) {
				if (_verbose) {
					cout << "Setting Show_Samples to true" << endl;
				}
				_showSampledDistribution = true;
			} else if (!strcmp(argv[i], "trace")) {
				if (_verbose) {
					cout << "Setting Save a Trace to true" << endl;
				}
				_saveTrace = true;
			}
		} else if (!strcmp(argv[i], "--sat-spread")) {
			++i;
			SATempSpread = atof(argv[i]);
			if (_verbose) {
				cout << "Setting Simulated Annealing \"spread\" parameter to " << SATempSpread << endl;
			}
		} else if (!strcmp(argv[i], "--sat-decay")) {
			++i;
			SATempDecay = atof(argv[i]);
			if (_verbose) {
				cout << "Setting Simulated Annealing \"decay\" parameter to " << SATempDecay << endl;
			}
		} else if (!strcmp(argv[i], "--verbose")) {
			_verbose = true;
		} else if (!strcmp(argv[i], "--silent")) {
			_silent = true;
		}
	}
	CMM.doPageReconciliation();
	map<string, int> sampledDistribution;
	Algorithm1(CMM, sampledDistribution);
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

