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
bool _saveTrace(false);
bool _saveFinal(false);
bool _showSampledDistribution(false);
int outputInterval(1);
int nSteps(1000);
double Tinitial(10.0);
double Tfinal(0.0);
int nFinal(0);

extern std::mt19937 generator;
extern unsigned seed;

std::ofstream summaryfile;

double lossCost(defLossCost);			// XXX
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
	int h = CMM.calcCombinedDuplicationHeight(S["v4"]);
	cout << "dh(v4) = " << h << endl;
	h = CMM.calcCombinedDuplicationHeight(S["v3"]);
	cout << "dh(v3) = " << h << endl;
	h = CMM.calcCombinedDuplicationHeight(S["v5"]);
	cout << "dh(v5) = " << h << endl;
	cout << "Counting events with segmental duplications" << endl << CMM.countEvents() << endl;
}

//void Algorithm1(CophyMultiMap& CMM, map<string, int>& sampledDistribution) {
//	bool _debugging(true);
//	DEBUG(cout << hline << "Algorithm1" << endl << hline);
//	DEBUG(cout << "Input multi-map:" << endl << CMM);
//	CMM.doPageReconciliation();
//	DEBUG(cout << "Initial reconciliation complete:" << endl << CMM);
//
//	CophyMultiMap nuMM(CMM);
//
//	vector<pair<Node*,CophyMap*>> allMoveableNodes;
//	CMM.putAllMoveableNodes(allMoveableNodes);
//
//	double T(0.0);
//	std::set<Contender> neighbours;
//	string mapDescription;
//	CMM.toCompactString(mapDescription);
////	cout << "ORIGINAL MAP:" << endl << CMM << "Events\tScore\tMap\n" << ecoriginal << '\t' << CSD(ecoriginal) << '\t' << mapDescription << endl;
//	string bestMMap;
//	EventCount bestEventCount;
//	double bestCost(1e100);	// 10^100 should be enough!!
//	ofstream ftrace;
//	int sampleNumber(0);
//	if (_saveTrace) {
//		ftrace.open("segdup-trace.csv", std::ofstream::out);
//		ftrace << "i,c,d,l,s,pr" << endl;
//	}
//	ostringstream bestPrettyMap;
//	unsigned int numNeighbours(0);
//	EventCount currentEventCount = CMM.countEvents();
//
//	/**
//	 * Choose a node p at uniform random from all internal nodes of all gene trees.
//	 * Collect all nodes from all maps into one big array
//	 * Given p, calculate all potential new locations and note "no change" option too
//	 * Sample from that set according to the relative scores
//	 */
//	double minScore(0.0);
//	for (int t(1); t <= nSteps; ++t) {
//		minScore = 1e10;	// start with a large number
//		if (t % 10 == 0) {
//			advance_cursor();
//		}
//		if (t % 100 == 0) {
//			update_message(" steps = " + to_string(t) + "/" + to_string(nSteps)
//					+ "; log score = "
//					+ to_string(-CSD(currentEventCount) / (1.0*T)) + "; ec = "
//					+ to_string(currentEventCount.codivs) + "C + "
//					+ to_string(currentEventCount.dups) + "D + "
//					+ to_string(currentEventCount.losses) + "L; cost = "
//					+ to_string(CSD(currentEventCount))
//					+ "; best cost = " + to_string(bestCost));
//			DEBUG(cout << endl << CMM);
////			break;
//		}
//		int nullMoves(0);
//		EventCount ec;
//		T = (Tinitial-Tfinal)*(1.0 - (1.0 * (t-1.0) / nSteps)) + Tfinal;
////		double x = 1.0 * t / nSteps;
////		T = Tinitial*1.5/(1.0 + 0.2 * sqrt(x*x*t+1.0));	// function 4 / sqrts	// XXX looks good on ybc-case3 and on Guido et al..
//		neighbours.clear();
//		numNeighbours = 0;
//		DEBUG(
//			cout << hline << "currentEventCount = " << currentEventCount << endl << hline;
//		);
//#define ChooseByNode
//#ifdef ChooseByNode
//		auto mpr = allMoveableNodes[iran(allMoveableNodes.size())];	// pick at random from allMoveableNodes
//		Node* p(mpr.first);
//		CophyMap *M = mpr.second;
//#else
//		for (auto mpr : CMM.getMaps()) {	// get all maps of the gene trees into single species tree
//			CophyMap* M = mpr.second;
//			for (Node* p : iV[M]) {
//#endif
//				set<pair<Node*, eventType>> nextImages = M->calcAvailableNewHosts(p);
//				DEBUG(
//					cout << "Node " << p->getLabel() << " has possible images { ";
//					for (auto a : nextImages) {
//						cout << eventSymbol[a.second] << a.first->getLabel() << " ";
//					}
//					cout << "};" << endl;
//				);
//				for (auto a : nextImages) {
//					Node* nuHost = a.first;
//					eventType nuEvent = a.second;
//					if (nuHost == M->getHost(p) && nuEvent == M->getEvent(p)) {
//						double relativeSamplingProbability = exp(-CSD(currentEventCount) / (1.0*T));
////						cerr << t << ',' << T << ',' << p << ',' << CSD(currentEventCount) << endl;
//						Contender noChange( currentEventCount, CSD(currentEventCount), p, nuHost, nuEvent, M );
////						Contender noChange( currentEventCount, relativeSamplingProbability, p, nuHost, nuEvent, M );
////						DEBUG(cout << "XXX entering CSD for no move: " << CSD(currentEventCount) << endl);
//						minScore = min(minScore, CSD(currentEventCount));
//
//						noChange._noMove = true;
//						noChange.setLabel("leaving " + p->getLabel() + " on [" + eventSymbol[nuEvent] + "]" + nuHost->getLabel());
//
//						neighbours.insert(noChange);
//						++nullMoves;
//						continue;
//					}
//					ec.clear();
//					Node* oldHost = M->getHost(p);
//					eventType oldEvent = M->getEvent(p);
//					// Get cost of move, without recalculating everything:
//
//					// Any change in the number of codivergences?
//					DEBUG(cout << "p=" << p->getLabel() << "; oldHost=" << oldHost->getLabel()
//							<< "; newHost=" << a.first->getLabel() << endl);
//					ec.codivs = (nuEvent == codivergence) ? 1 : 0;
//					ec.codivs -= (oldEvent == codivergence) ? 1 : 0;
//
//					// Any change in number of losses?
//					int lossFactor = (p->hasParent()) ? 1 : 2;	// root of gene tree has no ancestral branch!
//					if (oldHost->isAncestralTo(nuHost)) {
//						DEBUG(cout << "old host is ancestral to new host" << endl);
//						ec.losses = -lossFactor*(nuHost->getTree()->getDistUp(nuHost, oldHost));
////						if (p->hasParent()) {
////							ec.losses = -(nuHost->getTree()->getDistUp(nuHost, oldHost));
////						} else {
////							ec.losses = -2*(nuHost->getTree()->getDistUp(nuHost, oldHost));
////						}
//					} else if (nuHost->isAncestralTo(oldHost)) {
//						DEBUG(cout << "new host is ancestral to old host" << endl);
//						ec.losses = lossFactor*(nuHost->getTree()->getDistUp(oldHost, nuHost));
////						if (p->hasParent()) {
////							ec.losses = oldHost->getTree()->getDistUp(oldHost, nuHost);
////						} else {
////							ec.losses = 2*oldHost->getTree()->getDistUp(oldHost, nuHost);
////						}
//					} else {
//						DEBUG(cout << "old and new hosts are not ancestrally comparable" << endl);
//					}
//					if (oldEvent == duplication && nuEvent == codivergence) {
//						ec.losses -= 2;
//					} else if (nuEvent == duplication && oldEvent == codivergence) {
//						ec.losses += 2;
//					}
//
//					// Any change in number of duplications?
//					int oldJointDupHeightHere(CMM.calcCombinedDuplicationHeight(oldHost));
//					int oldJointDupHeightThere(CMM.calcCombinedDuplicationHeight(nuHost));
//
//
//					DEBUG(cout << "oldJointDupHeightHere = " << oldJointDupHeightHere << endl);
//					DEBUG(cout << "oldJointDupHeightThere = " << oldJointDupHeightThere << endl);
//					M->moveToHost(p, nuHost, nuEvent);
//					CMM.movePToHost(p,oldHost,nuHost);	// XXX are both these function calls necessary?
//					int nuJointDupHeightHere(CMM.calcCombinedDuplicationHeight(oldHost));
//					int nuJointDupHeightThere(CMM.calcCombinedDuplicationHeight(nuHost));
//					DEBUG(cout << "nuJointDupHeightHere = " << nuJointDupHeightHere << endl);
//					DEBUG(cout << "nuJointDupHeightThere = " << nuJointDupHeightThere << endl);
//					M->moveToHost(p, oldHost, oldEvent);
//					CMM.movePToHost(p, nuHost, oldHost);	// XXX are both these function calls necessary?
//					// return p to its old place
////					DEBUG(
////							cout << "currentJointDupHeightHere = " << currentJointDupHeightHere << endl;
////							cout << "currentJointDupHeightThere = " << currentJointDupHeightThere << endl;
////							cout << "nuJointDupHeightHere = " << nuJointDupHeightHere << endl;
////							cout << "nuJointDupHeightThere = " << nuJointDupHeightThere << endl;
////							);
//					ec.dups = nuJointDupHeightThere - oldJointDupHeightThere;
//					if (nuHost != oldHost) {
//						ec.dups += nuJointDupHeightHere - oldJointDupHeightHere;
//					}
////					int oldHostDuplicationHeight(M->calcDuplicationHeight(p)); // the dup height of this p on its original host
////					int destinationDuplicationHeight(CMM.calcCombinedDuplicationHeight(nuHost));	// the current JOINT dup height on the destination
////					M->moveToHost(p, nuHost, nuEvent);
////					// also need new JOINT dup heights
////					// TODO think of a nice way that I don't have to move p twice..
////					int dupHeightOfPOnNewHost = M->calcDuplicationHeight(p);
////					ec.dups = 0;
////					if (dupHeightOfPOnNewHost > destinationDuplicationHeight) {
////						ec.dups = 1;	// one more duplication required (can only go up by one)
////					}
////					int remainingHostDupHeightOnOldHost(CMM.calcCombinedDuplicationHeight(oldHost));
////					if (oldHostDuplicationHeight > remainingHostDupHeightOnOldHost) {
////						ec.dups -= 1;	// dup height on original host/species was due to this parasite/gene
////					}
//
//					// Store this contender and its eventcount:
////					DEBUG(
////					);
//					EventCount dec(ec);
//					ec += currentEventCount;
//					if (ec.dups < 0) {
//						cout << "step " << t << ": CRITICAL FAILURE!" << endl;
//						cout << "previous event count = " << currentEventCount << endl;
//						cout << "delta ec = " << dec << ";\t";
//						cout << "ec+originalEventCount = " << ec << ";\t";
//						exit(-1);
//					}
//					Association ass(p, nuHost, nuEvent);
////					double relativeSamplingProbability = exp(-CSD(ec) / (1.0*T));	// XXX test just using dec not ec here
//					minScore = min(minScore, CSD(ec));
////					DEBUG(
////							cout << "SCORE = " << relativeSamplingProbability << endl;
////					);
////					Contender con( ec, relativeSamplingProbability, p, nuHost, nuEvent, M );
//					Contender con( ec, CSD(ec), p, nuHost, nuEvent, M );
////					DEBUG(cout << "XXX setting contender score to " << CSD(ec) << endl);
//					if (con.getEventCount().dups < 0) {
//						cout << "step " << t << ": CRITICAL FAILURE!" << endl;
//						cout << "Contender event count = " << con.getEventCount().dups << endl;
//						exit(-1);
//					}
//					string label = "moving " + p->getLabel() + " to [" + eventSymbol[nuEvent] + "]" + nuHost->getLabel();
//					con.setLabel(label);
//					neighbours.insert(con);
////					cerr << t << ',' << T << ',' << p << ',' << con.getScore() << endl;
//#ifdef ChooseByNode
//#else
//				}
//			}
//#endif
//		}
//		double d;
//		_debugging = true;
//		double total(0.0);
//		set<Contender> adjustedNeighbours;
//		for (auto nei : neighbours) {
//			d = nei.getScore() - minScore; // XXX This is a fudge...
////			cerr << "d=" << d << endl;
//			DEBUG(cout << "relative CSD for this contender = " << d << endl);
//			nei.setScore(exp(-d / (1.0*T)));
//			total += nei.getScore();
//			adjustedNeighbours.insert(nei);
//			DEBUG(cout << "T = " << T << "; new score = " << nei.getScore() << endl);
//			DEBUG(cout << "new total " << total << endl);
//		}
////		for (auto nei : neighbours) {
////		}
//		numNeighbours += neighbours.size();
//		DEBUG(cout << "total before dran() = " << total << endl);
//		double r = dran(total);
//		DEBUG(cout << "total after dran() = " << total << endl);
//		/***********************************
//		 * Now the sampling!
//		 ***********************************/
//		DEBUG(cout << "initial r = " << r << " from U[0, " << total << "]" << endl);
//		for (auto nei : adjustedNeighbours) {
//			if (r <= nei.getScore()) {
//				DEBUG(cout << "r=" << r << "; score=" << nei.getScore() << endl);
//				if (nei._noMove) {
//					DEBUG(cout << "Selected move: No change (probability = " << (nei.getScore()/total) << ")" << endl);
//				} else {
//					// sample this one
//					DEBUG(cout << "Selected move: " << nei.getLabel() << " (rel. probability = " << nei.getScore() << ")" << endl);
//					nei.getMap()->moveToHost(nei.getParasite(), nei.getHost(), nei.getEvent());
//					DEBUG(nei.getMap()->checkValidHostOrdering());
//				}
//				if (_showSampledDistribution) {
//					CMM.toCompactString(mapDescription);
//#define UseLongMapDescription
//#ifdef UseLongMapDescription
//					EventCount totalEC = CMM.countEvents();
//					mapDescription += "-D" + to_string(totalEC.dups) + "L" + to_string(totalEC.losses);
//#else
//					mapDescription += "-D" + to_string(nei.getEventCount().dups)
//							+ "L" + to_string(nei.getEventCount().losses);
//#endif
//					sampledDistribution[mapDescription] += 1;
//				}
//				if (_cacheEventCounts) {
//					CMM.toCompactString(mapDescription);
//					sampledDistribution[mapDescription] += 1;
//					ec = CMM.getEventCount(mapDescription);
//				} else {
//					ec = nei.getEventCount();
//				}
//				DEBUG(cout << "Retrieving event count from neighbour: " << ec << endl);
//				DEBUG(cout << CMM);
////				currentEventCount = ec;
//				currentEventCount = CMM.countEvents();
//				double cost = CSD(ec);
//				if (cost < bestCost) {
//					bestCost = cost;
//					DEBUG(cout << "Best cost = " << bestCost << endl);
//					bestEventCount = ec;
//					DEBUG(if (bestCost < 0) {
//						cout << "step " << t << ": best cost is NEGATIVE!" << endl;
//						cout << "\tbest cost = " << bestCost << endl;
//						cout << "\tbest event count = " << ec << endl;
//						bestPrettyMap.str("");
//						bestPrettyMap << CMM;
//						exit(-1);
//					});
//					bestMMap = mapDescription;
//					bestPrettyMap.str("");
//					bestPrettyMap << CMM;
//					DEBUG(
//							cout << CMM << endl;
//					);
//					break;
//				}
//				if (_saveTrace) {
//					++sampleNumber;
////					ftrace << sampleNumber << ",\"" << mapDescription << "\"," << to_string(CSD(ec)) << endl;
//					ftrace << sampleNumber << ',' << ec.codivs << ',' << ec.dups << ','
//							<< ec.losses << ',' << to_string(CSD(ec)) << ','
//							<< nei.getScore()
//							<< endl;
//				}
//				DEBUG(cout << nei.getLabel() << '\t' << nei.getScore() << endl);
//				break;
//			}
//			DEBUG(cout << r << ' ');
//			r -= nei.getScore();
//		}
//		DEBUG(cout << endl);
//		// XXX
//	}
//	cout << endl;
//	if (_showSampledDistribution) {
//		ofstream fout("segdup-samples.csv");
//		fout << hline << "Sampled Distribution of Solutions:" << endl
//				<< "Map-Events\tSamples" << endl;
//		for (auto dis : sampledDistribution) {
//			fout << dis.first << "\t" << dis.second << endl;
//		}
//		fout.close();
////		cout << hline << "Sampled Distribution of Solutions:" << endl
////				<< "Map-Events\tSamples" << endl;
////		for (auto dis : sampledDistribution) {
////			cout << dis.first << "\t" << dis.second << endl;
////		}
//	}
////	cout << "FINAL Multiple CophyMap found by Algorithm 1:" << endl;
////	EventCount ecFinal(CMM.countEvents());
////	cout << CMM << ecFinal << endl << hline << endl;
////	cout << "final CSD: " << CSD(ecFinal) << endl;
//	cout << hline << "BEST Multiple CophyMap found by Algorithm 1:" << endl;
//	cout << bestMMap << '\t' << endl;
//	cout << bestEventCount << '\t' << bestCost << endl;
////	cout << bestPrettyMap.str();
//	cout << hline << endl;
//	if (_saveTrace) {
//		ftrace.close();
//	}
//	summaryfile << bestEventCount.codivs << ',' << bestEventCount.dups << ',' << bestEventCount.losses << ',' << bestCost << endl;
//}

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

void Algorithm2(CophyMultiMap& CMM, vector<DupMove*> moves, vector<double> probs, map<string, int>* sampledDistribution = nullptr) {
	string bestMMap;
	EventCount bestEventCount;
	double bestCost(1e100);	// 10^100 should be enough!!
	string mapDescription;				

	bool _debugging(false);
	DEBUG(cout << hline << "Algorithm1" << endl << hline);
	//DEBUG(cout << "Input multi-map:" << endl << CMM);
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
		ftrace << "i,c,d,l,s,t" << endl;
	}
	ostringstream bestPrettyMap;
	CMM.calcEventCount();
	bestEventCount = CMM.countEvents();
	bestCost = CSD(bestEventCount);

	DEBUG(cout << "number of movable nodes: " << CMM.getAllMoveableNodes().size() << endl);

	if (CMM.getAllMoveableNodes().size() != 0) {
		for (int t(1); t <= nSteps + nFinal; ++t) {
			int nullMoves(0);
			EventCount ec;
			if (t <= nSteps)
				T = (Tinitial-Tfinal)*(1.0 - (1.0 * (t-1.0) / nSteps)) + Tfinal;
			else
				T = Tfinal;

			//if (t == 2752238) 

			SelectNextConfiguration(CMM, T, moves, probs);	// modifies CMM
													
			if (_showSampledDistribution && (!_saveFinal || t > nSteps)) {
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
	if (_showSampledDistribution) {
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

	//some ybc debugging
	
	/*CMM.calcEventCount();
	cout << CMM.countEvents() << endl;*/

	/*Tree* s = CMM.getHostTree();
	inversenodemap& invMap = CMM.getInverseMap();*/
	/*for (Node* p : invMap[s->getRoot()]) {
		CophyMap* M = CMM.getMap(p);
		cout << p->getLabel() << '\t' << M->getDuplicationHeight(p) << endl;
	}*/
	/*for (auto v : s->getVertices()) {
		cout << v.first << '\t' << CMM.calcCombinedDuplicationHeight(v.second) << endl;
	}*/

	if (_saveTrace) {
		ftrace.close();
	}
	summaryfile << bestEventCount.codivs << ',' << bestEventCount.dups << ',' << bestEventCount.losses << ',' << bestCost << endl;
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
	cout << "Duplication height at root of tree, should be 2: " << MCM.calcCombinedDuplicationHeight(h) << endl;
	cout << "Duplication height at (B,C), should be 0: " << MCM.calcCombinedDuplicationHeight(S.LCA("A","B")) << endl;
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
	vector<DupMove*> moves;
	moves.push_back(new SingleNodeMove());
	vector<double> probs;
	probs.push_back(1.0);
	Algorithm2(CMM, moves, probs);
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
	vector<DupMove*> moves;
	moves.push_back(new SingleNodeMove());
	vector<double> probs;
	probs.push_back(1.0);
	Algorithm2(CMM, moves, probs);
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
	vector<DupMove*> moves;
	moves.push_back(new SingleNodeMove());
	vector<double> probs;
	probs.push_back(1.0);
	Algorithm2(CMM, moves, probs);
}

void doContenderTest() {
	/**
	 * the Contender object should be stored in decreasing order of score in containers: checking this!
	 */
	set<Contender> C;
	Node p("p"), h("h");
	CophyMap M;
	EventCount ec;
	for (uint i(0); i < 100; ++i) {
		Contender c(ec, dran(1.0), &p, &h, noevent, &M);
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
		"\t-n <int>\n\t\tto supply the number of steps for Algorithm 1 (default value " + to_string(nSteps) + ")\n"
		"\t-d <float>\n\t\tto set the duplication event cost (default value " + to_string(duplicationCost) + ")\n"
		"\t-l <float>\n\t\tto set the loss event cost (default value " + to_string(lossCost) + ")\n"
		"\t-o (samples|trace|interval <int>|final)\n"
		"\t\tsamples to save the sampled distribution of maps (default value false)\n"
		"\t\ttrace to save a trace of the progress of the search (default value false)\n"
		"\t\tinterval to save only every few samples (default value 1)\n"
		"\t\tfinal to save only at the final temperature (default value false)\n"
		"\t-Tinit <float>\n\t\tto supply the initial temperature (default value " + to_string(Tinitial) + ")\n"
		"\t-Tfinal <float>\n\t\tto supply the final temperature (default value " + to_string(Tfinal) + ")\n"
		"\t-nfinal <float>\n\t\tto supply the number of steps at the final temperature (default value " + to_string(nFinal) + ")\n"
		"\t--seed <int>\n\t\tto set the random number generator seed\n"
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
			cout << "Input Species tree:" << endl << (*S) << endl;
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
			cout << "Input Gene tree:" << endl << (*P) << endl;
			++i;
			string assoc(argv[i]);
			NodeMap* A = new NodeMap(S, P, assoc);
			cout << "Input Associations:" << endl << (*A);
			CophyMap* M = new CophyMap(*A);
			P->setInfo(M->getInfo());
			P->setShowInfo(true);
			CMM.addCophyMap(M);
		} else if (!strcmp(argv[i], "-n")) {
			++i;
			nSteps = atoi(argv[i]);
			cout << "Setting Number of steps to " << nSteps << endl;
		} else if (!strcmp(argv[i], "-d")) {
			++i;
			duplicationCost = atof(argv[i]);
			cout << "Setting DuplicationCost to " << duplicationCost << endl;
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
				cout << "Setting Show_Samples to true" << endl;
				_showSampledDistribution = true;
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

