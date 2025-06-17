/*
 * Tests.cpp
 *
 *  Created on: 17 June 2025
 *      Author: mac
 */

#include <set>

#include "../utility/debugging.h"
#include "../utility/myrandom.h"

#include "Contender.h"
#include "CophyMap.h"
#include "CophyMultiMap.h"
#include "DupMove.h"
#include "Node.h"
#include "SegDup.h"
#include "SingleNodeMove.h"
#include "Tests.h"
#include "Tree.h"

using namespace segdup;
using namespace std;

extern string hline;
extern double lossCost;
extern double duplicationCost;

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
