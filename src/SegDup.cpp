/*
 * SegDup.cpp
 *
 *  Created on: 2 Aug 2022
 *      Author: mac
 */

#include <map>
#include <stdio.h>

#include "Node.h"
#include "Tree.h"
#include "CophyMap.h"
#include "CophyMultiMap.h"

#include "../utility/debugging.h"

/**
 * Input to an instance includes whether each node is mapped to a D or S
 */

using namespace std;
using namespace segdup;

bool _debugging(false);
bool _silent(false);

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
//	br->setAssociate(T["r"]);
//	bx->setAssociate(x);
//	bc->setAssociate(c);
//	by->setAssociate(a);
//	ba1->setAssociate(a);
//	ba2->setAssociate(a);
//	bb->setAssociate(b);
	cout << G1 << endl;

//	Node* g2 = new Node("redRoot");
//	cout << "all done" << endl;
//	cout << "depth of node a=" << a->getLabel() << " is " << a->getDepth() << endl;
//	T.gatherVertices();
//	for (auto v : T.getVertices()) {
//		cout << "Vertex " << v.first << " has depth " << v.second->getDepth() << endl;
//	}
//	T.compressTraverseWrite(cout);

	G1.compressTraverseWrite(cout);
//	NodeMap leafAssociation({("a","b")});
//	CophyMap M(&T, &G1);
//	M.doPageReconciliation();
//	nodemap* phi = M.getPhi();
//	for (auto iter = phi->begin(); iter != phi->end(); ++iter) {
//		cout << iter->first->getLabel() << "\t->\t" << iter->second.first->getLabel() << ' ';
//		if (iter->second.second == 0) {
//			cout << "codivergence";
//		} else {
//			cout << "duplication";
//		}
//		cout << endl;
//	}
	G1.compressTraverseWrite(cout);
	// TODO Create a reconciliation and count the events on it.
	// TODO helper methods:
	// TODO read a Newick format tree?
//	cout << "Ancestry of nodes:" << endl;
//	G1.calcAncestry();
//	for (auto upper : G1.getVertices()) {
//		cout << "first vertex = " << upper.first << endl;
//		if (upper.second->isLeaf()) {
//			continue;
//		}
//		for (auto lower : G1.getVertices()) {
//			cout << "\tsecond vertex = " << lower.first << endl;
//			if (G1.isAncestralTo(upper.second, lower.second)) {
//				cout << upper.first << " is ancestral to " << lower.first << " and is " << G1.getDistUp(lower.second, upper.second) << " edge(s) above it" << endl;
//			}
//		}
//	}
	Tree H("(a,(b,c));");
	H.getRoot()->setLabel("h0");
	Tree P("(p,(q,r))");
	std::map<std::string, Node*>& VH = H.getVertices();
	std::map<std::string, Node*>& VP = P.getVertices();
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

void testPerfectMatchReconciliation() {
	Tree H("(((a,b),c),(d,e))");
	H.getRoot()->setLabel("h0");
	Tree P("(((a,b),c),(d,e))");
	P.getRoot()->setLabel("p0");
	cout << "P = " << endl << P;
	cout << "H = " << endl << H;
	NodeMap Assocs(&H, &P, "a:a, b:b, c:c, d:d, e:e");
	CophyMap M(Assocs);
	M.doPageReconciliation();
	P.setShowInfo(true);
	cout << "P = " << endl << P;
	cout << "H = " << endl << H;
	EventCount E = M.countEvents();
	cout << "Event counts: " << E << endl;
}

void ybcTestCases() {
	string hline("============================================================================\n");
	cout << "Test case 1 (trivial match)" << endl;
	Tree T1("(A,(B,C))");
	Tree S1("(a,(b,c))");
	T1.setLabel("T1");
	S1.setLabel("S1");
	S1.setShowInfo(true);
	NodeMap assocs1(&T1, &S1, "a:A b:B c:C");
	CophyMap M1(assocs1);
	M1.doPageReconciliation();
	M1.storeHostInfo();
	cout << M1 << M1.countEvents() << endl << hline;

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
	G.setShowInfo(true);
	cout << "G:" << endl << G;
	EventCount E = M.countEvents();
	cout << "Event counts: " << E << endl << hline;

	cout << "Test case 3: LCA, 2 gene trees" << endl;
	Tree G2("((a1,a2),((b,c1),c2)))");
	G2.setLabel("G2");
	G2.setShowInfo(true);
	NodeMap Assocs2(&S, &G2, "a1:A, a2:A, b:B, c1:C, c2:C");
	CophyMap M2(Assocs2);
	M2.doPageReconciliation();
	cout << "Species tree " << S.getLabel() << endl << S;
	cout << "Gene tree " << G.getLabel() << endl << G;
	cout << "Gene tree " << G2.getLabel() << endl << G2;
	EventCount E2 = M2.countEvents();
	cout << "Event counts for " << G2.getLabel() << ":" << endl << E2 << endl;
	E += E2;
	cout << "Total event counts = " << E << endl << hline;

	CophyMultiMap MCM;
	MCM.addCophyMap(&M);
	MCM.addCophyMap(&M2);
}

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
	M.storeHostInfo();
	cout << M;
	M.calcAvailableNewHosts(G["a"]);
	cout << G;
	cout << M;
	cout << M.countEvents() << endl;
	std::set<Node*> avail = M.calcAvailableNewHosts(G["p"]);
	for (Node* n : avail) {
		cout << "available host : " << *n << endl;
	}
	M.moveToHost(G["p"], H["X"]);
	M.inferEvents();
	M.storeHostInfo();
	cout << M;
	cout << M.countEvents() << endl;
}

void testDuplicationHeight() {
	Tree S("(((A,B),(C,D)),E)");
	Tree G("((((a1,b1),((a2,b2),(a3,b3))),(((c1,d1),c2),d2)),e)");
	cout << S;
	cout << G;
	NodeMap A1(&S, &G, "a1:A a2:A a3:A b1:B b2:B b3:B c1:C c2:C d1:D d2:D e:E");
	CophyMap M(A1);
	M.doPageReconciliation();
	M.storeHostInfo();
	G.setShowInfo(true);
	cout << M << M.countEvents();
	Tree G2("(((((a1,b1),(a4,b4)),((a2,b2),(a3,b3))),((c1,d1),c2)),e)");
	NodeMap A2(&S, &G2, "a1:A a2:A a3:A a4:A b1:B b2:B b3:B b4:B c1:C c2:C d1:D e:E");
	CophyMap M2(A2);
	M2.doPageReconciliation();
	M2.storeHostInfo();
	G2.setShowInfo(true);
	cout << M2 << M2.countEvents() << endl;
	CophyMultiMap CMM;
	CMM.addCophyMap(&M);
	CMM.addCophyMap(&M2);
	int h = CMM.calcDuplicationHeight(S["v4"]);
	cout << "dh(v4) = " << h << endl;
	h = CMM.calcDuplicationHeight(S["v3"]);
	cout << "dh(v3) = " << h << endl;
	h = CMM.calcDuplicationHeight(S["v5"]);
	cout << "dh(v5) = " << h << endl;
}

void Algorithm1(CophyMultiMap* M) {

}

/**
 * TODO count events on a single cophymap XXX test
 * TODO sum events across multiple cophymaps ("Reconciliation")
 * TODO move mapped nodes within their feasible range (cannot move something above the image of its ancestor, for example)
 * TODO move all nodes above the root
 * TODO do minimal cost cophymap (image of parents mapped to LCA of images of all children)
 */
int main(int argc, char** argv) {
//	bool _debugging(true);

//	testPageReconciliation2();
//	testPageReconciliation3();
//	testPerfectMatchReconciliation();
//	testAvailableHosts();
	testDuplicationHeight();
//	ybcTestCases();

	return 0;
}

