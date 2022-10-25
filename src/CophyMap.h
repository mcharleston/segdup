/*
 * CophyMap.h
 *
 *  Created on: 10 Oct 2022
 *      Author: mac
 */

#ifndef SRC_COPHYMAP_H_
#define SRC_COPHYMAP_H_

#include <map>
#include "Node.h"
#include "Tree.h"
#include "NodeMap.h"
#include "EventCount.h"

namespace segdup {

/**
 * A nodemap maps nodes (in a parasite/gene tree) to (node, associationtype) pairs where the nodes are in a
 * host/species tree and the associationtype is
 * 	0 for "on the actual vertex" and
 * 	1 for "on the parent edge".
 * Thus nodemap Phi[x] = (a,1) means x is mapped to the edge above host node a in H so has a duplication event
 * plus 0 or more losses; and Phy[x] = (a,0) means x is mapped to the hode node a in H, so has a codivergence
 * event associated with it (and 0 or more losses in its descendant lineages).
 */

class CophyMap {
private:
	Tree* H;	// the species or host tree
	Tree* P;	// the gene or parasite tree
	NodeMap phi;	// nodes in P mapped to nodes in H -- starting with the leaves
	inversenodemap invPhi;
	// Note that in genetree/species tree problems we talk about a gene tree G and a species tree S.
public:
	CophyMap() : H(nullptr), P(nullptr) {}
	CophyMap(Tree *host, Tree *para) : H(host), P(para) {}
	CophyMap(NodeMap inphi);
	virtual ~CophyMap();
	CophyMap(const CophyMap &other);

	std::set<Node*> calcAvailableNewHosts(Node* p);
//	int calcDuplicationHeight(const std::string& str);

	EventCount countEvents();
	void doPageReconciliation();

	NodeMap& getPhi() { return phi; }
	Node* getHost(Node* p) { return phi[p].first; }
	Tree* getHostTree() { return H; }
	Tree* getParasiteTree() { return P; }
	short getType(Node* p) { return phi[p].second; }

	void inferEvents();
	void inferEvents(Node *p);
	bool isValid();

	Node* mapToLCAofChildren(Node* p);
	void moveToHost(Node* p, Node *h);

	CophyMap& operator=(const CophyMap &other);

	void setPhi(std::string pstr, std::string hstr);
	void setPhi(Node* p, Node* h);
	void storeHostInfo();

};

std::ostream& operator<<(std::ostream& os, CophyMap& M);

} /* namespace segdup */

#endif /* SRC_COPHYMAP_H_ */
