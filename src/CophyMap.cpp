/*
 * CophyMap.cpp
 *
 *  Created on: 10 Oct 2022
 *      Author: mac
 */

#include <assert.h>
#include <stdio.h>
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "CophyMap.h"
#include "EventCount.h"

using namespace std;

extern bool _debugging;

namespace segdup {

CophyMap::CophyMap(NodeMap inphi) {
	H = inphi.getHostTree();
	P = inphi.getParasiteTree();
	phi = inphi;
}

CophyMap::~CophyMap() {
}

CophyMap::CophyMap(const CophyMap &other) {
	H = other.H;
	P = other.P;
	phi = other.phi;
}

//int CophyMap::calcDuplicationHeight(const std::string& hstr) {
//	// for each parasite / gene on host node with this label, find the maximum run of p->parent->parent... on the same host
//	Node* h = (*H)[hstr];
//	return calcDuplicationHeight(h);
//}
//int CophyMap::calcDuplicationHeight(Node* h) {
//	// for each parasite / gene on host node with this label, find the maximum run of p->parent->parent... on the same host
//	int dupheight(0);
//	for (auto x : invPhi[h]) {
//		int dh = 0;
//		while ()
//	}
//}

std::set<Node*> CophyMap::calcAvailableNewHosts(Node* p) {
	// can go up to same vertex to which parent is mapped, and down to the LCA of there the children are mapped.
	bool _debugging(true);
	DEBUG(cout << "Calculating available hosts for node " << p->getLabel() << ":" << endl);
	std::set<Node*> avail;
	Node* bottom;
	Node* top = getHost(p->getParent()); //->getParent()->getAssociate();
	if (p->isLeaf()) {
		bottom = getHost(p);	//->getAssociate();
	} else {
		std::set<Node*> childHosts;
		for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
			childHosts.insert(getHost(c));	//->getAssociate());
		}
		bottom = H->LCA(childHosts);
	}
	DEBUG(cout << "bottom location = " << bottom->getLabel() << endl);
	DEBUG(cout << "top location = " << top->getLabel() << endl);
	for (Node* n = bottom; n != top; n = n->getParent()) {
		avail.insert(n);
	}
	avail.insert(top);
	return avail;
}

// Recall:
// typedef std::pair<Node*, short> associationtype;	// host x association type
// typedef std::map<Node*, associationtype> nodemaptype;	// parasite x (host x association type)

EventCount CophyMap::countEvents() {
	EventCount E;
	bool _debugging(true);
	for (auto m : phi.getData()) {
//		E.codivs += (m.second.second == 0);	// m is of form (p, (h, type))
//		E.dups += (m.second.second == 1);
		if (m.first->isLeaf()) {
			DEBUG(cout << m.first->getLabel() << " is a leaf: do not count duplication or codivergence events!" << endl);
		} else {
			DEBUG(cout << "counting events for " << m.first->getLabel() << ':' << m.second.first->getLabel() << ": this is a ");
			if (getType(m.first) == 0) {
				++E.codivs;
				DEBUG(cout << "CODIVERGENCE" << endl);
			} else {
				++E.dups;
				DEBUG(cout << "DUPLICATION" << endl);
			}
		}
		if (m.first->hasParent()) {
			DEBUG(cout << "counting losses : number of nodes between "
					<< getHost(m.first)->getLabel() << " and "
					<< getHost(m.first->getParent())->getLabel()
					<< " is " << H->getDistUp(getHost(m.first), getHost(m.first->getParent()))
					<< endl);
			E.losses += H->getDistUp(getHost(m.first), getHost(m.first->getParent()));
			E.losses -= (getType(m.first->getParent()) == 0) ? 1 : 0;
				// if the parent is mapped to a codivergence event then don't count that node as a loss!
		}
	}
	return E;
}

void CophyMap::doPageReconciliation() {
	bool _debugging(false);
//	phi.clear();
	// copy the associations from the CophyMap into the associated parasite vertices
//	for (auto x : phi.getData()) {
//		x.first->setAssociate(x.second.first);
//	}
	DEBUG(cout << "doPageReconciliation" << endl);
	DEBUG(cout << "Root of associate tree = " << P->getRoot()->getLabel() << endl);
	for (Node * n = P->getRoot(); n != nullptr; n = n->next()) {
		if (n->isLeaf()) {
//			DEBUG(cout << n->getLabel() << endl);
//			phi[n].first = n->getAssociate();
			if (getHost(n) != nullptr) { //->getAssociate() != nullptr) {
				DEBUG(cout << n->getLabel() << " is mapped to " << getHost(n)->getLabel() << endl);
			}
		}
	}
	Node* p = P->getRoot();	// root of the parasite / gene tree
	mapToLCAofChildren(p);
	DEBUG(cout << "reconciliation is done" << endl);
}

void CophyMap::inferEvents() {
	for (auto a : P->getVertices()) {
		inferEvents(a.second);	// seems so clunky to not just have a list of vertices in the Tree but hey.
	}
}

void CophyMap::inferEvents(Node *p) {
	/**
	 * If this is a leaf, then no event
	 * Else, this node has children: if they are both mapped to the same node in S, this must be a duplication;
	 * Else, this must be a loss...
	 *
	 * 1. Fix leaves: these are non-events
	 * 2. For each internal node
	 * 	* if t
	 * 	2.1 if at least one child is mapped to the same host (species) node then this must be a duplication
	 * 	2.2
	 */

	DEBUG(cout << "Counting events for node " << p->getLabel() << endl);
//	short event(0);	// 0 for codivergence
	if (p->getDepth() < 0) {
		p->calcDepth();
	}
	if (p->isLeaf()) {
		phi[p].second = noevent;
		DEBUG(cout << "This is a leaf: no event" << endl);
		return;
	}
	std::set<Node*> children;
	p->putChildren(children);
	if (P->LCA(children) == phi[p].first) {
		phi[p].second = codivergence;
	} else {
		phi[p].second = duplication;
	}
//	std::set<Node*> locations;
//	locations.insert(getHost(p));
//	Node* child = p->getFirstChild();
//	uint numLocations(1);
//	while (child != nullptr) {
//		// find child location: which child of the host is the most recent lineage on which the child p is
//		locations.insert(getHost(child));
//		++numLocations;
//		child = child->getSibling();
//	}
//	if (locations.size() < numLocations) {
//		// then there are at least two nodes mapped to the same place so this is a duplication
//		DEBUG(cout << "This is a leaf: no event" << endl);
//		phi[p].second = duplication;
//		return;
//	}
//	phi[p].second = event;
	// Now need to check to see if this node is mapped *above* the highest mapped child: this is also a duplication.
	// XXX did I do the above?
}


bool CophyMap::isValid() {
	bool _valid(true);
	// check all nodes in P are mapped to nodes in H
	Node* image, *childimage;
	for (Node* n = P->getRoot(); n != nullptr; n = n->next()) {
		image = phi.at(n).first;
		if (image->getTree() != H) {
			cerr << "node n=" << n->getLabel() << " is not mapped to the correct tree " << H->getLabel() << endl;
			return false;
		}
		// check no child node is mapped to a node earlier (closer to the root than) the image of the parent (no time travel!)
		for (Node* c = n->getFirstChild(); c != nullptr; c = c->getSibling()) {
			childimage = phi.at(c).first;
			if (H->isAncestralTo(childimage, image)) {
				cerr << "child node c=" << c->getLabel() << " is mapped ancestrally to its parent n=" << n->getLabel() << endl;
				return false;
			}
		}
	}
	return _valid;
}

// recursive function to map node p to the LCA of the images of its children
Node* CophyMap::mapToLCAofChildren(Node* p) {
	bool _debugging(false);
	assert(p);	// just check that this pointer isn't null!
	DEBUG(cout << "mapToLCAofChildren(" << p->getLabel() << ")" << endl);
	if (p->isLeaf()) {
		phi[p].second = 0;	// 0 for codivergence normally
		DEBUG(cout << "This is a leaf: setting image to known associate and returning." << endl);
		DEBUG(cout << p->getLabel() << " maps to " << phi[p].first->getLabel() << endl);
		return phi[p].first;	// return the host/species node
	}
	// first, find to which node will p be mapped:
	assert(p->getFirstChild() != nullptr);
	set<Node*> S;	// this will be where we collect the images of the children
	for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
		S.insert(mapToLCAofChildren(c));
	}
//	Node *lca = S.begin();// = phi[p->getFirstChild()].first;
//	lca = phi.getImage(p->getFirstChild());
	set<Node*>::iterator iter = S.begin();
	Node *lca = *iter;
	++iter;
	while (iter != S.end()) {
		lca = H->LCA(lca, *iter);
		++iter;
	}
	phi[p].first = lca;
//	p->setAssociate(lca);
//	DEBUG(cout << "\t\t" << phi[p->getFirstChild()].first << endl);
//	DEBUG(if (lca == nullptr) {
//		cout << "\t" << p->getFirstChild()->getLabel() << " is not mapped." << endl;
//		cout << "Checking children..." << endl;
//	}
//	);
//	DEBUG(cout << "Now mapping children of " << p->getLabel() << ":" << endl);
//	for (Node * c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
//		DEBUG(cout << "\tchild = " << c->getLabel() << endl);
////		lca = H->LCA(lca, mapToLCAofChildren(c));
//
//		DEBUG(cout << "\taddress of next child: " << c->getSibling() << endl);
//	}
//	DEBUG(cout << "Have mapped children: lca of images is " << lca->getLabel() << endl);
//	phi[p].first = lca;
//	p->setAssociate(lca);
	DEBUG(cout << "Associate " << p->getLabel() << " is mapped to " << lca->getLabel() << endl);
//	DEBUG(
//		for (Node * c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
//			cout << "\t" << c->getLabel() << " -> " << phi[c].second << endl;
//		}
//		cout << "\t" << p->getLabel() << " -> " << lca->getLabel() << endl;
//	);
	// now determine the event type:
	// if all the children of phi[p] are occupied by the children of p then it's a codivergence; else it's duplication.
	std::set<Node*> occupied;
	uint numChildren(0);
	for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
		++numChildren;
		if (phi[c].first == lca) {
			break;
		}
		Node *h = phi[c].first;
		while (h->getParent() != lca) {
			h = h->getParent();
			if (h == nullptr) {
				throw new app_exception("Life is pain, Princess.");
			}
		}
		occupied.insert(h);	// this should be a child node from the parent node p
	}
	phi[p].second = (numChildren == occupied.size()) ? 0 : 1;	// should generalise to multifurcations
	DEBUG(cout << p->getLabel() << ':' << phi[p].first->getLabel() << " is an event of type " << phi[p].second << endl);
	return lca;
}

void CophyMap::moveToHost(Node* p, Node *h) {
	try {
		// have to change the host of p stored in p, the nodemap AND the inversenodemap
		// first remove p from the inverseNodeMap of h:
		DEBUG(cout << "Moving " << p->getLabel() << " from " << phi[p].first->getLabel() << " to " << h->getLabel() << endl);
		Node* currentHost = phi[p].first;
		phi[p].first = h;
		invPhi[currentHost].erase(p);
		invPhi[h].insert(p);
	} catch (const exception& e) {
		cout << e.what();
	}
}

CophyMap& CophyMap::operator=(const CophyMap &other) {
	return *this;
}

void CophyMap::setPhi(string pstr, string hstr) {
	auto VP = P->getVertices();
	auto VH = H->getVertices();
	setPhi(VP[pstr], VH[hstr]);
//	VP[pstr]->setAssociate(VH[hstr]);
}
void CophyMap::setPhi(Node* p, Node* h) {
	phi[p].first = h;
	phi[p].second = -1;	// meaning it has to be recalculated
//	p->setAssociate(h);
	invPhi[h].insert(p);
}

void CophyMap::storeHostInfo() {
	map<Node*, string>& info = P->getInfo();
	info.clear();
	for (auto d : phi.getData()) {
		info[d.first] = d.second.first->getLabel();
		cout << d.first->getLabel() << ":" << d.second.first->getLabel() << endl;
	}
}

ostream& operator<<(ostream& os, CophyMap& T) {
	Tree* H = T.getHostTree();
	Tree* P = T.getParasiteTree();
	os << "Host / independent tree " << H->getLabel() << ":" << endl << *H;
	os << "Parasite / dependent tree " << P->getLabel() << ":" << endl << *P;
	return os;
}

} /* namespace segdup */
