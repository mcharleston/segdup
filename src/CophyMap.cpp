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
	invPhi = other.invPhi;
	event = other.event;
}

std::set<pair<Node*, eventType>> CophyMap::calcAvailableNewHosts(Node* p) {
	// can go up to same vertex to which parent is mapped, and down to the LCA of there the children are mapped.
	//
	bool _debugging(false);
	DEBUG(cout << "Calculating available hosts for node " << p->getLabel() << ":" << endl);
	std::set<pair<Node*, eventType>> avail;
	Node* bottom;
	Node* top;
	bool _addCodivergence(false);
	avail.insert(make_pair(getHost(p), event[p]));
	pair<Node*, eventType> nuCodiv;
	if (p->hasParent()) {
		DEBUG(cout << "p has parent with node label " << p->parent->getLabel() << endl);
		top = getHost(p->getParent()); //->getParent()->getAssociate();
	} else {
		top = getHost(p)->getTree()->getRoot();
	}
	if (p->isLeaf()) {
		bottom = getHost(p);	//->getAssociate();
	} else {
		std::set<Node*> childHosts;
		for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
			childHosts.insert(getHost(c));	//->getAssociate());
		}
		bottom = H->LCA(childHosts);
		if (childHosts.count(bottom) == 0) {
			_addCodivergence = true;
			nuCodiv = make_pair(bottom, codivergence);
		}
	}
	DEBUG(cout << "bottom location = " << bottom->getLabel() << endl);
	DEBUG(cout << "top location = " << top->getLabel() << endl);
	for (Node* n = bottom; n != top; n = n->getParent()) {
		avail.insert(make_pair(n, duplication));
	}
	if (p->hasParent()) {
		if (event[p->getParent()] == duplication) {
			avail.insert(make_pair(top, duplication));
			DEBUG(cout << "Adding top potential location: duplication at " << top->getLabel() << endl);
		}
	} else {
//		if (getHost(p->getParent()) == top && event[p->getParent()] == duplication) {
			avail.insert(make_pair(top, duplication));
//		}
	}
	if (_addCodivergence) {
		avail.insert(nuCodiv);
	}
	return avail;
}

void CophyMap::checkValidHostOrdering() {
	std::map<std::string, Node*>& V = P->getVertices();
	for (auto vpr : V) {
		Node *p = vpr.second;
		if (p->hasParent()) {
			if (getHost(p)->isAncestralTo(getHost(p->getParent()))) {
				throw new app_exception("Invalid map! " );//+ p->getLabel() + " has host " + getHost(p)
//						+ " but its parent node " + p->getParent()->getLabel() + " has host " + getHost(p->getParent()));
			}
		}
	}
}


// Recall:
//typedef Node* associationtype;	// host x association type
//typedef std::map<Node*, associationtype> nodemaptype;	// parasite x (host x association type)
//typedef std::map<Node*, std::set<associationtype> > wideassociationmap;
//typedef std::map<Node*, std::set<Node*>> inversenodemap;

EventCount CophyMap::countEvents() {
	EventCount E;
	bool _debugging(false);
	for (auto m : phi.getData()) {
//		E.codivs += (m.second.second == 0);	// m is of form (p, (h, type))
//		E.dups += (m.second.second == 1);
		if (m.first->isLeaf()) {
			DEBUG(cout << m.first->getLabel() << " is a leaf: do not count duplication or codivergence events!" << endl);
		} else {
			DEBUG(cout << "counting events for " << m.first->getLabel() << ':' << m.second->getLabel() << ": this is a ");
			if (getEvent(m.first) == codivergence) {
				++E.codivs;
				DEBUG(cout << "CODIVERGENCE" << endl);
			} else {
				++E.dups;
				DEBUG(cout << "DUPLICATION" << endl);
			}
		}
		if (m.first->hasParent()) {
			DEBUG(cout << "counting losses : number of nodes between "
					<< getHost(m.first->getParent())->getLabel() << " and "
					<< getHost(m.first)->getLabel()
					<< " is " << H->getDistUp(getHost(m.first), getHost(m.first->getParent()))
					<< endl);
			E.losses += H->getDistUp(getHost(m.first), getHost(m.first->getParent()));
			E.losses -= (getEvent(m.first->getParent()) == 0) ? 1 : 0;
				// if the parent is mapped to a codivergence event then don't count that node as a loss!
		}
	}
	return E;
}

std::string CophyMap::describeEvent(eventType e) {
	switch (e) {
		case codivergence: return "codivergence";
		case duplication: return "duplication";
		case loss: return "loss";
		case noevent: return "noevent";
		default: break;
	}
	return "";
}

void CophyMap::doPageReconciliation() {
	bool _debugging(false);
	DEBUG(cout << "doPageReconciliation" << endl);
	DEBUG(cout << "Root of associate tree = " << P->getRoot()->getLabel() << endl);
	for (Node * n = P->getRoot(); n != nullptr; n = n->next()) {
		if (n->isLeaf()) {
			if (getHost(n) != nullptr) { //->getAssociate() != nullptr) {
				DEBUG(cout << n->getLabel() << " is mapped to " << getHost(n)->getLabel() << endl);
			}
		}
	}
	Node* p = P->getRoot();	// root of the parasite / gene tree
	mapToLCAofChildren(p);
	storeAssociationInfo();
	DEBUG(cout << "reconciliation is done" << endl);
}

void CophyMap::inferEvents() {
	for (auto a : P->getVertices()) {
		inferEvents(a.second);	// seems so clunky to not just have a list of vertices in the Tree but hey.
	}
	storeAssociationInfo();
}

void CophyMap::inferEvents(Node *p) {
	bool _debugging(false);
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

	DEBUG(cout << "Counting events for node " << p->getLabel() << ":" << phi[p]->getLabel() << endl);
	if (p->getDepth() < 0) {
		p->calcDepth();
	}
	if (p->isLeaf()) {
		event[p] = noevent;
		DEBUG(cout << "This is a leaf: no event" << endl);
		return;
	}
	std::set<Node*> pchildren, himages;
	p->putChildren(pchildren);
	for (Node* q : pchildren) {
		himages.insert(phi[q]);
	}
	DEBUG(cout << "\tHosts of children: { "; for (Node* hc : himages) cout << hc->getLabel() << " "; cout << "}" << endl );
	DEBUG(cout << "\tLCA(<children>) = " << H->LCA(himages)->getLabel() << "; phi[p] = " << phi[p]->getLabel() << endl);
	if ((H->LCA(himages) == phi[p]) && (himages.count(phi[p]) == 0)) {
		event[p] = codivergence;
	} else {
		event[p] = duplication;
	}
	DEBUG(cout << "\tEvent for node " << p->getLabel() << " is " << describeEvent(p) << endl);
	// Now need to check to see if this node is mapped *above* the highest mapped child: this is also a duplication.
	// XXX did I do the above?
}


bool CophyMap::isValid() {
	bool _valid(true);
	// check all nodes in P are mapped to nodes in H
	Node* image, *childimage;
	for (Node* n = P->getRoot(); n != nullptr; n = n->next()) {
		image = phi.at(n);
		if (image->getTree() != H) {
			cerr << "node n=" << n->getLabel() << " is not mapped to the correct tree " << H->getLabel() << endl;
			return false;
		}
		// check no child node is mapped to a node earlier (closer to the root than) the image of the parent (no time travel!)
		for (Node* c = n->getFirstChild(); c != nullptr; c = c->getSibling()) {
			childimage = phi.at(c);
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
//		p->event = codivergence;	// 0 for codivergence normally
		DEBUG(cout << "This is a leaf: setting image to known associate and returning." << endl);
		DEBUG(cout << p->getLabel() << " maps to " << phi[p]->getLabel() << endl);
		return phi[p];	// return the host/species node
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
	phi[p] = lca;
	DEBUG(cout << "Associate " << p->getLabel() << " is mapped to " << lca->getLabel() << endl);
	// now determine the event type:
	// if all the children of phi[p] are occupied by the children of p then it's a codivergence; else it's duplication.
	std::set<Node*> occupied;
	uint numChildren(0);
	for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
		++numChildren;
		if (phi[c] == lca) {
			break;
		}
		Node *h = phi[c];
		while (h->getParent() != lca) {
			h = h->getParent();
			if (h == nullptr) {
				throw new app_exception("Life is pain, Princess.");
			}
		}
		occupied.insert(h);	// this should be a child node from the parent node p
	}
//	p->event = (numChildren == occupied.size()) ? codivergence : duplication;	// should generalise to multifurcations
	DEBUG(cout << p->getLabel() << ':' << phi[p]->getLabel() << " is an event of type " << getEvent(p) << endl);
	return lca;
}

void CophyMap::moveToHost(Node* p, Node *h) {
	try {
		// have to change the host of p stored in phi, the nodemap AND the inversenodemap
		// first remove p from the inverseNodeMap of h:
		Node* currentHost = phi[p];
//		for (unsigned int bug(0); bug < 10000; ++bug) {
//			phi[p] = h;
//			invPhi[currentHost].erase(p);
//			invPhi[h].insert(p);
//			phi[p] = currentHost;
//			invPhi[h].erase(p);
//			invPhi[currentHost].insert(p);
//		}
		phi[p] = h; // XXX is phi, being a map of Node* to Node*, just getting bigger and bigger?
		invPhi[currentHost].erase(p);
		invPhi[h].insert(p);
//		cerr << invPhi.size() << ", " << phi.getData().size() << endl;	// apparently not.
//		for (auto a : invPhi) {
//			cerr << a.second.size() << ',';
//		}
//		cerr << endl;
//		P->getInfo().at(p) = h->getLabel();
	} catch (const exception& e) {
		cout << e.what();
	}
}
void CophyMap::moveToHost(Node* p, Node *h, eventType e) {
	// XXX MEMORY LEAK IN HERE?
//	DEBUG(cout << "Moving " << p->getLabel() << " from " << eventSymbol[event[p]] << phi[p]->getLabel()
//			<< " to " << eventSymbol[e] << h->getLabel() << endl);
	moveToHost(p, h);
	event[p] = e;
	storeAssociationInfo(p, h, e);
}

CophyMap& CophyMap::operator=(const CophyMap &other) {
	H = other.H;	// the species or host tree
	P = other.P;	// the gene or parasite tree
	phi = other.phi;	// nodes in P mapped to nodes in H -- starting with the leaves
	invPhi = other.invPhi;
	event = other.event;
	return *this;
}

void CophyMap::setPhi(string pstr, string hstr) {
	auto VP = P->getVertices();
	auto VH = H->getVertices();
	setPhi(VP[pstr], VH[hstr]);
}
void CophyMap::setPhi(Node* p, Node* h) {
	phi[p] = h;
	event[p] = noevent;	// meaning it has to be recalculated
	invPhi[h].insert(p);
}

void CophyMap::storeAssociationInfo() {
	/**
	 * Put info into the info map for each parasite/gene node:
	 * [<]p:h to show p codiverges on node h
	 * [=]p:h to show p duplicates on node h
	 * [x]p:h to show p has a loss (extinction, sampling failure, or miss-the-boat) on h
	 * p:h to show no event, such as at leaves
	 */
	bool _debugging(false);
	info.clear();
	for (auto d : phi.getData()) {
		storeAssociationInfo(d.first, d.second, event[d.first]);
	}
	P->setInfo(&info);
}

void CophyMap::storeAssociationInfo(Node* p, Node* h, eventType e) {
	info[p] = "";
	if (!p->isLeaf()) {
		string str = "[";
		str += eventSymbol[event[p]];
		str += "]";
		info[p] = str;
	}
	info[p] += p->getLabel() + ":" + h->getLabel();
//	info[d.first] = d.second->getLabel();
//	DEBUG(cout << info[p] << endl);
}

ostream& operator<<(ostream& os, CophyMap& T) {
	Tree* H = T.getHostTree();
	Tree* P = T.getParasiteTree();
//	os << "Host / independent tree " << H->getLabel() << ":" << endl << *H;
	os << "Parasite / dependent tree " << P->getLabel() << ":" << endl << *P;
	return os;
}

} /* namespace segdup */
