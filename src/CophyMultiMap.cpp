/*
 * CophyMultiMap.cpp
 *
 *  Created on: 24 Oct 2022
 *      Author: mac
 */

#include "../utility/debugging.h"
#include "CophyMultiMap.h"

using namespace std;

extern bool _debugging;

namespace segdup {

CophyMultiMap::CophyMultiMap() {
	// TODO Auto-generated constructor stub

}

CophyMultiMap::~CophyMultiMap() {
	// TODO Auto-generated destructor stub
}

int CophyMultiMap::calcDuplicationHeight(Node *h) {
	/**
	 * Find each node in each parasite/gene tree mapped to node h.
	 * For each, calculate the height of each gene subtree that is mapped to h
	 * The duplication height is the maximum of these.
	 */
	if (invMap.size() == 0) {
		calcInverseMap();
	}
	int dupHeight(0);
	for (Node* g : invMap[h]) {
		int height = (g->event == duplication) ? 1 : 0;	// this is probably dodgy because it treats "noevent" and "loss" as "codivergence"
		while (g->getParent() != nullptr) {
			g = g->getParent();
			if (invMap[h].count(g) > 0) {
				++height;
			} else {
				dupHeight = max<int>(dupHeight, height);
				break;
			}
		}
	}
	return dupHeight;
}

void CophyMultiMap::calcInverseMap() {
	invMap.clear();
	for (CophyMap* M : maps) {
		nodemaptype& phi = M->getPhi().getData();
		for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
			invMap[iter->second].insert(iter->first);
			// iter is the association, that is, mapping Node* p to pair (Node* h, assoc type).
		}
	}
}

EventCount CophyMultiMap::countEvents() {
	/**
	 *
	 */
	bool _debugging(true);
	EventCount E;
	Tree *H;
	Node *p, *h;
	int nLosses;
	for (CophyMap* M : maps) {
		NodeMap& phi = M->getPhi();
		// Count duplications and losses for separate maps:
		H = phi.getHostTree();
		for (auto m : phi.getData()) {
			p = m.first;
			h = m.second;
	//		E.codivs += (m.second.second == 0);	// m is of form (p, (h, type))
	//		E.dups += (m.second.second == 1);
			if (m.first->isLeaf()) {
				DEBUG(cout << p->getLabel() << " is a leaf: do not count duplication or codivergence events!" << endl);
			} else {
				if (p->getEvent() == codivergence) {
					++E.codivs;
					DEBUG(cout << "counting events for " << p->getLabel() << ':' << h->getLabel() << ": this is a ");
					DEBUG(cout << "CODIVERGENCE" << endl);
				}
			}
			if (m.first->hasParent()) {
				nLosses = H->getDistUp(h, phi[p->getParent()]);
				nLosses -= (p->getParent()->getEvent() == codivergence) ? 1 : 0;
				nLosses = max(0, nLosses);
				E.losses += nLosses;
				DEBUG(cout << "counting losses leading to " << p->getLabel() << ":"
						<< h->getLabel() << " -- number of nodes between "
						<< h->getLabel() << " and "
						<< h->getParent()->getLabel()
						<< " is " << H->getDistUp(h, phi[p->getParent()])
						<< endl);
				DEBUG(cout << "\tadding " << nLosses << " losses." << endl);
					// if the parent is mapped to a codivergence event then don't count that node as a loss!
			}
		}
	}
	// now count segmental duplications (which may be shared by multiple tree-maps):
	map<string, Node*>& V = H->getVertices();
	for (auto iter : V) {
		E.dups += calcDuplicationHeight(iter.second);
	}
	return E;
}

void CophyMultiMap::doPageReconciliation() {
	for (CophyMap *Phi : maps) {
		Phi->doPageReconciliation();
	}
}


} /* namespace segdup */
