/*
 * CophyMultiMap.cpp
 *
 *  Created on: 24 Oct 2022
 *      Author: mac
 */

#include "CophyMultiMap.h"

using namespace std;

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
		int height(1);
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
			invMap[iter->second.first].insert(iter->first);
			// iter is the association, that is, mapping Node* p to pair (Node* h, assoc type).
		}
	}
}

EventCount CophyMultiMap::countEvents() {
	EventCount EC;
	for (CophyMap *Phi : maps) {
		EC += Phi->countEvents();
	}
	return EC;
}

void CophyMultiMap::doPageReconciliation() {
	for (CophyMap *Phi : maps) {
		Phi->doPageReconciliation();
	}
}


} /* namespace segdup */
