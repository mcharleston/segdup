/*
 * Reconciliation.cpp
 *
 *  Created on: 5 Sept 2022
 *      Author: mac
 */

#include <map>

#include "../utility/appexception.h"
#include "Reconciliation.h"
#include "Tree.h"
#include "NodeMap.h"

using namespace std;

namespace segdup {

Reconciliation::~Reconciliation() {
	// TODO Auto-generated destructor stub
}

void Reconciliation::countEvents(eventCounts& C) {
//	for (CophyMap M : AllMappedTrees) {
//		// phi maps a tree to a cophymap
//		// for each mapped node, if it's a codivergence then count it as one;
//		// else it's a duplication; also add losses for how many edges from it to the image of its parent
//		NodeMap& phi = M.getPhi();
//		nodemaptype& NM = phi.getData();
//		for (auto im : NM) {
//			if (im.first->event == codivergence) {
//				C.nCodivs++;
//			} else {
//				C.nDups++;
//			}
//			Node* par = im.second.first->getParent();
//			if (par != nullptr) {
//				C.nLosses += H->getDistUp(im.second.first, par);
//				C.nLosses -= (phi[par].second == duplication);
//			}
//		}
//	}
}

void Reconciliation::move(Tree* T, Node* n, Node* to) {
//	Tree *H = ;
//	if (n->getAssociate()->getTree() != to->getTree()) {
//		throw new app_exception("Cannot move image of node " + n->getLabel() + " to a different tree!");
//	}
//	if (!n->isLeaf()) {
//		// check children are not mapped strictly BEFORE the new destination: if any are, then throw an exception.
//	}
//	if (n->hasParent()) {
//		// check parent is not mapped strictly AFTER the new destination: if it is, throw an exception.
////		auto M = AllMappedTrees.at(T);
//
//	}
}

} /* namespace segdup */
