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
extern bool _cacheEventCounts;

namespace segdup {

int CophyMultiMap::calcCombinedDuplicationHeight(Node *h) {
	/**
	 * Find each node in each parasite/gene tree mapped to node h.
	 * For each, calculate the height of each gene subtree that is mapped to h
	 * The duplication height is the maximum of these.
	 *
	 * XXX should not be re-counting for any node that hasn't moved
	 *
	 */
	bool _debugging(false);
	DEBUG(cout << "Calculating duplication height for host node " << h->getLabel() << endl);
//	if (invMap.size() == 0) {
//		DEBUG(cout << "calculating inverse map since it is empty" << endl);
		calcInverseMap();	// TODO find out why this isn't being updated properly by moving p-nodes around.
		// XXX I've turned OFF the test to see if the invMap.size() is zero and it SEEMS to be working.. (MAC 2022/12/05)
//	}
	int dupHeight(0);
	for (Node* p : invMap[h]) {
//		DEBUG(cout << p->getLabel() << " is on host " << h->getLabel() << endl);
		CophyMap* M = maps[p->getTree()];
//		int height = (M->getEvent(p) == duplication) ? 1 : 0;	// this is probably dodgy because it treats "noevent" and "loss" as "codivergence"
//		DEBUG(cout << "\tinitial height for this lineage is " << height << endl);
//		while (p->getParent() != nullptr) {
//			p = p->getParent();	// XXX This is a good thing to optimise for speed
//			if (invMap[h].count(p) > 0) {
//				++height;
//			} else {
//				DEBUG(cout << "\theight for this p = " << height << endl);
//				DEBUG(cout << "\tdupHeight = " << dupHeight << endl);
//				break;
//			}
//		}
		DEBUG(cout << "Considering node " << p->getLabel() << " on tree " << p->getTree()->getLabel() << endl);
		dupHeight = max<int>(dupHeight, M->getDuplicationHeight(p));
	}
	DEBUG(cout << "\tduplication height is " << dupHeight << endl);
	DEBUG(cout << *this << endl);
	return dupHeight;
}

void CophyMultiMap::calcInverseMap() {
	invMap.clear();
	for (auto mpr : maps) {
		CophyMap* M = mpr.second;
		nodemaptype& phi = M->getPhi().getData();
		for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
			invMap[iter->second].insert(iter->first);
			// iter is the association, that is, mapping Node* p to pair (Node* h, assoc type).
		}
	}
}

void CophyMultiMap::clear() {
	maps.clear();
	invMap.clear();
	duplicationCost = defDuplicationCost;
	lossCost = defLossCost;
}

EventCount CophyMultiMap::countEvents() {
	/**
	 *
	 */
	bool _debugging(false);
	EventCount E;
	Tree *H;
	Node *p, *h;
	string description;
	if	(_cacheEventCounts) {
		toCompactString(description);
		if (mmEventCounts.find(description) != mmEventCounts.end()) {
//			cerr << "size of score table = " << mmEventCounts.size() << endl;
			return mmEventCounts[description];
		}
	}
	int nLosses;
	DEBUG(cout << "CophyMultiMap:countEvents()\n");
	for (auto mpr : maps) {
		CophyMap* M = mpr.second;
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
				if (M->getEvent(p) == codivergence) {
					++E.codivs;
					DEBUG(cout << "counting events for " << p->getLabel() << ':' << h->getLabel() << ": this is a ");
					DEBUG(cout << "CODIVERGENCE" << endl);
				}
			}
			if (m.first->hasParent()) {
				nLosses = H->getDistUp(h, phi[p->getParent()]);
				nLosses -= (M->getEvent(p->getParent()) == codivergence) ? 1 : 0;
				nLosses = max(0, nLosses);
				E.losses += nLosses;
				DEBUG(cout << "counting losses leading to " << p->getLabel() << ":"
						<< h->getLabel() << " -- number of nodes between "
						<< h->getParent()->getLabel() << " and "
						<< h->getLabel()
						<< " is " << H->getDistUp(h, phi[p->getParent()])
						<< endl);
				DEBUG(cout << "\tadding " << nLosses << " losses." << endl);
					// if the parent is mapped to a codivergence event then don't count that node as a loss!
			} else {
//				E.losses += calcDuplicationHeight(m.second);
			}
		}
	}
	// now count segmental duplications (which may be shared by multiple tree-maps):
	map<string, Node*>& V = H->getVertices();
	set<Node*> countedHosts;
	for (auto iter : V) {
		if (countedHosts.count(iter.second) > 0) {
			continue;
		}
		E.dups += calcCombinedDuplicationHeight(iter.second);
		countedHosts.insert(iter.second);
		// XXX assumes all segmental duplications are permitted -- everything is adjacent!
	}
	if (_cacheEventCounts) {
		storeEventCount(description, E);
	}
	return E;
}

EventCount CophyMultiMap::getEventCount(const std::string& mapDescription) {
	return mmEventCounts.at(mapDescription);
}

void CophyMultiMap::doPageReconciliation() {
	bool _debugging(true);
	for (auto mpr : maps) {
		CophyMap *Phi = mpr.second;
		Phi->doPageReconciliation();
		Phi->inferEvents();
		DEBUG(cout << *Phi << endl);
	}
}

void CophyMultiMap::toCompactString(string & str) {
	str = "";
//	EventCount ec = countEvents();
//	str = to_string(ec.dups) + " dups + " + to_string(ec.losses) + " losses; score=";
//	str += to_string(ec.dups * duplicationCost + ec.losses*lossCost);
//	str += "\t";
	Node *p;
	Node *h;
	for (auto mpr : maps) {
		CophyMap *M = mpr.second;
		NodeMap& phi = M->getPhi();
		for (auto m : phi.getData()) {
			p = m.first;
			if (p->isLeaf()) {
				continue;
			}
			h = m.second;
			str += p->getLabel() + ":[" + eventSymbol[M->getEvent(p)] + "]" + h->getLabel() + ",";
		}
	}
	str.pop_back();
}

ostream& operator<<(ostream& os, CophyMultiMap& CMM) {
	os << "CophyMultiMap of " << endl;
	Tree* H(nullptr);
	auto mitr = CMM.getMaps().begin();
	H = mitr->second->getHostTree();
	for (auto mpr : CMM.getMaps()) {
		os << '\t' << mpr.first->getLabel() << "->" << mpr.second->getHostTree()->getLabel() << endl;
	}
	os << H->getLabel() << endl << *H;
	for (auto mpr : CMM.getMaps()) {
		mpr.first->setShowInfo(true);
		os << mpr.first->getLabel() << endl << *(mpr.first);
		os << "Event count: " << mpr.second->countEvents() << endl;
	}
	os << "{Key: duplication marked by [" << eventSymbol[duplication] << "]; codivergence by ["
			<< eventSymbol[codivergence] << "].}" << endl;
	return os;
}


} /* namespace segdup */
