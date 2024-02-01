/*
 * CophyMultiMap.h
 *
 *  Created on: 24 Oct 2022
 *      Author: mac
 */

#ifndef SRC_COPHYMULTIMAP_H_
#define SRC_COPHYMULTIMAP_H_

#include <vector>
#include "CophyMap.h"

namespace segdup {

const double defDuplicationCost(5.0);
const double defLossCost(1.0);

class CophyMultiMap {
private:
	std::map<Tree*, CophyMap*> maps; // Parasite tree -> Map
	std::map<std::string, EventCount> mmEventCounts;
	std::vector<std::pair<Node*,CophyMap*>> allMoveableNodes;
	inversenodemap invMap;
	EventCount currentEventCount;
	double duplicationCost;
	double lossCost;
	Tree* H;
public:
	CophyMultiMap(double dup = defDuplicationCost, double loss = defLossCost) : duplicationCost(dup), lossCost(loss), H(nullptr) {}
	virtual ~CophyMultiMap() {}

	void addCophyMap(CophyMap* M) { maps[M->getParasiteTree()] = M; H = M->getHostTree(); }
	int calcCombinedDuplicationHeight(Node *h);
	void calcInverseMap();
	inversenodemap& getInverseMap();
	void clear();

	EventCount& countEvents(const std::string& mapDescription);
	EventCount& countEvents();
	void setCurrentEventCount(const EventCount& ec);
	void calcEventCount();	// XXX change to count events
	void doEarlyReconciliation();
	void doPageReconciliation();

//	CophyMap*& operator[](Tree* P) { return maps[P]; }

//	void toShortDescription(std::string& str);

	std::map<Tree*, CophyMap*>& getMaps() { return maps; }
	CophyMap* getMap(Node* n) { return maps[n->getTree()]; }
	Tree* getHostTree() { return H; }

	void movePToHost(Node* p, Node *oldHost, Node *nuHost);

	void putAllMoveableNodes();
	std::vector<std::pair<Node*,CophyMap*>>& getAllMoveableNodes();
	inline void setDuplicationCost(double d) { duplicationCost = d; }
	inline void storeEventCount(const std::string& label, const EventCount& ec) { mmEventCounts[label] = ec; }
	inline void setLossCost(double l) { lossCost = l; }
	void toCompactString(std::string& str);

};

std::ostream& operator<<(std::ostream& os, CophyMultiMap& CMM);

} /* namespace segdup */

#endif /* SRC_COPHYMULTIMAP_H_ */
