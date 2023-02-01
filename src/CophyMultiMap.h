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

const double defDuplicationCost(10.0);
const double defLossCost(1.0);

class CophyMultiMap {
private:
	std::map<Tree*, CophyMap*> maps; // Parasite tree -> Map
	std::map<std::string, EventCount> mmEventCounts;
	double duplicationCost;
	double lossCost;
public:
	inversenodemap invMap;
	CophyMultiMap(double dup = defDuplicationCost, double loss = defLossCost) : duplicationCost(dup), lossCost(loss) {}
	virtual ~CophyMultiMap() {}

	void addCophyMap(CophyMap* M) { maps[M->getParasiteTree()] = M; }
	int calcDuplicationHeight(Node *h);
	void calcInverseMap();
	void clear();

	EventCount countEvents();
	void doPageReconciliation();

//	CophyMap*& operator[](Tree* P) { return maps[P]; }

//	void toShortDescription(std::string& str);

	EventCount getEventCount(const std::string& mapDescription);
	std::map<Tree*, CophyMap*>& getMaps() { return maps; }

	inline void setDuplicationCost(double d) { duplicationCost = d; }
	inline void storeEventCount(const std::string& label, const EventCount& ec) { mmEventCounts[label] = ec; }
	inline void setLossCost(double l) { lossCost = l; }
	void toCompactString(std::string& str);

};

std::ostream& operator<<(std::ostream& os, CophyMultiMap& CMM);

} /* namespace segdup */

#endif /* SRC_COPHYMULTIMAP_H_ */
