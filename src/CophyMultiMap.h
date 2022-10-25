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

class CophyMultiMap {
private:
	std::vector<CophyMap*> maps;
	inversenodemap invMap;
public:
	CophyMultiMap();
	virtual ~CophyMultiMap();

	void addCophyMap(CophyMap* M) { maps.push_back(M); }
	int calcDuplicationHeight(Node *h);
	void calcInverseMap();

	EventCount countEvents();
	void doPageReconciliation();

};

} /* namespace segdup */

#endif /* SRC_COPHYMULTIMAP_H_ */
