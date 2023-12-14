/*
 * SingleNodeMove.h
 *
 *  Created on: 14 Dec 2023
 *      Author: mac
 */

#ifndef SRC_SINGLENODEMOVE_H_
#define SRC_SINGLENODEMOVE_H_

#include "DupMove.h"
#include "CophyMultiMap.h"

namespace segdup {

class SingleNodeMove: public segdup::DupMove {
public:
	SingleNodeMove();
	virtual ~SingleNodeMove();

	void apply(CophyMultiMap& CMM, double T);
};

}
#endif /* SRC_SINGLENODEMOVE_H_ */
