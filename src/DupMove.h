/*
 * DupMove.h
 *
 *  Created on: 6 Dec 2023
 *      Author: mac
 */

#ifndef SRC_DUPMOVE_H_
#define SRC_DUPMOVE_H_

#include "CophyMultiMap.h"

namespace segdup {

/**
 *
 */
class DupMove {
public:
	DupMove();
	virtual ~DupMove();

	virtual void apply(CophyMultiMap& CMM, double T);
};

} /* namespace segdup */

#endif /* SRC_DUPMOVE_H_ */
