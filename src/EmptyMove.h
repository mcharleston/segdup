/*
 * EmptyMove.h
 *
 *  Created on: 14 Dec 2023
 *      Author: mac
 */

#ifndef SRC_EMPTYMOVE_H_
#define SRC_EMPTYMOVE_H_

#include "SegDup.h"
#include "DupMove.h"
#include "CophyMultiMap.h"

namespace segdup {

class EmptyMove: public segdup::DupMove {
public:
	EmptyMove();
	virtual ~EmptyMove();

	void apply(CophyMultiMap& CMM, double T);
};

}
#endif /* SRC_EMPTYMOVE_H_ */
