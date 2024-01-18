/*
 * SingleVertexMove.h
 *
 *  Created on: 18 Jan 20234
 *      Author: mac + ybc
 */

#ifndef SRC_SINGLEVERTEXMOVE_H_
#define SRC_SINGLEVERTEXMOVE_H_

#include "SegDup.h"
#include "DupMove.h"
#include "CophyMultiMap.h"

namespace segdup {

class SingleVertexMove: public segdup::DupMove {
public:
	SingleVertexMove();
	virtual ~SingleVertexMove();

	void apply(CophyMultiMap& CMM, double T);
};

}
#endif /* SRC_SINGLEVERTEXMOVE_H_ */
