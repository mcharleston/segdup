/*
 * Adjacency.h
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#ifndef ADJACENCY_H_
#define ADJACENCY_H_

#include <map>

#include "Node.h"

namespace segdup {

/**
 * Stores the collection of which gene tree branches are adjacent to which others.
 *
 */

class Adjacency {
private:
	// ?
public:
	Adjacency();
	virtual ~Adjacency();


	// XXX Override this when NOT everything is adjacent. For now, it's fine.
	bool adjacent(Node* x, Node* u, Node* v) { return true; }
};

} /* namespace segdup */

#endif /* ADJACENCY_H_ */
