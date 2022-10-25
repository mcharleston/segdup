/*
 * Reconciliation.h
 *
 *  Created on: 5 Sept 2022
 *      Author: mac
 */

#ifndef SRC_RECONCILIATION_H_
#define SRC_RECONCILIATION_H_

#include <vector>

#include "EventCount.h"
#include "Tree.h"
#include "CophyMap.h"

/************************************************************************************************************
 * Notes
 *
 * SDA := Segmental Duplication with Adjacency
 * c_SDA :=
 */

namespace segdup {

struct eventCounts {
public:
	int nCodivs;
	int nDups;
	int nLosses;
};
// maps Gene/Parasite Node -> (Species/Host Node, event)

class Reconciliation {
private:
	std::vector<CophyMap> AllMappedTrees;	// map of gene trees to cophylogeny maps: Phi[T] : V(T) -> V(H) x {eventtypes}

	Tree* H;	// the species/host tree
	float dupCost;	// ("delta")
	float lossCost; // ("lambda")
public:
	Reconciliation() : H(nullptr), dupCost(1.0), lossCost(1.0) {}
	virtual ~Reconciliation();

	void countEvents(eventCounts& C); // count and return the combined event counts across all trees

	void move(Tree* T, Node* n, Node* to); // move node N in tree T to a new destination in the mapped-to tree

	Reconciliation& operator=(const Reconciliation &other);
};

} /* namespace segdup */

#endif /* SRC_RECONCILIATION_H_ */
