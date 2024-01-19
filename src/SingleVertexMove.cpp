/*
 * SingleVertexMove.cpp
 *
 *  Created on: 18 Jan 2024
 *      Author: mac + ybc
 */

#include <iostream>
#include <cmath>

#include "../utility/debugging.h"
#include "../utility/myrandom.h"

#include "Contender.h"
#include "Node.h"
#include "SegDup.h"
#include "SingleVertexMove.h"

using namespace std;

namespace segdup {

double binom(int n, int k) {
    return 1/((n+1)*std::beta(n-k+1,k+1));
}

SingleVertexMove::SingleVertexMove() {
	// TODO Auto-generated constructor stub

}

SingleVertexMove::~SingleVertexMove() {
	// TODO Auto-generated destructor stub
}

void calculateFromVertices(CophyMultiMap& CMM, inversenodemap& invMap, Tree* s, bool up, vector<Node*>& fromVertices) {
	for (auto v : s->getVertices()) {
		if (up && !v.second->hasParent() || !up && v.second->isLeaf())
			continue;

		int dupsAtV = 0;
		for (auto n : invMap[v.second]) 
			if (CMM.getMap(n)->getEvent(n) == duplication)
				dupsAtV++;

		if (dupsAtV >= 2)
			fromVertices.push_back(v.second);
	}
}

void calculateMovableNodes(CophyMultiMap& CMM, inversenodemap& invMap, Node* fromVertex, Node* toVertex, bool up, vector<Node*>& movableNodes) {
	bool _debugging(false);
	DEBUG(
		cout << "invMap[fromVertex]: ";
		for (auto v : invMap[fromVertex]) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)

	for (auto n : invMap[fromVertex]) {
		if (CMM.getMap(n)->getEvent(n) != duplication) 
			continue;

		if (up) {

			if (!n->hasParent()) {
				movableNodes.push_back(n);
				continue;
			}

			if (CMM.getMap(n)->getHost(n->getParent()) != fromVertex) {
				if (CMM.getMap(n)->getHost(n->getParent()) == toVertex && CMM.getMap(n)->getEvent(n->getParent()) != duplication)
					continue;

				movableNodes.push_back(n);
			}
		}
		else {
			bool movable = true;

			for (Node* c = n->getFirstChild(); c != nullptr; c = c->getSibling()) {
				if (!toVertex->isAncestralTo(CMM.getMap(c)->getHost(c)) && toVertex != CMM.getMap(c)->getHost(c))
					movable = false;
			}

			if (movable)
				movableNodes.push_back(n);
		}
	}
}

void SingleVertexMove::apply(CophyMultiMap& CMM, double T) {
	bool _debugging(true);

	DEBUG(cout << "Attempting single vertex move" << endl;)
	DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)

	//choose up or down
	bool up = (iran(2) == 0 ? true : false);
	DEBUG(cout << "Moving " << (up ? "up" : "down") << endl;)

	//find all species tree vertices with >= 2 duplications mapped to them and are non-root/leaf
	Tree* s = CMM.getHostTree();
	inversenodemap& invMap = CMM.getInverseMap();

	vector<Node*> fromVertices;
	calculateFromVertices(CMM, invMap, s, up, fromVertices);
	DEBUG(
		cout << "fromVertices: ";
		for (auto v : fromVertices) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)

	if (fromVertices.size() == 0) {
		DEBUG(cout << "Nothing to move" << endl;)
		return;
	}

	//choose vertex
	Node* fromVertex = fromVertices[iran(fromVertices.size())];
	DEBUG(cout << "Moving from vertex " << fromVertex->getLabel() << endl;)

	//if moving down, choose child
	Node* toVertex;
	if (up)
		toVertex = fromVertex->getParent();
	else {
		toVertex = fromVertex->getFirstChild();
		if (iran(2) == 1)
			toVertex = toVertex->getSibling();
	}
	DEBUG(cout << "Moving to vertex " << toVertex->getLabel() << endl;)
	
	//find all movable nodes for the chosen vertex
	vector<Node*> movableNodes;
	calculateMovableNodes(CMM, invMap, fromVertex, toVertex, up, movableNodes);
	DEBUG(
		cout << "movableNodes: ";
		for (auto v : movableNodes) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)

	//reject if not enough movable nodes
	if (movableNodes.size() < 2) {
		DEBUG(cout << "Not enough movable nodes at chosen vertex" << endl;)
		return;
	}

	//choose number of nodes to move
	int movableNodesSize = movableNodes.size();
	int choose = iran((movableNodesSize*(movableNodesSize-1))/2);
	int numberToMove = 2;
	while (choose >= numberToMove) {
		choose -= numberToMove;
		numberToMove++;
	}
	DEBUG(cout << "Moving " << numberToMove << " nodes" << endl;)

	//choose nodes to move
	//changes movableNodes!
	vector<Node*> nodesToMove;
	for (int i = 0; i < numberToMove; i++) {
		int j = iran(movableNodes.size());
		nodesToMove.push_back(movableNodes[j]);
		movableNodes.erase(movableNodes.begin()+j);
	}
	DEBUG(
		cout << "nodesToMove: ";
		for (auto v : nodesToMove) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)
	
	//move nodes + calculate cost
	EventCount ec;
	ec.losses = (up ? 1 : -1) * numberToMove;

	int oldJointDupHeightHere(CMM.calcCombinedDuplicationHeight(fromVertex));
	int oldJointDupHeightThere(CMM.calcCombinedDuplicationHeight(toVertex));

	for (auto n : nodesToMove) {
		CMM.getMap(n)->moveToHost(n, toVertex, duplication);
		CMM.movePToHost(n, fromVertex, toVertex);

		if (!n->hasParent())
			ec.losses += (up ? 1 : -1);
	}
	DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)

	int nuJointDupHeightHere(CMM.calcCombinedDuplicationHeight(fromVertex));
	int nuJointDupHeightThere(CMM.calcCombinedDuplicationHeight(toVertex));

	ec.dups = nuJointDupHeightThere - oldJointDupHeightThere;
	ec.dups += nuJointDupHeightHere - oldJointDupHeightHere;	//here is never there

	//calculate stuff for reverse transition probability
	vector<Node*> reverseFromVertices;
	vector<Node*> reverseMovableNodes;
	calculateFromVertices(CMM, invMap, s, !up, reverseFromVertices);
	calculateMovableNodes(CMM, invMap, toVertex, fromVertex, !up, reverseMovableNodes);
	DEBUG(
		cout << "reverseFromVertices: ";
		for (auto v : reverseFromVertices) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)
	DEBUG(
		cout << "reverseMovableNodes: ";
		for (auto v : reverseMovableNodes) {
			cout << v->getLabel() << " ";
		}
		cout << endl;
	)

	//accept/reject
	double MHProb = (up ? 0.5 : 2.0);
	MHProb *= fromVertices.size();
	MHProb /= reverseFromVertices.size();
	MHProb *= movableNodesSize*(movableNodesSize-1);
	MHProb /= reverseMovableNodes.size()*(reverseMovableNodes.size()-1);
	MHProb *= binom(movableNodesSize, numberToMove);
	MHProb /= binom(reverseMovableNodes.size(), numberToMove);
	MHProb *= exp(-CSD(ec)/T);

	DEBUG(cout << "M-H probability = " << (up ? 0.5 : 2.0) << " * " << fromVertices.size() << " / " << reverseFromVertices.size() << " * " << movableNodesSize*(movableNodesSize-1) << " / " << reverseMovableNodes.size()*(reverseMovableNodes.size()-1) << " * " << binom(movableNodesSize, numberToMove) << " / " << binom(reverseMovableNodes.size(), numberToMove) << " * " << exp(-CSD(ec)/T) << endl;)
	DEBUG(cout << "M-H probability = " << MHProb << endl;)
	
	if (dran(1) < MHProb) {
		//accept
		ec += CMM.countEvents();
		CMM.setCurrentEventCount(ec);
		DEBUG(cout << "SingleVertex move succeeded!" << endl;)
	}
	else {
		//reject, move back
		for (auto n : nodesToMove) {
			CMM.getMap(n)->moveToHost(n, fromVertex, duplication);
			CMM.movePToHost(n, toVertex, fromVertex);
		}
	}

	DEBUG(cout << endl;)
}

}; // end of namespace
