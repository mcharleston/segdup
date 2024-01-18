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

SingleVertexMove::SingleVertexMove() {
	// TODO Auto-generated constructor stub

}

SingleVertexMove::~SingleVertexMove() {
	// TODO Auto-generated destructor stub
}

void SingleVertexMove::apply(CophyMultiMap& CMM, double T) {
	bool _debugging(true);

	DEBUG(cout << "Attempting single vertex move" << endl;)
	DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)

	//choose up or down
	bool up = (iran(2) == 0 ? true : false);

	//find all species tree vertices with >= 2 duplications mapped to them and are non-root/leaf
	Tree* s = CMM.getHostTree();
	inversenodemap& invMap = CMM.getInverseMap();

	vector<Node*> fromVertices;
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

	if (fromVertices.size() == 0) {
		DEBUG(cout << "Nothing to move" << endl;)
		return;
	}

	//choose vertex
	Node* fromVertex = fromVertices[iran(fromVertices.size())];

	//if moving down, choose child
	Node* toVertex;
	if (up)
		toVertex = fromVertex->getParent();
	else {
		toVertex = fromVertex->getFirstChild();
		if (iran(2) == 1)
			toVertex = toVertex->getSibling();
	}
	
	//find all movable nodes for the chosen vertex
	vector<Node*> movableNodes;

	for (auto n : invMap[fromVertex]) {
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
				if (!toVertex->isAncestralTo(CMM.getMap(c)->getHost(c)))
					movable = false;
			}

			if (movable)
				movableNodes.push_back(n);
		}
	}

	//reject if not enough movable nodes
	if (movableNodes.size() < 2) {
		DEBUG(cout << "Not enough movable nodes at chosen vertex" << endl;)
		return;
	}

	//choose number of nodes to move
	//choose nodes to move
	//move nodes
	//accept/reject
	






	


	Node* top = s->getRoot();
	Node* restrict;
	Node* bottom;
	set<Node*> childHosts;
	//find all vertices with which can be images for the duplications
	for (auto n : invMap[fromVertex]) {
		if (CMM.getMap(n)->getEvent(n) == duplication) {
			//if parent of n is not at fromVertex get upper limit
			if (n->hasParent() && CMM.getMap(n->getParent())->getHost(n->getParent()) != fromVertex) {
				//cannot go higher than image of parent
				restrict = CMM.getMap(n->getParent())->getHost(n->getParent()); 
				//move down 1 if speciation
				if (CMM.getMap(n->getParent())->getEvent(n->getParent()) != duplication) {
					Node* temp = fromVertex;
					while (temp->getParent() != restrict)
						temp = temp->getParent();
					restrict = temp;
				}
				top = s->min(top, restrict);
			}

			//if children of n are not at fromVertex OR are speciations get lower limit
			for (Node* c = n->getFirstChild(); c != nullptr; c = c->getSibling()) {
				if (CMM.getMap(c)->getHost(c) != fromVertex || CMM.getMap(c)->getEvent(c) != duplication) {
					childHosts.insert(CMM.getMap(c)->getHost(c));
				}
			}
			
		}
	}
	bottom = s->LCA(childHosts);

	//calculate weights, total
	int dupsToMove = 0;
	for (auto n : invMap[fromVertex]) {
		if (CMM.getMap(n)->getEvent(n) == duplication) {
			dupsToMove++;

			if (!n->hasParent()) 
				dupsToMove++;
		}
	}

	EventCount ec;
	Node* temp;
	double total(1.0);

	//go from bottom to top INCLUSIVE
	vector<pair<Node*,double>> toVertices;
	toVertices.push_back(make_pair(fromVertex,1.0));
	for (Node* v = bottom; true; v = v->getParent()) {
		int dupsAtV = 0;
		for (auto n : invMap[v]) 
			if (CMM.getMap(n)->getEvent(n) == duplication)
				dupsAtV++;
		if (dupsAtV == 0) {
			ec.clear();
			if (v->isAncestralTo(fromVertex))
				ec.losses = dupsToMove*s->getDistUp(fromVertex,v);
			else if (fromVertex->isAncestralTo(v))
				ec.losses = -dupsToMove*s->getDistUp(v,fromVertex);

			toVertices.push_back(make_pair(v,exp(-CSD(ec)/T)));
			total += exp(-CSD(ec)/T);
		}

		if (v == top)
			break;
	}

	if (toVertices.size() == 1) {
		DEBUG(cout << "Nowhere to go" << endl;)
		return;
	}

	DEBUG(cout << "toVertices: "; for (auto v : toVertices) cout << v.first->getLabel() << "," << v.second << " "; cout << endl;)
	DEBUG(cout << "total is " << total << endl;)
	
	//resample vertex
	double r = dran(total);
	/***********************************
	 * Now the sampling!
	 ***********************************/
	for (auto v : toVertices) {
		DEBUG(cout << "r = " << r << ", weight = " << v.second << endl;)

		if (r <= v.second) {
			// sample this one

			//recount :(
			EventCount ec;
			if (v.first->isAncestralTo(fromVertex))
				ec.losses = dupsToMove*s->getDistUp(fromVertex,v.first);
			else if (fromVertex->isAncestralTo(v.first))
				ec.losses = -dupsToMove*s->getDistUp(v.first,fromVertex);
			DEBUG(cout << "Event count change is " << ec << endl;)
			DEBUG(cout << "Current event count is " << CMM.countEvents() << endl;)
			ec += CMM.countEvents();
			CMM.setCurrentEventCount(ec);
			//event count is wrong
			DEBUG(cout << "New event count is " << ec << endl;)

			DEBUG(cout << "Moving from " << fromVertex->getLabel() << " to " << v.first->getLabel() << endl;)

			set<Node*> toMove(invMap[fromVertex]);
			for (auto n : toMove) 
				if (CMM.getMap(n)->getEvent(n) == duplication) {
					DEBUG(cout << "Moving " << n->getLabel() << " to " << v.first->getLabel() << endl;)
					CMM.getMap(n)->moveToHost(n, v.first, duplication);
					CMM.movePToHost(n, fromVertex, v.first);
				}

			if (v.first != fromVertex)
				DEBUG(cout << "SingleVertex move succeeded!" << endl;)

			DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)
			DEBUG(for(auto v : s->getVertices()) { cout << v.second->getLabel() << ": "; for (auto n : invMap[v.second]) cout << n->getLabel() << " "; cout << endl; })

			break;
		}
		r -= v.second;
	}	

	DEBUG(cout << endl;)
}

}; // end of namespace
