/*
 * EmptyMove.cpp
 *
 *  Created on: 14 Dec 2023
 *      Author: mac
 */

#include <iostream>
#include <cmath>

#include "../utility/debugging.h"
#include "../utility/myrandom.h"

#include "Contender.h"
#include "Node.h"
#include "SegDup.h"
#include "EmptyMove.h"

using namespace std;

namespace segdup {

EmptyMove::EmptyMove() {
	// TODO Auto-generated constructor stub

}

EmptyMove::~EmptyMove() {
	// TODO Auto-generated destructor stub
}

void EmptyMove::apply(CophyMultiMap& CMM, double T) {
	bool _debugging(true);

	DEBUG(cout << "Attempting empty move" << endl;)
	DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)
	//find all species tree vertices with >1 duplications mapped to them

	Tree* s = CMM.getHostTree();
	inversenodemap& invMap = CMM.getInverseMap();
	DEBUG(
		for(auto v : s->getVertices()) {
			cout << v.second->getLabel() << ": ";
			for (auto n : invMap[v.second]) {
				cout << n->getLabel() << " ";
			}
			cout << endl;
		}
	)
	vector<Node*> fromVertices;
	for (auto v : s->getVertices()) {
		int dupsAtV = 0;
		for (auto n : invMap[v.second]) 
			if (CMM.getMap(n)->getEvent(n) == duplication)
				dupsAtV++;

		if (dupsAtV > 1)
			fromVertices.push_back(v.second);
	}
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
	
	//choose a vertex from fromVertices
	Node* fromVertex = fromVertices[iran(fromVertices.size())];

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
				DEBUG(cout << "Empty move succeeded!" << endl;)

			DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)
			DEBUG(for(auto v : s->getVertices()) { cout << v.second->getLabel() << ": "; for (auto n : invMap[v.second]) cout << n->getLabel() << " "; cout << endl; })

			break;
		}
		r -= v.second;
	}	

	DEBUG(cout << endl;)
}

}; // end of namespace
