/*
 * SingleNodeMove.cpp
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
#include "SingleNodeMove.h"

using namespace std;

namespace segdup {

SingleNodeMove::SingleNodeMove() {
	// TODO Auto-generated constructor stub

}

SingleNodeMove::~SingleNodeMove() {
	// TODO Auto-generated destructor stub
}

void SingleNodeMove::apply(CophyMultiMap& CMM, double T) {
	bool _debugging(false);
	int size = CMM.getAllMoveableNodes().size();
	auto mpr = (CMM.getAllMoveableNodes())[iran(size)];	// pick at random from allMoveableNodes

	Node* p(mpr.first);

	CophyMap *M = mpr.second;
	EventCount& currentEventCount = CMM.countEvents();
	EventCount ec;
	unsigned int nullMoves;
	std::set<Contender> neighbours;

	double minScore(1e10);	
	unsigned int numNeighbours(0);

	set<pair<Node*, eventType>> nextImages = M->calcAvailableNewHosts(p);
		DEBUG(
			cout << "Node " << p->getLabel() << " has possible images { ";
			for (auto a : nextImages) {
				cout << eventSymbol[a.second] << a.first->getLabel() << " ";
			}
			cout << "};" << endl;
		);
		for (auto a : nextImages) {
			Node* nuHost = a.first;
			eventType nuEvent = a.second;
			if (nuHost == M->getHost(p) && nuEvent == M->getEvent(p)) {
				double relativeSamplingProbability = exp(-CSD(currentEventCount) / (1.0*T));
				Contender noChange( currentEventCount, CSD(currentEventCount), p, nuHost, nuEvent, M );
//						DEBUG(cout << "XXX entering CSD for no move: " << CSD(currentEventCount) << endl);
				minScore = std::min(minScore, CSD(currentEventCount));

				noChange._noMove = true;
				noChange.setLabel("leaving " + p->getLabel() + " on [" + eventSymbol[nuEvent] + "]" + nuHost->getLabel());

				neighbours.insert(noChange);
				++nullMoves;
				continue;
			}
			ec.clear();
			Node* oldHost = M->getHost(p);
			eventType oldEvent = M->getEvent(p);
			// Get cost of move, without recalculating everything:

			// Any change in the number of codivergences?
			DEBUG(cout << "p=" << p->getLabel() << "; oldHost=" << oldHost->getLabel()
					<< "; newHost=" << a.first->getLabel() << endl);
			ec.codivs = (nuEvent == codivergence) ? 1 : 0;
			ec.codivs -= (oldEvent == codivergence) ? 1 : 0;

			// Any change in number of losses?
			int lossFactor = (p->hasParent()) ? 1 : 2;	// root of gene tree has no ancestral branch!
			if (oldHost->isAncestralTo(nuHost)) {
				DEBUG(cout << "old host is ancestral to new host" << endl);
				ec.losses = -lossFactor*(nuHost->getTree()->getDistUp(nuHost, oldHost));
			} else if (nuHost->isAncestralTo(oldHost)) {
				DEBUG(cout << "new host is ancestral to old host" << endl);
				ec.losses = lossFactor*(nuHost->getTree()->getDistUp(oldHost, nuHost));
			} else {
				DEBUG(cout << "old and new hosts are not ancestrally comparable" << endl);
			}
			if (oldEvent == duplication && nuEvent == codivergence) {
				ec.losses -= 2;
			} else if (nuEvent == duplication && oldEvent == codivergence) {
				ec.losses += 2;
			}

			// Any change in number of duplications?
			int oldJointDupHeightHere(CMM.calcCombinedDuplicationHeight(oldHost));
			int oldJointDupHeightThere(CMM.calcCombinedDuplicationHeight(nuHost));


			DEBUG(cout << "oldJointDupHeightHere = " << oldJointDupHeightHere << endl);
			DEBUG(cout << "oldJointDupHeightThere = " << oldJointDupHeightThere << endl);
			M->moveToHost(p, nuHost, nuEvent);
			CMM.movePToHost(p,oldHost,nuHost);	// XXX are both these function calls necessary?
			int nuJointDupHeightHere(CMM.calcCombinedDuplicationHeight(oldHost));
			int nuJointDupHeightThere(CMM.calcCombinedDuplicationHeight(nuHost));
			DEBUG(cout << "nuJointDupHeightHere = " << nuJointDupHeightHere << endl);
			DEBUG(cout << "nuJointDupHeightThere = " << nuJointDupHeightThere << endl);
			M->moveToHost(p, oldHost, oldEvent);
			CMM.movePToHost(p, nuHost, oldHost);	// XXX are both these function calls necessary?
			ec.dups = nuJointDupHeightThere - oldJointDupHeightThere;
			if (nuHost != oldHost) {
				ec.dups += nuJointDupHeightHere - oldJointDupHeightHere;
			}
			EventCount dec(ec);
			ec += currentEventCount;
			if (ec.dups < 0) {
				cout << "CRITICAL FAILURE!" << endl;
				cout << "previous event count = " << currentEventCount << endl;
				cout << "delta ec = " << dec << ";\t";
				cout << "ec+originalEventCount = " << ec << ";\t";
				exit(-1);
			}
			Association ass(p, nuHost, nuEvent);
//					double relativeSamplingProbability = exp(-CSD(ec) / (1.0*T));	// XXX test just using dec not ec here
			minScore = min(minScore, CSD(ec));
//					DEBUG(
//							cout << "SCORE = " << relativeSamplingProbability << endl;
//					);
			Contender con( ec, CSD(ec), p, nuHost, nuEvent, M );
//					DEBUG(cout << "XXX setting contender score to " << CSD(ec) << endl);
			if (con.getEventCount().dups < 0) {
				cout << "CRITICAL FAILURE!" << endl;
				cout << "Contender event count = " << con.getEventCount().dups << endl;
				exit(-1);
			}
			string label = "moving " + p->getLabel() + " to [" + eventSymbol[nuEvent] + "]" + nuHost->getLabel();
			con.setLabel(label);
			neighbours.insert(con);
//					cerr << t << ',' << T << ',' << p << ',' << con.getScore() << endl;
		}
		double d;
		double total(0.0);
		set<Contender> adjustedNeighbours;
		for (auto nei : neighbours) {
			d = nei.getScore() - minScore; // XXX This is a fudge...
//			cerr << "d=" << d << endl;
			DEBUG(cout << "relative CSD for this contender = " << d << endl);
			nei.setScore(exp(-d / (1.0*T)));
			total += nei.getScore();
			adjustedNeighbours.insert(nei);
			DEBUG(cout << "T = " << T << "; new score = " << nei.getScore() << endl);
			DEBUG(cout << "new total " << total << endl);
		}
//		for (auto nei : neighbours) {
//		}
		numNeighbours += neighbours.size();
		DEBUG(cout << "total before dran() = " << total << endl);
		double r = dran(total);
		DEBUG(cout << "total after dran() = " << total << endl);
		/***********************************
		 * Now the sampling!
		 ***********************************/
		DEBUG(cout << "initial r = " << r << " from U[0, " << total << "]" << endl);
		for (auto nei : adjustedNeighbours) {
			if (r <= nei.getScore()) {
				DEBUG(cout << "r=" << r << "; score=" << nei.getScore() << endl);
				if (nei._noMove) {
					DEBUG(cout << "Selected move: No change (probability = " << (nei.getScore()/total) << ")" << endl);
				} else {
					// sample this one
					DEBUG(cout << "Selected move: " << nei.getLabel() << " (rel. probability = " << nei.getScore() << ")" << endl);
					CMM.movePToHost(nei.getParasite(), CMM.getMap(nei.getParasite())->getHost(nei.getParasite()), nei.getHost());
					nei.getMap()->moveToHost(nei.getParasite(), nei.getHost(), nei.getEvent());
					DEBUG(nei.getMap()->checkValidHostOrdering());
				}

				CMM.setCurrentEventCount(nei.getEventCount());

				DEBUG(string str; CMM.toCompactString(str); cout << str << endl;)
				DEBUG(inversenodemap& invMap = CMM.getInverseMap(); for(auto v : CMM.getHostTree()->getVertices()) { cout << v.second->getLabel() << ": "; for (auto n : invMap[v.second]) cout << n->getLabel() << " "; cout << endl; })

				break;
			}
			DEBUG(cout << r << ' ');
			r -= nei.getScore();
		}
		DEBUG(cout << endl);


}

}; // end of namespace
