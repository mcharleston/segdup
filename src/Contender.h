/*
 * Contender.h
 *
 *  Created on: 29 Nov 2022
 *      Author: mac
 */

#ifndef SRC_CONTENDER_H_
#define SRC_CONTENDER_H_

#include <iostream>
#include <string>

#include "CophyMap.h"
//#include "NodeMap.h"

namespace segdup {

/**
 * A Contender is a potential move we could make on a CophyMap, with its probability of being selected.
 * Contenders are ranked in decreasing order of this probability to make sampling from a collection of
 * them efficient.
 */
class Contender {
private:
	EventCount ec;
	double score;
	std::string label;
	Association assoc;
	CophyMap* M;
public:
	bool _noMove;
//	Contender(double d, const Association& as) : score(d), label(""), assoc(as), M(nullptr) {}
	Contender(EventCount ecount, double d, Node* p, Node* h, eventType e, CophyMap* m)
			: ec(ecount), score(d), assoc(p, h, e), M(m), _noMove(false)  {}

//	Contender(double d, std::string s) : score(d), label(s) {}
	Contender(const Contender& c) : ec(c.ec), score(c.score), label(c.label), assoc(c.assoc), M(c.M), _noMove(c._noMove) {}
	virtual ~Contender();

	eventType getEvent() const { return assoc.e; }
	EventCount getEventCount() const { return ec; }
	Node* getHost() { return assoc.h; }
	Node* getHost() const { return assoc.h; }
	std::string getLabel() const { return label; }
	CophyMap* getMap() { return M; }
	Node* getParasite() { return assoc.p; }
	Node* getParasite() const { return assoc.p; }
	double getScore() { return score; }
	double getScore() const { return score; }

	void operator*=(double d) { score *= d; }
	void operator/=(double d) { score /= d; }

	friend bool operator<(const Contender& a, const Contender& b);
	friend bool operator>(const Contender& a, const Contender& b);
	friend bool operator==(const Contender& a, const Contender& b);
	friend std::ostream& operator<<(std::ostream &os, const Contender& c);

	void setLabel(const std::string& str) { label = str; }
	void setScore(double d) { score = d; }
};

} /* namespace segdup */

#endif /* SRC_CONTENDER_H_ */
