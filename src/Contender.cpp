/*
 * Contender.cpp
 *
 *  Created on: 29 Nov 2022
 *      Author: mac
 */

#include "Contender.h"

namespace segdup {

Contender::~Contender() {
	// TODO Auto-generated destructor stub
}

bool operator<(const Contender &a, const Contender &b) {
	if (a.score > b.score) {
		return true;
	}
	if (a.getParasite() < b.getParasite()) {
		return true;
	}
	if (a.getHost() < b.getHost()) {
		return true;
	}
	if (a.getParasite() < b.getParasite()) {
		return true;
	}
	if (a.getEvent() == b.getEvent()) {
		return (a.label < b.label);
	}
	return true;
}
bool operator>(const Contender &a, const Contender &b) {
	return ( b < a );
}

bool operator==(const Contender &a, const Contender &b) {
	if (a.score != b.score) {
		return false;
	}
	if (a.assoc.p != b.assoc.p) {
		return false;
	}
	if (a.assoc.h != b.assoc.h) {
		return false;
	}
	return (a.label == b.label);
}

std::ostream& operator<<(std::ostream& os, const Contender& c) {
	os << c.label << ':' << c.score << " for " << eventSymbol[c.assoc.e]
		<< c.assoc.p->getLabel() << ':' << c.assoc.h->getLabel() ;
	return os;
}

} /* namespace segdup */
