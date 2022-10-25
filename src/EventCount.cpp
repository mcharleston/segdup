/*
 * EventCount.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: mac
 */

#include "EventCount.h"

namespace segdup {

EventCount::~EventCount() {
}

EventCount::EventCount(const EventCount &other) {
	codivs = other.codivs;
	dups = other.dups;
	losses = other.losses;
}

EventCount& EventCount::operator+=(const EventCount &other) {
	codivs += other.codivs;
	dups += other.dups;
	losses += other.losses;
	return *this;
}

EventCount& EventCount::operator=(const EventCount &other) {
	codivs = other.codivs;
	dups = other.dups;
	losses = other.losses;
	return *this;
}

std::ostream& operator<<(std::ostream& os, const EventCount& ec) {
	os << "codivs=" << ec.codivs << "; dups=" << ec.dups << "; losses=" << ec.losses;
	return os;
}

} /* namespace segdup */
