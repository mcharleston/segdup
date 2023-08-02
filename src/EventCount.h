/*
 * EventCount.h
 *
 *  Created on: 14 Oct 2022
 *      Author: mac
 */

#ifndef SRC_EVENTCOUNT_H_
#define SRC_EVENTCOUNT_H_

#include <iostream>


namespace segdup {

class EventCount {
	friend class CophyMap;
public:
	int codivs;
	int dups;
	int losses;
	EventCount() : codivs(0), dups(0), losses(0) {}
	virtual ~EventCount();
	EventCount(const EventCount &other);
	void clear() { codivs = dups = losses = 0; }
	EventCount& operator+=(const EventCount &other);
	EventCount& operator=(const EventCount &other);
	friend std::ostream& operator<<(std::ostream& os, const EventCount& ec);
};


} /* namespace segdup */

#endif /* SRC_EVENTCOUNT_H_ */
