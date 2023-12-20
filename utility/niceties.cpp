/*
 * niceties.cpp
 *
 *  Created on: 2 Aug 2023
 *      Author: mac
 */

#include "niceties.h"

void advance_cursor() {
  static int pos = 0;
  char cursor[4] = {'/','-','\\','|'};
  printf("%c\b", cursor[pos]);
  fflush(stdout);
  pos = (pos+1) % 4;
}

void update_message(std::string str) {
	static unsigned int lastMessageLength(0);
	for (unsigned int i(0); i < lastMessageLength; ++i) {
		printf("\b");
	}
//	printf(str.c_str());
	lastMessageLength = str.length();
}
