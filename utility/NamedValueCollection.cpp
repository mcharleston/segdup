/*
 * NamedValueCollection.cpp
 *
 *  Created on: 27 June 2022
 *      Author: mac
 */

#include <sstream>
#include "NamedValueCollection.h"

NamedValueCollection::~NamedValueCollection() {
	// TODO Auto-generated destructor stub
}

NamedValueCollection::NamedValueCollection(const NamedValueCollection &other) {
	// TODO Auto-generated constructor stub
	keyValSep = ':';
	fieldSep = '\t';
}

void NamedValueCollection::decr(const std::string& str) {
	if (intValues.find(str) == intValues.end()) {
		std::stringstream ss;
		ss << "Cannot decrement this value \'" << str << "\' as it is not an integer in the collection.";
		throw new app_exception(ss.str());
	}
	setInt(str, intValues.at(str)-1);
}

bool NamedValueCollection::getBool(const std::string& s) {
	return boolValues.at(s);
}
int NamedValueCollection::getInt(const std::string& s) {
	return intValues.at(s);
}
double NamedValueCollection::getReal(const std::string& s) {
	return realValues.at(s);
}
std::string NamedValueCollection::getString(const std::string& s) {
	return strValues.at(s);
}

void NamedValueCollection::incr(const std::string& str) {
	if (intValues.find(str) == intValues.end()) {
		std::stringstream ss;
		ss << "Cannot increment this value \'" << str << "\' as it is not an integer in the collection.";
		throw new app_exception(ss.str());
	}
	setInt(str, intValues.at(str)+1);
}

NamedValueCollection& NamedValueCollection::operator=(
		const NamedValueCollection &other) {
	// TODO Auto-generated method stub
	return *this;
}

std::map<std::string, bool>& NamedValueCollection::theBools() {
	return boolValues;
}
const std::map<std::string, bool>& NamedValueCollection::theBools() const {
	return boolValues;
}
std::map<std::string, int>& NamedValueCollection::theInts() {
	return intValues;
}
const std::map<std::string, int>& NamedValueCollection::theInts() const {
	return intValues;
}
std::map<std::string, double>& NamedValueCollection::theReals() {
	return realValues;
}
const std::map<std::string, double>& NamedValueCollection::theReals() const {
	return realValues;
}
std::map<std::string, std::string>& NamedValueCollection::theStrings() {
	return strValues;
}
const std::map<std::string, std::string>& NamedValueCollection::theStrings() const {
	return strValues;
}

void NamedValueCollection::setBool(const std::string& s, bool b) {
	boolValues[s] = b;
}
void NamedValueCollection::setInt(const std::string& s, int i) {
	intValues[s] = i;
}
void NamedValueCollection::setReal(const std::string& s, double r) {
	realValues[s] = r;
}
void NamedValueCollection::setString(const std::string& s, std::string str) {
	strValues[s] = str;
}

void NamedValueCollection::toggle(const std::string& str) {
	if (boolValues.find(str) == boolValues.end()) {
		std::stringstream ss;
		ss << "Cannot toggle this value \'" << str << "\' as it is not a Boolean value in the collection.";
		throw new app_exception(ss.str());
	}
	setBool(str, !boolValues.at(str));
}

std::ostream& operator<<(std::ostream &os, NamedValueCollection& nvc) {
	std::map<std::string, bool>::iterator bi = nvc.theBools().begin(), bend = nvc.theBools().end();
	for ( ; bi != bend; ++bi) {
		os << bi->first << nvc.keyValSep << bi->second << nvc.fieldSep;
	}
	std::map<std::string, int>::iterator ii = nvc.theInts().begin(), iend = nvc.theInts().end();
	for ( ; ii != iend; ++ii) {
		os << ii->first << nvc.keyValSep << ii->second << nvc.fieldSep;
	}
	std::map<std::string, double>::iterator di = nvc.theReals().begin(), dend = nvc.theReals().end();
	for ( ; di != dend; ++di) {
		os << di->first << nvc.keyValSep << di->second << nvc.fieldSep;
	}
	std::map<std::string, std::string>::iterator si = nvc.theStrings().begin(), send = nvc.theStrings().end();
	for ( ; si != send; ++si) {
		os << si->first << nvc.keyValSep << si->second << nvc.fieldSep;
	}
	return os;
}

