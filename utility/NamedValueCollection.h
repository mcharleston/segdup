/*
 * NamedValueCollection.h
 *
 *  Created on: 27 June 2022
 *      Author: mac
 */

#ifndef UTILITY_NAMEDVALUECOLLECTION_H_
#define UTILITY_NAMEDVALUECOLLECTION_H_

#include <string>
#include <map>
#include <stdio.h>
#include <iostream>

#include <boost/range.hpp>

#include "../utility/appexception.h"
#include "../utility/debugging.h"

/**
 * I just don't need the complexity of a property_map thing; I'm happy with a simple collection that
 * is a collection of maps of string -> <type> with <type> being one of int, float, bool, string.
 */
class NamedValueCollection {
private:
	std::map<std::string, bool> boolValues;
	std::map<std::string, int> intValues;
	std::map<std::string, double> realValues;
	std::map<std::string, std::string> strValues;
public:
	char keyValSep, fieldSep;
	NamedValueCollection() : keyValSep(':'), fieldSep('\t') {}
	virtual ~NamedValueCollection();
	NamedValueCollection(const NamedValueCollection &other);

	void decr(const std::string& str);
	void incr(const std::string& str);

	bool getBool(const std::string& s);
	int getInt(const std::string& s);
	double getReal(const std::string& s);
	std::string getString(const std::string& s);

	NamedValueCollection& operator=(const NamedValueCollection &other);

	std::map<std::string, bool>& theBools();
	const std::map<std::string, bool>& theBools() const;
	std::map<std::string, int>& theInts();
	const std::map<std::string, int>& theInts() const;
	std::map<std::string, double>& theReals();
	const std::map<std::string, double>& theReals() const;
	std::map<std::string, std::string>& theStrings();
	const std::map<std::string, std::string>& theStrings() const;

	void setBool(const std::string& s, bool b);
	void setInt(const std::string& s, int i);
	void setReal(const std::string& s, double r);
	void setString(const std::string& s, std::string str);

	void toggle(const std::string& str);

};

std::ostream& operator<<(std::ostream &os, NamedValueCollection& nvc);

#endif /* UTILITY_NAMEDVALUECOLLECTION_H_ */
