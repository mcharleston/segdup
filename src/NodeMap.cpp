/*
 * NodeMap.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: mac
 */

#include <cassert>
#include <string>
#include <vector>

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split

#include "../utility/debugging.h"
#include "NodeMap.h"
#include "Tree.h"


extern bool _debugging;

namespace segdup {

using namespace std;

NodeMap::NodeMap(initializer_list<pair<Node*, Node*>> leaves) {
	try {
		P = leaves.begin()->first->getTree();
		H = leaves.begin()->second->getTree();
		DEBUG(cout << "NodeMap(initializer_list) constructor: setting Host tree = " << H
				<< " and P = " << P << endl);
		for (auto phi : leaves) {
			data[phi.first] = pair<Node*, short>(phi.second, 0);
		}
	} catch (const exception& e) {
		cout << e.what() << endl;
	}
}


NodeMap::NodeMap(Tree* inH, Tree* inP, string assocStr) : H(inH), P(inP) {
	// read in associations as node names.  Only works once the trees have been created!
	// NOTE that the Host and Parasite trees MUST be given in this order. They are the same type.
	// NOTE too that the associations are written as p:h -- yes this is the other way around gtom above but it's traditional...
	// format:
	// <assocs> ::= <paraname> ':' <hostname> { ',' <assocs> }

	bool _debugging(false);
	DEBUG(cout << "Constructing NodeMap from an input string\n");
	std::vector<std::string> words;
	boost::split(words, assocStr, boost::is_any_of(":, \t"), boost::token_compress_on);
	for (vector<string>::iterator sit = words.begin(); sit != words.end(); ++sit) {
		Node* p = (*P)[*sit];
		++sit;
		Node* h = (*H)[*sit];
		data[p] = pair<Node*, short>(h, 0);
		DEBUG(cout << "parasite " << p->getLabel() << " is on host " << h->getLabel() << endl)
	}
}

NodeMap::~NodeMap() {
}

Node* NodeMap::getImage(Node* p) {
	try {
		return data.at(p).first;
	} catch (const exception& e) {
		cout << e.what();
		return nullptr;
	}
}

Node* NodeMap::getImage(const std::string& pLabel) {
	try {
		return data.at((*P)[pLabel]).first;
	} catch (const exception& e) {
		cout << e.what();
		return nullptr;
	}
}

/**
 * Returns the (host node, association type) pair for the parasite
 */
associationtype& NodeMap::operator[](const std::string& paraLabel) {
	return data[(*P)[paraLabel]];
}

NodeMap& NodeMap::operator=(const NodeMap &other) {
	H = other.H;
	P = other.P;
	data = other.data;
	return *this;
}

NodeMap::NodeMap(const NodeMap &other) {
	data = other.data;
	H = other.H;	// host / species tree
	P = other.P;	// parasite / gene tree
}

void NodeMap::setPHAssociation(const std::string& pstr, const std::string& hstr, short a) {
	bool _debugging(false);
	DEBUG(cout << "setPHAssociation(pstr, hstr, a): P = " << P << "; H = " << H << endl);
	assert(P != nullptr);
	assert(H != nullptr);
	DEBUG(cout << "pstr = " << pstr <<"; hstr = " << hstr << endl);
	Node *p = P->at(pstr);
	Node *h = H->at(hstr);
	assert(p);
	assert(h);
	setPHAssociation(p, h, a);
}

ostream& operator<<(ostream& os, NodeMap& NM) {
	os << "P -> H" << endl;
	for (auto d : NM.getData()) {
		os << d.first->getLabel() << " -> (" << d.second.first->getLabel();
		os << ',' << d.second.second << ')' << endl;
	}
	return os;
}


} /* namespace segdup */
