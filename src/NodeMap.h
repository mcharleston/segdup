/*
 * NodeMap.h
 *
 *  Created on: 13 Oct 2022
 *      Author: mac
 */

#ifndef SRC_NODEMAP_H_
#define SRC_NODEMAP_H_

#include <initializer_list>
#include <map>
#include <set>
#include <string>

namespace segdup {

class Node;
class Tree;

typedef std::pair<Node*, short> associationtype;	// host x association type
typedef std::map<Node*, associationtype> nodemaptype;	// parasite x (host x association type)
typedef std::map<Node*, std::set<associationtype> > wideassociationmap;
typedef std::map<Node*, std::set<Node*>> inversenodemap;

/**
 * Note that this formulation precludes having more than one host for a parasite!
 */

class NodeMap {
private:
	nodemaptype data;
	Tree* H;	// host / species tree
	Tree* P;	// parasite / gene tree
public:
	NodeMap() : H(nullptr), P(nullptr) {}
	NodeMap(std::initializer_list<std::pair<Node*, Node*> > leaves);
	NodeMap(Tree* H, Tree* P, std::string assocStr);
	NodeMap(const NodeMap &other);
	virtual ~NodeMap();

	associationtype& at(Node* n) { return data.at(n); }
	inline void clear() { data.clear(); }

	nodemaptype& getData() { return data; }
	Tree* getHostTree() { return H; }	// should I also have a getSpeciesTree and a getGeneTree?
	Tree* getParasiteTree() { return P; }
	Node* getImage(Node* p);
	Node* getImage(const std::string& pLabel);

	inline associationtype& operator[](Node* n) { return data[n]; }
	associationtype& operator[](const std::string& paraLabel);
	NodeMap& operator=(const NodeMap &other);

	void setPHAssociation(Node* p, Node* h, short a) { data[p] = std::pair<Node*, short>(h, a); }
	void setPHAssociation(const std::string& pstr, const std::string& hstr, short a);
	void setHostTree(Tree* hptr) { H = hptr; }
	void setParasiteTree(Tree* pptr) { P = pptr; }

};

std::ostream& operator<<(std::ostream& os, NodeMap& NM);

} /* namespace segdup */

#endif /* SRC_NODEMAP_H_ */
