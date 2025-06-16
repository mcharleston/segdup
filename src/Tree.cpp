/*
 * Tree.cpp
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#include <assert.h>
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "SegDupParser.h"
#include "Tree.h"

using namespace std;

namespace segdup {

Tree::Tree(char pre, std::string str) : root(nullptr), labelSpace(0), numVertices(-1), prefix(pre) {
	if (str[0] == '(') {
		constructFromNewickString(str);
		calculateHeights(root);
		calcAncestry();
		_showInfo = false;
	} else {
		cerr << "Culprit:\n" << str << endl;
		throw new app_exception("Tree constructor[1] with a single string argument is expecting Newick format tree description.");
	}
}
Tree::Tree(std::string str) : root(nullptr), labelSpace(0), numVertices(-1) {
	prefix = 'v';
	if (str[0] == '(') {
		constructFromNewickString(str);
		calculateHeights(root);
		calcAncestry();
		_showInfo = false;
	} else {
		throw new app_exception("Tree constructor[2] with a single string argument is expecting Newick format tree description.");
	}
}

/**
 * @pre: tree is complete
 * @post: distUp<x,y> contains number of vertices UP from x to y.
 * distUp<x,y> < 0 iff they are unrelated
 * distUp<x,y> == 0 iff x==y
 * distUp<x,y> > 0 iff y is ancestral to x and is this number of EDGEs above it.
 */
void Tree::calcAncestry() {
	distUp.clear();
	if (V.size() == 0) {
		gatherVertices();
	}
	for (auto v : V) {
		Node* p = v.second->getParent();
		int d(0);
		while (p != nullptr) {
			distUp[{ v.second, p }] = ++d;
			p = p->getParent();
		}
	}
}

const int kMediumStringLength = 65536;
char treebuf[kMediumStringLength];	// XXX magic number

/**
 * Traverse the tree from the root node and print it in ASCII characters to the ostream
 */
void Tree::compressTraverseWrite(ostream& os) {
	bool _debugging(false);
	DEBUG(cout << "starting compressTraverseWrite" << endl);
	DEBUG(if (root==nullptr) { cout << "Oops" << endl; } );
	if (root->getHeight() < 0) {
		DEBUG(cout << "root->getHeight() negative so calculating heights" << endl);
		calculateHeights(root);
	}
	DEBUG(cout << "Address of info map = " << &info << endl);
	bool _oldShowInfo(_showInfo);
	if (info == nullptr) {
		_showInfo = false;
	}
	labelSpace = getMaxLabelWidth(root);
	compressTraverseWrite(os, root);
	_showInfo = _oldShowInfo;
}
void Tree::compressTraverseWrite(ostream & os, Node* p) {
	int
		spacing;
	char
		hline = '-',
		link[1024],
		vline = '|';
	string nodeLabel;
	if (p != nullptr) {
		Node* q = p->getParent();
		int startpos = 0, stoppos = 0;
		spacing = 3 + labelSpace;
/* q is the parent. Height of a vertex is the distance from the root as it is scaled
 * from the time of the vertex, with time(root) = 1 and time(tips) = 1 typically.
 * Thus we have startpos = q->height and stoppos = p->height.  Let's give that a try...
 */
		bool _debugging(false);
//		DEBUG(cout << "starting compressTraverseWrite" << endl);
		if (p != root) {
			startpos = spacing*(root->getHeight()+1 - q->getHeight());
			stoppos = spacing*(root->getHeight()+1 - p->getHeight());
		}
		int i;	// handy counter
		for (i = startpos + 1; i <= stoppos; ++i)
			treebuf[i] = hline;
		std::strncpy(link, "--------------------------------", spacing);
		DEBUG(cout << "p = " << p->getLabel() << endl);
		DEBUG(if (q != nullptr) { cout << "parent(p) = " << q->getLabel() << endl; } );

		if (p->isFirstChild()) {
			nodeLabel = "";
//			if (_showInfo && info.count(q) > 0) {
//				nodeLabel = (q->getEvent() == codivergence) ? "[<]" : "[=]"; // XXX This is a problem!
//				// XXX I need to store the event separately from the node, in the cophymap
//			}
//			nodeLabel += q->getLabel();
			if (_showInfo && (info->count(q) > 0)) {
				nodeLabel = info->at(q);
			} else {
				nodeLabel = q->getLabel();
			}
			std::memcpy(link, nodeLabel.c_str(), nodeLabel.length());
		} else {
			link[0] = '+';
			while (q != nullptr) {
				if (q->getSibling() != nullptr) {
					i = spacing*(root->getHeight()+1 - q->getParent()->getHeight());
					treebuf[i++] = vline;
					for ( ; i < spacing; i++) {
						treebuf[i] = ' ';
					}
				}
				q = q->getParent();
			}
		}
		for (int x = startpos; x < startpos + spacing; ++x) {
			treebuf[x] = link[x-startpos];
		}

		// if p is a leaf then we need to output the line buffer:
		if (p->isLeaf()) {
			treebuf[spacing*(root->getHeight()+1)] = '\0';
			if (_showInfo && (info->count(p) > 0)) {
				nodeLabel = info->at(p);
			} else {
				nodeLabel = p->getLabel();
			}
//			string label(p->getLabel());
//			if (_showInfo && (info.count(p) > 0)) {
//				label += ":";
//				label += info.at(p);
//			}
			os << treebuf << " " << nodeLabel << "\n";
			os.flush();
			// we've output the string now, so clear the buffer:
			for (i = 0; i < kMediumStringLength; i++)
				treebuf[i] = ' ';
		} else {
			compressTraverseWrite(os, p->getFirstChild());
		}
		if (p->getSibling() != nullptr)
			compressTraverseWrite(os, p->getSibling());
	}
}

void Tree::constructFromNewickString(std::string str) {
	parsing::TokenList TL;
	istringstream instr(str);
	TL.tokenize(instr);
	parsing::SegDupParser p(TL);
	root = new Node("*");
	TL.reset();
	numVertices = 1;
	p.parseNewickSubtree(root, prefix);
	gatherVertices();
	for (auto n : V) {
		n.second->setTree(this);
	}
}

int Tree::getDistUp(Node* lower, Node* upper) {
	if (distUp.size() == 0) {
		calcAncestry();
	}
	return distUp[std::pair<Node*, Node*>(lower, upper)];
}


void Tree::gatherVertices() {
	V.clear();
	if (root == nullptr) {
		return;	// there are no vertices!
	}
	root->addSubtreeVertices(V);
}

int Tree::getMaxLabelWidth(Node *v) {
	int length = v->getLabel().length();
	if (_showInfo && (info->count(v) > 0)) {
		length = info->at(v).length();
	}
	Node* child = v->getFirstChild();
	while (child != nullptr) {
		length= std::max<int>(length, getMaxLabelWidth(child));
		child = child->getSibling();
	}
	return length;
}

bool Tree::isAncestralTo(Node* x, Node*y) {
	if (x->getTree() != y->getTree()) {
		throw new app_exception("Comparing nodes from two different trees for ancestry -- this doesn't make sense.");
	}
	if (x->isLeaf()) {
		return false;
	}
	return (distUp[std::pair<Node*, Node*>(y, x)] > 0);
}

Node* Tree::min(Node* u, Node* v) {
	//ASSUMES u and v are comparable!
	if (LCA(u,v) == u)
		return v;
	else
		return u;
}

Node* Tree::LCA(std::set<Node*> V) {
	Node* lca = *(V.begin());
	for (Node* v : V) {
		lca = LCA(lca, v);
	}
	return lca;
}

Node* Tree::LCA(const string& ustr, const string& vstr) {
	if (V.count(ustr) == 0) {
		throw new app_exception("Tree::LCA(ustr, vstr); no vertex matching ustr");
	}
	if (V.count(vstr) == 0) {
		throw new app_exception("Tree::LCA(ustr, vstr); no vertex matching vstr");
	}
	return LCA(V[ustr], V[vstr]);
}

Node* Tree::LCA(Node* u, Node* v) {
	/**
	 * Not the fastest, but it should be fine for now.  Make it work BEFORE trying to make it faster!
	 */
	bool _debugging(false);
	assert(u != nullptr);
	assert(v != nullptr);
	DEBUG(cout << u->getLabel() << endl);
	DEBUG(cout << v->getLabel() << endl);
	DEBUG(cout << "LCA(" << u->getLabel() << ',' << v->getLabel() << ')' << endl);
	if (u == v) {
		return u;
	}
	std::set<Node*> ancestors;	// XXX could cache this
	Node* a = u;
	while (a != nullptr) {
		ancestors.insert(a);
		DEBUG(cout << "Adding node " << a->getLabel() << endl);
		a = a->getParent();
	}
	DEBUG(
			cout << "LCA finding ancestors of u=";
			cout << u->getLabel() << " = { ";
			for (Node* n : ancestors) {
				cout << n->getLabel() << " ";
			}
			cout << "}" << endl;
			);
	a = v;
	while (a != nullptr) {
		if (ancestors.find(a) != ancestors.end()) {
			DEBUG(cout << "LCA(" << u->getLabel() << ',' << v->getLabel() << ") = " << a->getLabel() << endl);
			return a;
		}
		a = a->getParent();
	}
	DEBUG(cout << "No LCA found!" << endl);
	throw new app_exception("Something has gone terribly wrong here: no LCA was found!");
}

Tree& Tree::operator=(const string& str) {
	throw new app_exception("Tree::operator=(const string& str) is not yet tested.");
	if (str[0] == '(') {
		*this = str;
		parsing::TokenList TL;
		istringstream instr(str);
		TL.tokenize(instr);
		parsing::SegDupParser p(TL);
		root = new Node("root");
		TL.reset();
		numVertices = 1;
		p.parseNewickSubtree(root, 'v');
		gatherVertices();
		for (auto n : V) {
			n.second->setTree(this);
		}
		calculateHeights(root);
		calcAncestry();
		_showInfo = false;
		return *this;
	} else {
		throw new app_exception("Tree copy constructor with a single string argument is expecting Newick format tree description.");
	}
}
ostream& operator<<(ostream& os, Tree& T) {
	if (T.getRoot() == nullptr) {
		return os;
	}
//	os << T.getLabel() << " = " << *(T.getRoot()) << std::endl;
	T.compressTraverseWrite(os);
	return os;
}

void Tree::putInternalVertices(std::set<Node*>& IV) {
	if (V.size() == 0) {
		gatherVertices();
	}
	for (auto v : V) {
		if (!v.second->isLeaf()) {
			IV.insert(v.second);
		}
	}
}

void Tree::setNodeLabel(const std::string& str, const std::string& newlabel) {
	Node* n(V[str]);
	setNodeLabel(n, newlabel);
}

void Tree::setNodeLabel(Node* n, const std::string& newlabel) {
	V.erase(n->getLabel());
	n->setLabel(newlabel);
	V[newlabel] = n;
}

} /* namespace segdup */
