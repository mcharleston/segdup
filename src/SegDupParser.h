/*
 * NEXUSParser.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef SEGDUPPARSER_H_
#define SEGDUPPARSER_H_

#include "Node.h"
#include "parser.h"
#include "Tree.h"

namespace parsing {

const std::vector<std::string> NEXUSSuffixes = { ".nx", ".nxs", ".nex", ".nexus" };

class SegDupParser : public Parser {
private:
	segdup::Project* proj;

public:
	virtual ~SegDupParser() { }
	SegDupParser(const std::string& fileName, segdup::Project* pr);
	SegDupParser(TokenList& tl) : Parser(tl), proj(nullptr) {}

	void parse();
	void parseBranchLength(segdup::Node* v);
//	void parseDataBlock();
//	void parseDataFormat();
//	void parseDimensionsComponent();
	void parseSegDupBlock();
//	void parseMatrix();
	void parseNEXUSBlock();
	void parseNewickTree(segdup::Tree* T);
	void parseNewickSubtree(segdup::Node* v, char prefix);
//	void parseSetsBlock();
	void parseTaxaBlock();

	void skipBlock();
};

} /* namespace parsing */

#endif /* SEGDUPPARSER_H_ */
