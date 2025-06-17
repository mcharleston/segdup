/*
 * SegDup.h
 *
 *  Created on: 14 Dec 2023
 *      Author: mac
 */

#ifndef SRC_SEGDUP_H_
#define SRC_SEGDUP_H_

#include <map>
#include <string>
#include <vector>

#include "CophyMultiMap.h"
#include "DupMove.h"
#include "EventCount.h"

namespace segdup {

double CSD(const EventCount& ec);

}

void SelectNextConfiguration(segdup::CophyMultiMap& CMM, double T, std::vector<segdup::DupMove*>& moves, std::vector<double> probs);
void Algorithm2(segdup::CophyMultiMap& CMM, std::vector<segdup::DupMove*> moves, std::vector<double> probs, std::map<std::string, int>* sampledDistribution = nullptr);


#endif /* SRC_SEGDUP_H_ */
