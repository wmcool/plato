//
// Created by wmcool on 2021/1/1.
//

#ifndef PLATO_OPERATOR_H
#define PLATO_OPERATOR_H
#include <vector>
#include "segment.h"

double average(const std::vector<Segment>& segments, double& error_guarantee);

double standard_deviation(std::vector<Segment>& segments, double& error_guarantee);

double correlation(std::vector<Segment>& segments1, std::vector<Segment>& segments2, double& error_guarantee);

double cross_correlation(std::vector<Segment>& segments1, std::vector<Segment>& segments2, int shift, double& error_guarantee);

double auto_correlation(std::vector<Segment>& segments1, int shift, double& error_guarantee);

#endif //PLATO_OPERATOR_H
