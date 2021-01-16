//
// Created by wmcool on 2021/1/7.
//

#ifndef PLATO_TEST_H
#define PLATO_TEST_H

#include <vector>

double true_average(const std::vector<double>& data);

double true_sigma(const std::vector<double>& data);

double true_corr(const std::vector<double>& data1, const std::vector<double>& data2);

double true_ccorr(const std::vector<double>& data1, const std::vector<double>& data2, int m);

#endif //PLATO_TEST_H
