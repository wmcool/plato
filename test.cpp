//
// Created by wmcool on 2021/1/7.
//

#include "test.h"
#include <cmath>
#include <iostream>

double true_average(const std::vector<double>& data) {
    double res = 0;
    int n = data.size();
    for(int i=0;i<n;i++) {
        res += data[i];
    }
    res /= n;
    return res;
}

double true_sigma(const std::vector<double>& data) {
    double res = 0;
    double mu = true_average(data);
    int n = data.size();
    for(int i=0;i<n;i++) {
        res += (data[i] - mu) * (data[i] - mu);
    }
    res = std::sqrt(res / n);
    return res;
}

double true_corr(const std::vector<double>& data1, const std::vector<double>& data2) {
    double mu1 = true_average(data1);
    double mu2 = true_average(data2);
    double sigma1 = true_sigma(data1);
    double sigma2 = true_sigma(data2);
    double sum = 0;
    int n = data1.size();
    for(int i=0;i<n;i++) {
        sum += (data1[i] - mu1) * (data2[i] - mu2);
    }
    sum /= n;
    return sum / (sigma1 * sigma2);
}

double true_ccorr(const std::vector<double>& data1, const std::vector<double>& data2, int m) {
    int begin = m;
    int end = data1.size();
    double mu1 = true_average(data1);
    double mu2 = true_average(data2);
    double sigma1 = true_sigma(data1);
    double sigma2 = true_sigma(data2);
    double sum = 0;
    for(int k=begin;k<end;k++) {
        sum += (data1[k] - mu1) * (data2[k-m] - mu2);
    }
    sum /= (end - begin);
    return sum / (sigma1 * sigma2);
}

