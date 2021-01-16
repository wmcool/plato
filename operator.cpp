//
// Created by wmcool on 2021/1/2.
//
#include <vector>
#include <iostream>
#include <cmath>
#include "operator.h"
#include "segment.h"

double average(const std::vector<Segment>& segments, double& error_guarantee) {
    double res = 0;
    error_guarantee = 0;
    int n = segments[segments.size()-1].end - segments[0].begin;
    for(Segment segment : segments) {
        for(int i=segment.begin;i<segment.end;i++) {
            res += segment.param[0] * i + segment.param[1];
        }
        error_guarantee += segment.gamma;
    }
    return res/n;
}

double standard_deviation(std::vector<Segment>& segments, double& error_guarantee) {
    double mu_error_guarantee = 0;
    double mu = average(segments, mu_error_guarantee);
    int n = segments[segments.size()-1].end - segments[0].begin;
    for(Segment segment : segments) {
        segment.param[1] -= mu;
    }
    double res = sum_of_times(segments, segments, error_guarantee);
    error_guarantee /= n;
    res = res / n;
    error_guarantee = std::sqrt(res + error_guarantee) - std::sqrt(res);
    res = std::sqrt(res);
    return res;
}

double correlation(std::vector<Segment>& segments1, std::vector<Segment>& segments2, double& error_guarantee) {
    double mu1_error_guarantee = 0;
    double mu2_error_guarantee = 0;
    double mu1 = average(segments1, mu1_error_guarantee);
    double mu2 = average(segments2, mu2_error_guarantee);
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    int n = end - begin;
//    std::cout << "mu1: " << mu1 << " mu2: " << mu2 << std::endl;
    for(int i=0;i<segments1.size();i++) {
        segments1[i].param[1] -= mu1;
    }
    for(int j=0;j<segments2.size();j++) {
        segments2[j].param[1] -= mu2;
    }
    double res = sum_of_times(segments1, segments2, error_guarantee);
    res /= n;
    error_guarantee /= n;
    double var1_error_guarantee = 0;
    double var2_error_guarantee = 0;
    double var1 = standard_deviation(segments1, var1_error_guarantee);
    double var2 = standard_deviation(segments2, var2_error_guarantee);
    std::cout << "error_guarantee: " << error_guarantee << " var1_error_guarantee: " << var1_error_guarantee
    << " var2_error_guarantee: " << var2_error_guarantee << std::endl;
    std::cout << "sigma1: " << var1 << " sigma2: " << var2 << std::endl;
    error_guarantee = (error_guarantee * var1 + res * var1_error_guarantee) / ((var1 - var1_error_guarantee) * var1);
    res /= var1;
    error_guarantee = (error_guarantee * var2 + res * var2_error_guarantee) / ((var2 - var2_error_guarantee) * var2);
    res /= var2;
    return res;
}

double cross_correlation(std::vector<Segment>& segments1, std::vector<Segment>& segments2, int m, double& error_guarantee) {
    std::vector<Segment> shift_segments2 = shift(segments2, m);
    return correlation(segments1, shift_segments2, error_guarantee);
}

double auto_correlation(std::vector<Segment>& segments1, int m, double& error_guarantee) {
    std::vector<Segment> shift_segments1 = shift(segments1, m);
    return correlation(segments1, shift_segments1, error_guarantee);
}