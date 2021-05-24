//
// Created by wmcool on 2020/12/23.
//
#include <vector>
#include <string>
#ifndef PLATO_SEGMENT_H
#define PLATO_SEGMENT_H

class Segment{
public:
    std::vector<double> param;
    int begin;
    int end;
    double epsilon;
    double f;
    double gamma;
    Segment(std::vector<double> param, int begin, int end, double epsilon, double f, double gamma);
};

class Estimator{
public:
    virtual std::vector<Segment> create_segments(std::vector<double> data, double max_error, double (*ff)(std::vector<double> data, std::vector<double> &param, int offset, double& epsilon, double& f, double& gamma)) = 0;
};

class SlidingWindow : public Estimator{
public:
    std::vector<Segment> create_segments(std::vector<double> data, double max_error, double (*ff)(std::vector<double> data, std::vector<double> &param, int offset, double& epsilon, double& f, double& gamma)) override;

//    std::vector<Segment> create_segments_opt(std::vector<double> data, double max_error);
};

std::string serialize(const Segment& segment);

Segment deserialize(const std::string& s);

double linear_regression(std::vector<double> data, std::vector<double> &param, int offset, double& epsilon, double& f, double& gamma);

// operator of segments
std::vector<Segment> constant(double value, int begin, int end, double error_guarantee);

std::vector<Segment> shift(const std::vector<Segment>& segments, int m);

double sum_of_add(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee);

double sum_of_subtract(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee);

double sum_of_times(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee);

double compute_gamma_f(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2);

double compute_gamma_gamma(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2);

#endif //PLATO_SEGMENT_H


