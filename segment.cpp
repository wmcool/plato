//
// Created by wmcool on 2020/12/24.
//
#include "segment.h"
#include <cmath>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <assert.h>

Segment::Segment(std::vector<double> param, int begin, int end, double epsilon, double f, double gamma) :
param(std::move(param)), begin(begin), end(end), epsilon(epsilon), f(f), gamma(gamma){}

Interval::Interval(int begin, int end) :
begin(begin), end(end){}

std::string serialize(const Segment& segment) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6) << segment.begin << " " << segment.end << " " <<
    segment.epsilon << " " << segment.f << " " << segment.gamma << " ";
    int n = segment.param.size();
    for(int i=0;i<n;i++) {
        ss << segment.param[i] << " ";
    }
    return ss.str();
}

Segment deserialize(const std::string& s) {
    std::stringstream ss(s);
    int begin;
    int end;
    double epsilon;
    double f;
    double gamma;
    std::vector<double> param;
    double cur_param;
    ss >> begin >> end >> epsilon >> f >> gamma;
    while(ss >> cur_param) {
        param.push_back(cur_param);
    }
//    std::cout << std::fixed << std::setprecision(15)<< "f: " << f << std::endl;
    return Segment(param, begin, end, epsilon, f, gamma);
}

double linear_regression(std::vector<double> data, std::vector<double> &param, int offset, double& epsilon, double& f, double& gamma) {

    // compute param
    int n = data.size();
    double x_bar = offset + (double)(n-1)/2;
    double y_bar = 0;
    double x_2_sum = 0;
    double y_sum = 0;
    double xy = 0;
    for(int i=0;i<n;i++) {
        xy += (i+offset) * data[i];
        x_2_sum += (double)(i+offset)*(i+offset);
        y_sum += data[i];
    }
    y_bar = y_sum / n;
    double a = (xy - n * x_bar * y_bar) / (x_2_sum - n * x_bar * x_bar);
    double b = y_bar - a * x_bar;
    param.clear();
    param.push_back(a);
    param.push_back(b);

    // compute error
    epsilon = 0;
    f = 0;
    gamma = 0;
    for(int i=0;i<n;i++) {
        epsilon += (data[i] - a*(i+offset) - b) * (data[i] - a*(i+offset) - b);
        f += (a*(i+offset) + b) * (a*(i+offset) + b);
        gamma += (a*(i+offset) + b);
    }
    epsilon = std::sqrt(epsilon);
    f = std::sqrt(f);
    gamma = std::abs(y_sum - gamma);

    // return error (default mse)
    return epsilon;
}

std::vector<Segment> SlidingWindow::create_segments(std::vector<double> data, double max_error, double (*ff)(std::vector<double> data, std::vector<double> &param, int offset, double& epsilon, double& f, double& gamma)) {
    std::vector<Segment> segments;
    int anchor = 0;
    bool finished = false;
    int n = data.size();
    while(!finished) {
        int i = 2;
        std::vector<double> sub_data(data.begin() + anchor, data.begin() + anchor + i);
        double epsilon = 0;
        double f = 0;
        double gamma = 0;
        std::vector<double> param;
        double pre_epsilon = epsilon;
        double pre_f = f;
        double pre_gamma = gamma;
        std::vector<double> pre_param = param;
        while(ff(sub_data, param, anchor, epsilon, f, gamma) < max_error && anchor + i <= n) {
            pre_epsilon = epsilon;
            pre_f = f;
            pre_gamma = gamma;
            pre_param = param;
            sub_data.push_back(data[anchor+i]);
            if((anchor+i) == (n-1)) i++;
            i++;
        }
        segments.emplace_back(pre_param, anchor, anchor + i - 1, pre_epsilon, pre_f, pre_gamma);
        anchor += (i-1);
        if(anchor >= n) finished = true;
    }
    return segments;
}

std::vector<Segment> SlidingWindow::create_segments_opt(std::vector<double> data, double max_error) {
    std::vector<Segment> segments;
    int anchor = 0;
    int n = data.size();
    while(true) {
        double x_sum = anchor + anchor + 1;
        double x_2_sum = anchor*anchor + (anchor+1)*(anchor+1);
        double y_sum = data[anchor] + data[anchor+1];
        double xy_sum = anchor * data[anchor] + (anchor+1) * data[anchor+1];
        int i = 2;
        double x_bar = x_sum / 2;
        double y_bar = y_sum / 2;
        double a = (xy_sum - i * x_bar * y_bar) / (x_2_sum - i * x_bar * x_bar);
        double b = y_bar - a * x_bar;
        double epsilon = 0.;
        double f = 0.;
        double gamma = 0;
        double pre_a;
        double pre_b;
        while(anchor + i < n && epsilon <= max_error) {
            x_sum += anchor+i;
            y_sum += data[anchor+i];
            x_2_sum += (anchor+i) * (anchor+i);
            xy_sum += (anchor+i) * data[anchor+i];
            x_bar = anchor + (double)i / 2;
            y_bar = y_sum / (i+1);
            pre_a = a;
            pre_b = b;
            a = (xy_sum - (i+1) * x_bar * y_bar) / (x_2_sum - (i+1) * x_bar * x_bar);
            b = y_bar - a * x_bar;
            if(anchor + i == (n-1)) {
                i++;
                break;
            }
            epsilon += (data[anchor+i] - (a*(anchor+i) + b)) * (data[anchor+i] - (a*(anchor+i) + b));
            if(epsilon > max_error) break;
            i++;
        }
        epsilon = 0.;
        a = pre_a;
        b = pre_b;
        for(int k=anchor;k<anchor+i;k++) {
            epsilon += (data[k] - (a*k+b)) * (data[k] - (a*k+b));
            f += (a*k+b) * (a*k+b);
            gamma += (a*k+b);
        }
        epsilon = sqrt(epsilon);
        f = sqrt(f);
        if(i == n) gamma = std::abs(y_sum - gamma);
        else gamma = std::abs(y_sum - data[anchor+i] - gamma);
        std::vector<double> param;
        param.push_back(a);
        param.push_back(b);
        segments.emplace_back(param, anchor, anchor + i, epsilon, f, gamma);
        anchor += i;
        if(anchor >= n) break;
    }
    return segments;
}

std::vector<Segment> constant(double value, int begin, int end, double error_guarantee) {
    std::vector<double> param;
    param.push_back(0);
    param.push_back(value);
    Segment only(param, begin, end, 0, std::sqrt(end - begin) * value, (end - begin) * error_guarantee);
    std::vector<Segment> segments;
    segments.push_back(only);
    return segments;
}

std::vector<Segment> shift(const std::vector<Segment>& segments, int m) {
    std::vector<Segment> res;
    for(Segment segment : segments) {
        std::vector<double> param;
        param.push_back(segment.param[0]);
        param.push_back(segment.param[1] - param[0] * m);
        Segment to_add(param, segment.begin + m, segment.end + m, segment.epsilon, segment.f, segment.gamma);
        res.push_back(to_add);
    }
    return res;
}

double sum_of_add(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    double res = 0;
    error_guarantee = 0;
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    for(int k=begin;k<end;k++) {
        double value1 = segments1[i].param[0] * k + segments1[i].param[1];
        double value2 = segments2[j].param[0] * k + segments2[j].param[1];
        res += value1 + value2;
        if(k == segments1[i].end-1 || k == end-1) {
            error_guarantee += segments1[i].gamma;
            i++;
        }
        if(k == segments2[j].end-1 || k == end-1) {
            error_guarantee += segments2[j].gamma;
            j++;
        }
    }
    return res;
}

double sum_of_subtract(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    double res = 0;
    error_guarantee = 0;
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    for(int k=begin;k<end;k++) {
        double value1 = segments1[i].param[0] * k + segments1[i].param[1];
        double value2 = segments2[j].param[0] * k + segments2[j].param[1];
        res += value1 - value2;
        if(k == segments1[i].end-1 || k == end-1) {
            error_guarantee += segments1[i].gamma;
            i++;
        }
        if(k == segments2[j].end-1 || k == end-1) {
            error_guarantee += segments2[j].gamma;
            j++;
        }
    }
    return res;
}

double sum_of_times(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2, double& error_guarantee) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    double res = 0;
    error_guarantee = 0;
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    for(int k=begin;k<end;k++) {
        double value1 = segments1[i].param[0] * k + segments1[i].param[1];
        double value2 = segments2[j].param[0] * k + segments2[j].param[1];
        res += value1 * value2;
        if(k == segments1[i].end-1 || k == end-1) i++;
        if(k == segments2[j].end-1 || k == end-1) j++;
    }
    error_guarantee += compute_gamma_f_opt(segments1, segments2);
    error_guarantee += compute_gamma_f_opt(segments2, segments1);
//    error_guarantee += compute_gamma_gamma_opt(segments1, segments2);
    return res;
}

double compute_gamma_f_opt(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    double res = 0;
//    double sos = 0; // sum of square
    std::vector<double> fuck;
    bool overlap = false;
    for(int k=begin;k<end;k++) {
//        double value1 = segments1[i].param[0] * k + segments1[i].param[1];
        double value2 = segments2[j].param[0] * k + segments2[j].param[1];
        fuck.push_back(value2);
//        sos += (value2 - value1) * (value2 - value1);
//        sos += value2 * value2;
        if(k == segments2[j].end-1 || k == end-1) {
            if(k != segments1[i].end - 1 && k != end-1) overlap = true;
            j++;
        }
        if(k == segments1[i].end-1 || k == end-1) {
//            res += segments1[i].epsilon * std::sqrt(sos);
//            sos = 0;
//            i++;
            if(overlap) {
                std::vector<double> param;
                double epsilon;
                double f;
                double gamma;
                double sos = linear_regression(fuck, param, k, epsilon, f, gamma);
                res += segments1[i].epsilon * sos;
            }
            overlap = false;
            i++;
            fuck.clear();
        }
    }
    return res;
}

double compute_gamma_f(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    double res = 0;
    double sos = 0; // sum of square
//    std::vector<double> fuck;
    bool overlap = false;
    for(int k=begin;k<end;k++) {
        double value1 = segments1[i].param[0] * k + segments1[i].param[1];
        double value2 = segments2[j].param[0] * k + segments2[j].param[1];
//        fuck.push_back(value2);
        sos += (value2 - value1) * (value2 - value1);
//        sos += value2 * value2;
        if(k == segments2[j].end-1 || k == end-1) {
            if(k != segments1[i].end - 1 && k != end-1) overlap = true;
            j++;
        }
        if(k == segments1[i].end-1 || k == end-1) {
            res += segments1[i].epsilon * std::sqrt(sos);
            sos = 0;
            i++;
//            if(overlap) {
//                std::vector<double> param;
//                double epsilon;
//                double f;
//                double gamma;
//                double sos = linear_regression(fuck, param, k, epsilon, f, gamma);
//                res += segments1[i].epsilon * sos;
//            }
//            overlap = false;
//            i++;
//            fuck.clear();
        }
    }
    return res;
}

double compute_gamma_gamma(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    double e1 = 0;
    double e2 = 0;
    while(segments1[i].end < end) {
        e1 += segments1[i].epsilon * segments1[i].epsilon;
        i++;
    }
    e1 += segments1[i].epsilon * segments1[i].epsilon;
    e1 = std::sqrt(e1);
    while(segments2[j].end < end) {
        e2 += segments2[j].epsilon * segments2[j].epsilon;
        j++;
    }
    e2 += segments2[j].epsilon * segments2[j].epsilon;
    e2 = std::sqrt(e2);
    return e1 * e2;
}

double compute_gamma_gamma_opt(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2) {
    int begin = std::max(segments1[0].begin, segments2[0].begin);
    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
    int i = 0;
    while(segments1[i].end <= begin) i++;
    int j = 0;
    while(segments2[j].end <= begin) j++;
    int local_start = begin;
    std::vector<Interval> opt;
    while(segments1[i].end < end && segments2[j].end < end) {
        int local_end = std::max(segments1[i].end, segments2[j].end);
        opt.emplace_back(local_start, local_end);
        local_start = local_end;
        while(segments1[i].end < end && segments1[i].end <= local_end) i++;
        while(segments2[j].end < end && segments2[j].end <= local_end) j++;
    }
    opt.emplace_back(local_start, end);
    double e = 0.;
    i = 0;
    while(segments1[i].end <= begin) i++;
    j = 0;
    while(segments2[j].end <= begin) j++;
    for(int k=0;k<opt.size();k++) {
        double e1 = 0.;
        while(i < segments1.size() && segments1[i].begin < opt[k].end) {
            e1 += segments1[i].epsilon * segments1[i].epsilon;
            if(segments1[i].end > opt[k].end) break;
            i++;
        }
        e1 = sqrt(e1);
        double e2 = 0;
        while(j < segments2.size() && segments2[j].begin < opt[k].end) {
            e2 += segments2[j].epsilon * segments2[j].epsilon;
            if(segments2[j].end > opt[k].end) break;
            j++;
        }
        e2 = sqrt(e2);
        e += e1 * e2;
    }
    return e;
}

//double compute_gamma_gamma(const std::vector<Segment>& segments1, const std::vector<Segment>& segments2) {
//    std::vector<int> opt;
//    int begin = std::max(segments1[0].begin, segments2[0].begin);
//    int end = std::min(segments1[segments1.size()-1].end, segments2[segments2.size()-1].end);
//    int i = 0;
//    while(segments1[i].end <= begin) i++;
//    int j = 0;
//    while(segments2[j].end <= begin) j++;
//    int anchor = begin;
//    while(anchor < end) {
//        if(segments1[i].end > end || segments2[j].end > end) anchor = end;
//        else anchor = std::max(segments1[i].end, segments2[j].end);
//        opt.push_back(anchor);
//        while(i < segments1.size() && segments1[i].end <= anchor) i++;
//        while(j < segments2.size() && segments2[j].end <= anchor) j++;
//    }
//    double e1 = 0;
//    double e2 = 0;
//    i = 0;
//    j = 0;
//    double res = 0;
//    for(int k=0;k<opt.size();k++) {
//        e1 = 0;
//        e2 = 0;
//        anchor = opt[k];
//        while(i < segments1.size() && segments1[i].begin < anchor) {
//            e1 += segments1[i].epsilon * segments1[i].epsilon;
//            if(segments1[i].end > anchor) break;
//            i++;
//        }
//        while(j < segments2.size() && segments2[j].begin < anchor) {
//            e2 += segments2[j].epsilon * segments2[j].epsilon;
//            if(segments2[j].end > anchor) break;
//            j++;
//        }
//        res += std::sqrt(e1) * std::sqrt(e2);
//    }
//    return res;
//}

//int main() {
//    std::vector<double> fuck;
//    fuck.push_back(78.44);
//    fuck.push_back(79.34);
//    std::vector<double> param;
//    double epsilon;
//    double f;
//    double gamma;
//    std::cout << linear_regression(fuck, param, 52423, epsilon, f, gamma) << std::endl;
//    return 0;
//}
