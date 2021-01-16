//
// Created by wmcool on 2020/12/30.
//

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "segment.h"
#include "utils.h"

std::vector<double> readfile(const std::string& file_name) {
    std::ifstream in(file_name);
    std::vector<double> data;
    double cur;
    double avg = 0;
    int n = 0;
    while(in >> cur) {
        data.push_back(cur);
        avg += cur;
        n++;
    }
    in.close();
//    std::cout<< std::fixed << std::setprecision(10) << avg/n << std::endl;
    return data;
}

void write_segments(const std::vector<Segment>& segments, const std::string& file_name) {
    std::ofstream out(file_name);
    int n = segments.size();
    for(int i=0;i<n;i++) {
        out << serialize(segments[i]) << "\n";
    }
    out.close();
}

std::vector<Segment> read_segments(const std::string& file_name) {
    std::ifstream in(file_name);
    std::string s;
    std::vector<Segment> segments;
    while(std::getline(in, s)) {
        segments.push_back(deserialize(s));
    }
    in.close();
    return segments;
}

//int main() {
////    std::vector<Segment> segments = read_segments("temperature.segs");
////    for(int i=0;i<segments.size();i++) {
////        std::cout<<"segment"<< i << " a: " << segments[i].param[0] << " b: " << segments[i].param[1] << std::endl;
////        std::cout<<"scope: [" << segments[i].begin << ", " << segments[i].end-1 << "]"<<std::endl;
////    }
//    std::vector<double> data = readfile("data");
//    return 0;
//}

