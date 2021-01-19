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
    std::ofstream out(file_name, std::ios::out | std::ios::binary);
    for(Segment s : segments) {
        out.write((char*)&s.begin, sizeof(int));
        out.write((char*)&s.end, sizeof(int));
        out.write((char*)&s.epsilon, sizeof(double));
        out.write((char*)&s.f, sizeof(double));
        out.write((char*)&s.gamma, sizeof(double));
        for(double param : s.param) {
            out.write((char*)&param, sizeof(double));
        }
    }
    out.close();
}

std::vector<Segment> read_segments(const std::string& file_name) {
    std::ifstream in(file_name, std::ios::in | std::ios::binary);
    std::vector<Segment> segments;
    in.seekg(0, std::ios::end);
    int size = in.tellg();
    in.seekg(0, std::ios::beg);
    while(in.tellg() < size) {
        int begin;
        int end;
        double epsilon;
        double f;
        double gamma;
        std::vector<double> param;
        double a;
        double b;
        in.read((char*)(&begin), sizeof(int));
        in.read((char*)(&end), sizeof(int));
        in.read((char*)(&epsilon), sizeof(double));
        in.read((char*)(&f), sizeof(double));
        in.read((char*)(&gamma), sizeof(double));
        in.read((char*)(&a), sizeof(double));
        in.read((char*)(&b), sizeof(double));
        param.push_back((a));
        param.push_back((b));
        Segment segment(param, begin, end, epsilon, f, gamma);
        segments.push_back(segment);
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

