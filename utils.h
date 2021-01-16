//
// Created by wmcool on 2021/1/2.
//

#ifndef PLATO_UTILS_H
#define PLATO_UTILS_H

#include <vector>
#include <string>

std::vector<double> readfile(const std::string& file_name);

void write_segments(const std::vector<Segment>& segments, const std::string& file_name);

std::vector<Segment> read_segments(const std::string& file_name);

#endif //PLATO_UTILS_H
