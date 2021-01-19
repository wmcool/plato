#include <iostream>
#include "cmdline.h"
#include "operator.h"
#include "segment.h"
#include "utils.h"
#include "test.h"

using namespace std;

int main(int argc, char *argv[]) {
    cmdline::parser parser;
    parser.add<string>("action", 'a', "action", true, "",
            cmdline::oneof<string>("import", "average", "sigma", "corr", "ccorr", "acorr"));
    parser.add<string>("filename", 'f', "filename", false, "");
    parser.add<string>("output_name", 'o', "output file name", false, "");
    parser.add<string>("table", 't', "table name", false, "");
    parser.add<string>("table1", '\0', "table1 name", false, "");
    parser.add<string>("table2", '\0', "table2 name", false, "");
    parser.add<int>("shift", 's', "shift value", false, 0);
    parser.parse_check(argc, argv);
    string action = parser.get<string>("action");

    double max_error = 200;
    double error_guarantee = 0.;
    cout.precision(10);
    if(action == "import") {
        string filename = parser.get<string>("filename");
        string output_name = parser.get<string>("output_name");
        vector<double> data = readfile(filename);
        SlidingWindow sw;
        std::vector<Segment> segments = sw.create_segments(data, max_error, linear_regression);
        write_segments(segments, output_name);
    }else if(action == "average") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        double avg = average(segments, error_guarantee);
        cout << fixed << "average: " << avg << " error guarantee: " << error_guarantee << endl;
        vector<double> data = readfile("data1");
        cout<< fixed << "true average: " << true_average(data) << endl;
    }else if(action == "sigma") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        double sigma = standard_deviation(segments, error_guarantee);
        cout << fixed << "standard deviation: " << sigma << " error guarantee: " << error_guarantee << endl;
//        vector<double> data = readfile("data");
//        cout<< fixed << "true sigma: " << true_sigma(data) << endl;
    }else if(action == "corr") {
        string table1_name = parser.get<string>("table1");
        string table2_name = parser.get<string>("table2");
        vector<Segment> segments1 = read_segments(table1_name);
        vector<Segment> segments2 = read_segments(table2_name);
        double corr = correlation(segments1, segments2, error_guarantee);
        cout << fixed << "correlation: " << corr << " error guarantee: " << error_guarantee << endl;
        vector<double> data1 = readfile("Open.txt");
        vector<double> data2 = readfile("Close.txt");
        cout<< fixed << "true correlation: " << true_corr(data1, data2) << endl;
    }else if(action == "ccorr") {
        string table1_name = parser.get<string>("table1");
        string table2_name = parser.get<string>("table2");
        vector<Segment> segments1 = read_segments(table1_name);
        vector<Segment> segments2 = read_segments(table2_name);
        int shift = parser.get<int>("shift");
        double ccorr = cross_correlation(segments1, segments2, shift, error_guarantee);
        cout << fixed << "cross correlation: " << ccorr << " error guarantee: " << error_guarantee << endl;
        vector<double> data1 = readfile("data1");
        vector<double> data2 = readfile("data2");
        cout<< fixed << "true cross-correlation: " << true_ccorr(data1, data2, shift) << endl;
    }else if(action == "acorr") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        int shift = parser.get<int>("shift");
        double acorr = auto_correlation(segments, shift, error_guarantee);
        cout << fixed << "auto correlation: " << acorr << " error guarantee: " << endl;
    }
    return 0;
}
