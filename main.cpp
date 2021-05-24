#include <iostream>
#include "include/cmdline.h"
#include "operator.h"
#include "segment.h"
#include "utils.h"
#include "test.h"
#include "time.h"

using namespace std;

int main(int argc, char *argv[]) {
    cmdline::parser parser;
    parser.add<string>("action", 'a', "action", true, "",
            cmdline::oneof<string>("import", "average", "sigma", "corr", "ccorr", "acorr"));
    parser.add<string>("filename", 'f', "filename", false, "");
    parser.add<string>("output_name", 'o', "output file name", false, "");
    parser.add<double>("max_error", 'e', "max error of segmenting algorithm", false, 200);
    parser.add<string>("table", 't', "table name", false, "");
    parser.add<string>("table1", '\0', "table1 name", false, "");
    parser.add<string>("table2", '\0', "table2 name", false, "");
    parser.add<string>("data", 'd', "original data", false, "");
    parser.add<string>("data1", '\0', "original data1", false, "");
    parser.add<string>("data2", '\0', "original data2", false, "");
    parser.add<int>("shift", 's', "shift value", false, 0);
    parser.parse_check(argc, argv);
    string action = parser.get<string>("action");
    double error_guarantee = 0.;
    clock_t start,end;
    cout.precision(10);
    if(action == "import") {
        double max_error = max_error = parser.get<double>("max_error");
        string filename = parser.get<string>("filename");
        string output_name = parser.get<string>("output_name");
        vector<double> data = readfile(filename);
        SlidingWindow sw;
        std::vector<Segment> segments = sw.create_segments(data, max_error, linear_regression);
        write_segments(segments, output_name);
    }else if(action == "average") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        start = clock();
        double avg = average(segments, error_guarantee);
        end = clock();
        cout << fixed << "average: " << avg << " error guarantee: "
        << error_guarantee << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
        if(parser.exist("data")) {
            string data_name = parser.get<string>("data");
            vector<double> data = readfile(data_name);
            start = clock();
            double ta = true_average(data);
            end = clock();
            cout << fixed << "true average: " << ta << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
            cout << fixed << "true error: " << abs(ta - avg) << endl;
            cout << fixed << "ratio: " << abs((ta - avg) / ta) << endl;
        }
    }else if(action == "sigma") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        start = clock();
        double sigma = standard_deviation(segments, error_guarantee);
        end = clock();
        cout << fixed << "standard deviation: " << sigma << " error guarantee: " << error_guarantee
        << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
        if(parser.exist("data")) {
            string data_name = parser.get<string>("data");
            vector<double> data = readfile(data_name);
            start = clock();
            double ts = true_sigma(data);
            end = clock();
            cout << fixed << "true sigma: " << ts << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
            cout << fixed << "true error:" << abs(ts-sigma) << endl;
            cout << fixed << "ratio:" << abs((ts-sigma) / ts) << endl;
        }
    }else if(action == "corr") {
        string table1_name = parser.get<string>("table1");
        string table2_name = parser.get<string>("table2");
        vector<Segment> segments1 = read_segments(table1_name);
        vector<Segment> segments2 = read_segments(table2_name);
        start = clock();
        double corr = correlation(segments1, segments2, error_guarantee);
        end = clock();
        cout << fixed << "correlation: " << corr << " error guarantee: " << error_guarantee
        << " time cost: "<< (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
        if(parser.exist("data1") && parser.exist("data2")) {
            string data1_name = parser.get<string>("data1");
            string data2_name = parser.get<string>("data2");
            vector<double> data1 = readfile(data1_name);
            vector<double> data2 = readfile(data2_name);
            start = clock();
            double tc = true_corr(data1, data2);
            end = clock();
            cout << fixed << "true correlation: " << tc << " time cost: "<< (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
            cout << fixed << "true error: " << abs(tc - corr) << endl;
            cout << fixed << "ratio: " << abs((tc - corr) / tc) << endl;
        }
    }else if(action == "ccorr") {
        string table1_name = parser.get<string>("table1");
        string table2_name = parser.get<string>("table2");
        vector<Segment> segments1 = read_segments(table1_name);
        vector<Segment> segments2 = read_segments(table2_name);
        int shift = parser.get<int>("shift");
        start = clock();
        double ccorr = cross_correlation(segments1, segments2, shift, error_guarantee);
        end = clock();
        cout << fixed << "cross correlation: " << ccorr << " error guarantee: " << error_guarantee
        << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
        if(parser.exist("data1") && parser.exist("data2")) {
            string data1_name = parser.get<string>("data1");
            string data2_name = parser.get<string>("data2");
            vector<double> data1 = readfile(data1_name);
            vector<double> data2 = readfile(data2_name);
            start = clock();
            double tcc = true_ccorr(data1, data2, shift);
            end = clock();
            cout << fixed << "true cross-correlation: " << tcc << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
            cout << fixed << "true error: " << abs(tcc - ccorr) << endl;
            cout << fixed << "ratio: " << abs((tcc - ccorr) / tcc) << endl;
        }
    }else if(action == "acorr") {
        string table_name = parser.get<string>("table");
        vector<Segment> segments = read_segments(table_name);
        int shift = parser.get<int>("shift");
        start = clock();
        double acorr = auto_correlation(segments, shift, error_guarantee);
        end = clock();
        cout << fixed << "auto correlation: " << acorr << " error guarantee: " << error_guarantee
        << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
        if(parser.exist("data")) {
            string data_name = parser.get<string>("data");
            vector<double> data = readfile(data_name);
            start = clock();
            double ta = true_ccorr(data, data, shift);
            end = clock();
            cout << fixed << "true auto-correlation: " << ta << " time cost: "<< (end - start) / CLOCKS_PER_SEC << "s" << endl;
            cout << fixed << "true error: " << abs(ta - acorr) << endl;
            cout << fixed << "ratio: " << abs((ta - acorr) / ta) << endl;
        }
    }
    return 0;
}
