#pragma once

#include <string>
#include <vector>
//a struct to save data point by point
struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

//Calculate DNA sequence's energy
namespace util {
    std::string creat_path_out(const std::string & path_in);
    std::string readtxt_firstline(const std::string & path);
    std::vector<double> calculate_sequence_energy(const std::string & sequence);
}

