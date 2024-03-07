#pragma once

#include <string>
#include <vector>

//Calculate DNA sequence's energy
namespace parser {
    std::string creat_path_out(const std::string & path_in);
    std::string readtxt_firstline(const std::string & path);
    std::vector<double> calculate_sequence_energy(const std::string & sequence);
}

