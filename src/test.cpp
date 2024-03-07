// command to build:
// g++ test.cpp -o test.exe -L. -l_parser

#include <iostream>
#include <string>
#include <vector>
#include "../include/constants.h"
#include "../include/temp_data.h"
#include "../include/parser.h"

int main() {

    //std::string seq {temp_data::data};
    std::string path {"../test_data/NEB_H5alpha_Accessory_colonization_factor_AcfD.txt"};
    std::string seq = parser::readtxt_firstline(path);
    std::vector<double> seq_energy = parser::calculate_sequence_energy(seq);
    //test, the "43048.29409015185" is from my python code
    std::cout << "Sequence energy differs from python program by "<< KT * (*(seq_energy.cend() - 1)) - 43048.29409015185 << std::endl;// no difference, good
    

    return 0;
}