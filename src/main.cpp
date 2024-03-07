// compile cmd:
// g++ parser.cpp main.cpp -o main.exe && main.exe

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <chrono>
#include <windows.h>

// constants.h define "NTHREAD", this value is used here and the .dll source code
#include "../include/constants.h"
#include "../include/parser.h"
#include "../include/ThreadPool.h"//https://github.com/progschj/ThreadPool

#define JJ 100
#define ZZ 800

using Clock = std::chrono::high_resolution_clock;

struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

dp calculate_array(int z, const std::vector<double> & seq_energy, const float * force, const float * energy) {
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    int X_DIM = 10;

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        temp_f.at(j) = force[j + z * X_DIM];
        temp_e.at(j) = energy[j + z * X_DIM];
        
        temp_e.at(j) += seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    double prob = 0;
    double Fprob = 0;
    double FFprob = 0;
    double Jprob = 0;
    double JJprob = 0;
    double temp1,temp2,temp3;
    for (int j = 0; j < seq_energy.size(); ++j) {

        temp1 = temp_e.at(j) - min_e;
        temp1 = temp1 > ENERGY_THRESHOLD ? 0.0 : std::exp(-temp1);
        temp2 = temp_f.at(j);

        prob += temp1;

        temp3 = temp2 * temp1;
        Fprob += temp3;
        FFprob += temp3 * temp2;

        temp3 = j * temp1;
        Jprob += temp3;
        JJprob += j * temp3;

    }

    dp point;
    point.extension_total = z;
    point.force_average = Fprob/prob;
    point.extension_DNA = z - point.force_average/PILLARSTIFFNESS;
    point.force_SD = std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob));
    point.junzipped_average = Jprob/prob;
    point.junzipped_SD = std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob));
    return point;

}

int main() {

    constexpr int max_j_unzipped = 8000;
    int X_DIM = (max_j_unzipped + NTHREAD - 1);
    X_DIM -= X_DIM % NTHREAD;
    int Y_DIM = (int)((max_j_unzipped + ARMLENGTH) * L0SS * 2) + NTHREAD - 1;
    Y_DIM -= Y_DIM % NTHREAD;
    
    float *force = new float [X_DIM * Y_DIM];
    float *energy = new float [X_DIM * Y_DIM];

    {
        HINSTANCE hDLL = LoadLibrary("cuLUT.dll"); // cast to "LPCWSTR" to suppress the error
        if (hDLL == NULL) {
            std::cerr << "cuUNZIP.dll open failed.";
            std::cerr << "Error code:" << GetLastError();
            return 1;
        }

        //void get_LUT(int X_DIM, int Y_DIM, float * force_out, float * energy_out)
        auto get_LUT = (void (*)(int, int, float*, float*))GetProcAddress(hDLL, "get_LUT");
        if (get_LUT == NULL) {
            std::cerr << "get_LUT() doesn't exist.";
            std::cerr << "Error code:" << GetLastError();
            return 1;
        }

        get_LUT(X_DIM, Y_DIM, force, energy);
        FreeLibrary(hDLL);
        printf("force[%d][%d] = %f\n", JJ,ZZ, force[JJ + ZZ * X_DIM]);
        printf("energy[%d][%d] = %f\n", JJ,ZZ, energy[JJ + ZZ * X_DIM]);
    }

    std::string path {"../data/NEB_H5alpha_Accessory_colonization_factor_AcfD.txt"};
    std::string seq = parser::readtxt_firstline(path);
    std::vector<double> seq_energy = parser::calculate_sequence_energy(seq);
    //test, the "43048.29409015185" is from my python code
    std::cout << "Sequence energy differs from python program by "<< KT * (*(seq_energy.cend() - 1)) - 43048.29409015185 << std::endl;// no difference, good

    auto start = Clock::now();

    int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads: " << numThreads << std::endl;
    ThreadPool pool(numThreads);

    std::vector<std::future<dp>> results;

    for(int extension = 1; extension < static_cast<int>(1.2 * seq_energy.size()); ++extension) {
        results.emplace_back(
            pool.enqueue([extension, &seq_energy, force, energy]{
                return calculate_array(extension, seq_energy, force, energy);
            })
        );
    }

    std::vector<dp> result_array;
    for(auto && result: results)
        result_array.emplace_back(result.get());

    std::ofstream fout("res.csv");
    if (!fout) {
        std::cerr << "Error opening file. Cannot save result" << std::endl;
        return 1;
    }

    fout << "total extension (nm),DNA extension (nm),average force (pN),sd force (pN),average bp unzipped,sd bp unzipped" << std::endl;

    char delimit = ',';
    std::for_each(result_array.cbegin(), result_array.cend(), [&fout, delimit](const dp & point) {
        fout << point.extension_total << delimit;
        fout << point.extension_DNA << delimit;
        fout << point.force_average << delimit;
        fout << point.force_SD << delimit;
        fout << point.junzipped_average << delimit;
        fout << point.junzipped_SD << std::endl;;
    });

    // Close the file
    fout.close();
    std::cout << "Result has been written to res.csv" << "'." <<std::endl;
    
    auto end = Clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;



    delete[] force, energy;
    return 0;
}