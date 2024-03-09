#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <windows.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <chrono>
#include <filesystem>

using Clock = std::chrono::high_resolution_clock;

// constants.h define "NTHREAD", this value is used here and the .dll source code
#include "../include/constants.h"
/*
// There's no need to polish further, I am brain bleeding now.
// Another way is using concurrent kernels, may do it later?
*/

#ifndef NOTEST
#define JJ 100
#define ZZ 800
// f should be 16.3295 at 1600
#define EXT1 1600
// f should be 13.8639 at 4402
#define EXT2 2139
#define EXT3 2140
#endif

// =============================================util functions==========================================

// single read. todo: read the fasta directly, multiple line read.
std::string readtxt_firstline(const std::string & path) {
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return "0";
    }

    std::string line;
    if (getline(file, line)) {
        // std::cout << "Trunk sequence '" + path + "' read successfully." << std::endl;
    } else {
        std::cerr << "Failed to read a line from the file." << std::endl;
    }

    file.close();

    return line;
}
// generate a file path for output result from the input file path.
std::string creat_path_out(const std::string & path_in) {
    size_t lastSlashIndex = path_in.rfind('/');

    std::string parentPath;
    if (lastSlashIndex != std::string::npos) {
        parentPath = path_in.substr(0, lastSlashIndex) + "/";
    }

    std::string fullname = path_in.substr(lastSlashIndex + 1, std::string::npos);

    return parentPath + "out_" + fullname.substr(0, fullname.rfind('.')) + ".csv";;
}
// sequence's bp energy
std::vector<double> calculate_sequence_energy(const std::string & sequence) {
    //from the DNA sequence, calculate the energy at every j_unzipped.

    std::vector<double> sequenceEnergy;
    if (sequence.size() < 2) {
        std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
        return sequenceEnergy;
    }

    double accum = 0.0;
    *std::back_inserter(sequenceEnergy) = 0.0;
    std::transform(sequence.cbegin() + 1, sequence.cend(), std::back_inserter(sequenceEnergy), 
    [&accum](const char& bp2) {
        accum += BPEnergy::lookup_bp_energy(*(&bp2 - 1), bp2);
        return accum;
        }
    );//equivalent to std::partial_sum()

    return sequenceEnergy;
}

// main
int main(int argc, char** argv) {
    if (argc == 1) {
        return 0;
    }
    const std::string folder_path {argv[1]};

    // LUT's size along j and z directions, and:
    //      j = threadIdx.x + blockIdx.x * blockDim.x
    //      z = threadIdx.y + blockIdx.y * blockDim.y
    constexpr int LUT_j_dim = 16384;
    constexpr int LUT_z_dim = 16384;
    float *force_LUT;
    float *energy_LUT;
    force_LUT = new float [LUT_j_dim * LUT_z_dim];
    energy_LUT = new float [LUT_j_dim * LUT_z_dim];

    // this scope is common, will not be a thread
    {
        HINSTANCE hDLL = LoadLibrary("cuLUT.dll"); 
        if (hDLL == NULL) {
            std::cerr << "cuLUT.dll open failed.";
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

        get_LUT(LUT_j_dim, LUT_z_dim, force_LUT, energy_LUT);
        FreeLibrary(hDLL);
    }

    float *d_force_LUT;
    float *d_energy_LUT;
    cudaMalloc((void**) &d_force_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim);
    cudaMalloc((void**) &d_energy_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim);
    cudaMemcpy(d_force_LUT, force_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim, cudaMemcpyHostToDevice);
    cudaMemcpy(d_energy_LUT, energy_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim, cudaMemcpyHostToDevice);
    const LUT lut {LUT_j_dim, LUT_z_dim, d_force_LUT, d_energy_LUT};

    HINSTANCE hDLL = LoadLibrary("cuUnzip.dll"); 
    if (hDLL == NULL) {
        std::cerr << "cuUnzip.dll open failed.";
        std::cerr << "Error code:" << GetLastError();
        exit(0);
    }
    //data unzip(int seq_len, const float * d_seq_energy, LUT lut)
    auto unzip = (data (*)(int, const float * , LUT))GetProcAddress(hDLL, "unzip");
    if (unzip == NULL) {
        std::cerr << "unzip() doesn't exist.";
        std::cerr << "Error code:" << GetLastError();
        exit(0);
    }

    std::string fullFilePath {};
    for (const auto& entry : (std::filesystem::directory_iterator(folder_path))) {
        auto start = Clock::now();

        fullFilePath = folder_path + "/" + entry.path().filename().string();
        if (!entry.is_regular_file()) {
            std::cerr << "Error: irregular file: " << fullFilePath << std::endl;
        }
        
        std::string seq = readtxt_firstline(fullFilePath);
        //from the DNA sequence, calculate the energy at every j_unzipped.
        std::vector<double> seq_energy = calculate_sequence_energy(seq);
        const int seq_len = seq_energy.size();

        // convert the std::vector of double to an array of float, then copy to device
        float * seq_energy_arr = new float[seq_energy.size()];
        for (size_t i = 0; i < seq_energy.size(); ++i) {
            seq_energy_arr[i] = static_cast<float>(seq_energy[i]);
        }

        float * d_seq_energy;
        cudaMalloc((void**) &d_seq_energy, sizeof(float) * seq_energy.size());
        cudaMemcpy(d_seq_energy, seq_energy_arr, sizeof(float) * seq_energy.size(), cudaMemcpyHostToDevice);
        delete[] seq_energy_arr;

        data p = unzip(seq_len, d_seq_energy, lut); 
        cudaFree(d_seq_energy);
        
        std::string path_out = creat_path_out(fullFilePath);
        std::ofstream fout(path_out);
        if (!fout) {
            std::cerr << "Error in creating file '" << path_out << "'. Cannot save result" << std::endl;
            break;
        }

        fout << "total_extension_nm,DNA_extension_avg_nm,force_avg_pN,force_sd_pN,average_bp_unzipped_avg,bp_unzipped_sd" << std::endl;
        char delimit = ',';
        for(int i = 0; i < p.zmax; ++i) {
            fout << (p.res_arr)[i * SZ] << delimit;
            fout << (p.res_arr)[i * SZ + 1] << delimit;
            fout << (p.res_arr)[i * SZ + 2] << delimit;
            fout << (p.res_arr)[i * SZ + 3] << delimit;
            fout << (p.res_arr)[i * SZ + 4] << delimit;
            fout << (p.res_arr)[i * SZ + 5] << std::endl;
        }
        // Close the file
        fout.close();
        std::cout << "Result has been written to '" << path_out << "'." <<std::endl;
        delete[] (p.res_arr);
            
        auto end = Clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s for file '" << fullFilePath <<"'."<< std::endl;
    }
    
    FreeLibrary(hDLL);
    delete[] force_LUT;
    delete[] energy_LUT;
    cudaFree(lut.d_force_LUT);
    cudaFree(lut.d_energy_LUT);
    return 0;
}