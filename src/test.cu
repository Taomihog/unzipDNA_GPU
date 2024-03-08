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
// constants.h define "NTHREAD", this value is used here and the .dll source code
#include "../include/constants.h"
/*
// There's no need to polish further, I am brain bleeding now.
#include "../include/ThreadPool.h"//https://github.com/progschj/ThreadPool
*/

#define JJ 100
#define ZZ 800


//=============================================util functions==========================================
// single read. todo: read the fasta directly, multiple line read.
std::string readtxt_firstline(const std::string & path) {
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return "0";
    }

    std::string line;
    if (getline(file, line)) {
        std::cout << "Trunk sequence '" + path + "' read successfully." << std::endl;
    } else {
        std::cerr << "Failed to read a line from the file." << std::endl;
    }

    file.close();

    return line;
}
// generate a file path for output result from the input file path.
std::string creat_path_out(const std::string & path_in)
{
    size_t lastSlashIndex = path_in.rfind('/');

    std::string parentPath;
    if (lastSlashIndex != std::string::npos) {
        parentPath = path_in.substr(0, lastSlashIndex) + "/";
    }

    std::string fullname = path_in.substr(lastSlashIndex + 1, std::string::npos);

    return parentPath + "out_" + fullname.substr(0, fullname.rfind('.')) + ".csv";;
}
// energy
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
// CUDA kernel
__global__ void kernel(int seq_len, const float * seq_energy, int LUT_j_DIM, int LUT_z_DIM, const float * force_LUT, const float * energy_LUT, float * res, int zmax) {
    // z (= extension) is on x-axis only
    // In the future, we may expand along y-axis the calculation of difference sequences, like:
    // int seq_index = threadIdx.y + blockIdx.y * blockDim.y;
    // int offset = z + seq_index * blockDim.x * gridDim.x;
    int z = threadIdx.x + blockIdx.x * blockDim.x;
    if (z >= LUT_z_DIM || z >= zmax) {
        return;
    }

    float  temp_e, temp_f;
    float min_e = 1.0e20;
    // remember that the LUT_j_DIM is larger than j's range of this sequence, i.e., seq_len.
    for (int j = 0; j < seq_len; ++j) {
        temp_e = energy_LUT[j + z * LUT_j_DIM] + seq_energy[j];
        if (min_e > temp_e) {
            min_e = temp_e;
        }
    }

    float prob = 0;
    float Fprob = 0;
    float FFprob = 0;
    float Jprob = 0;
    float JJprob = 0;
    float temp_p,temp;
    for (int j = 0; j < seq_len; ++j) {
        temp_e = energy_LUT[j + z * LUT_j_DIM] + seq_energy[j] - min_e;
        temp_f = force_LUT[j + z * LUT_j_DIM];

        temp_p = temp_e > ENERGY_THRESHOLD ? 0.0f : exp(-temp_e);

        prob += temp_p;

        temp = temp_f * temp_p;
        Fprob += temp;
        FFprob += temp * temp_f;

        temp = j * temp_p;
        Jprob += temp;
        JJprob += j * temp;
    }

    res[z * 6    ] = z; // extension_total
    float f_avg = Fprob/prob; // force_average
    res[z * 6 + 1] = f_avg;
    res[z * 6 + 2] = z - f_avg/PILLARSTIFFNESS; // extension_DNA
    res[z * 6 + 3] = sqrt(FFprob/prob - f_avg * f_avg); // force_SD
    float j_avg = Jprob/prob;
    res[z * 6 + 4] = j_avg; // junzipped_average
    res[z * 6 + 5] = sqrt(JJprob/prob - j_avg * j_avg); // junzipped_SD
    // printf("z = %d, min_e = %f, prob = %f, Fprob = % f, F_avg = %f.\n", z, min_e, prob, Fprob, Fprob/prob);
    if (z == 1600 ) {
        printf("At z = %d, f_avg = %f.\n", z, f_avg);
    }
}
// main
int main() {
    // LUT's size along j and z directions, def:
    //      j: threadIdx.x + blockIdx.x * blockDim.x
    //      z: threadIdx.y + blockIdx.y * blockDim.y
    constexpr int LUT_J_DIM = 8000;
    constexpr int LUT_Z_DIM = (LUT_J_DIM + ARMLENGTH) * L0SS * 2;
    
    float *force_LUT;
    float *energy_LUT;
    force_LUT = new float [LUT_J_DIM * LUT_Z_DIM];
    energy_LUT = new float [LUT_J_DIM * LUT_Z_DIM];

    HINSTANCE hDLL = LoadLibrary("cuLUT.dll");
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

    get_LUT(LUT_J_DIM, LUT_Z_DIM, force_LUT, energy_LUT);
    FreeLibrary(hDLL);
    // make some test
    printf("force[%d][%d] = %f\n", JJ,ZZ, force_LUT[JJ + ZZ * LUT_J_DIM]);
    printf("energy[%d][%d] = %f\n", JJ,ZZ, energy_LUT[JJ + ZZ * LUT_J_DIM]);

    float *d_force_LUT;
    float *d_energy_LUT;
    cudaMalloc((void**) &d_force_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM);
    cudaMalloc((void**) &d_energy_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM);
    cudaMemcpy(d_force_LUT, force_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(d_energy_LUT, energy_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM, cudaMemcpyHostToDevice);

    std::string path {"../test_data/NEB_H5alpha_Accessory_colonization_factor_AcfD.txt"};
    std::string seq = readtxt_firstline(path);
    //from the DNA sequence, calculate the energy at every j_unzipped.
    std::vector<double> seq_energy = calculate_sequence_energy(seq);
    //test, the "43048.29409015185" is from my old python code
    std::cout << "Sequence energy differs from python program by "<< KT * (*(seq_energy.cend() - 1)) - 43048.29409015185 << std::endl;// no difference, good
    
    // convert the std::vector of double to an array of float, then copy to device
    float * seq_energy_arr = new float[seq_energy.size()];
    for (size_t i = 0; i < seq_energy.size(); ++i) {
        seq_energy_arr[i] = static_cast<float>(seq_energy[i]);
    }
    float * d_seq_energy_arr;
    cudaMalloc((void**) &d_seq_energy_arr, sizeof(float) * seq_energy.size());
    cudaMemcpy(d_seq_energy_arr, seq_energy_arr, sizeof(float) * seq_energy.size(), cudaMemcpyHostToDevice);
    delete[] seq_energy_arr;

    //allocate mem for the result array
    float *d_res_arr;
    const int zmax = 1.2 * (seq_energy.size() + ARMLENGTH);
    cudaMalloc((void**) &d_res_arr, sizeof(float) * 6 * zmax);

    // run kernel then push the result to queue
    int max_thread = 1024;
    int max_block = 10;//(zmax+ max_thread -1) / max_thread;
    kernel<<<max_block, max_thread>>>(seq_energy.size() - 1, d_seq_energy_arr, LUT_J_DIM, LUT_Z_DIM, d_force_LUT, d_energy_LUT, d_res_arr, zmax);
    cudaFree(d_seq_energy_arr);
    float *res_arr = new float [6 * zmax];
    cudaMemcpy(res_arr, d_res_arr, sizeof(float) * 6 * zmax, cudaMemcpyDeviceToHost);
    cudaFree(d_res_arr);

    delete[] force_LUT;
    delete[] energy_LUT;
    cudaFree(d_force_LUT);
    cudaFree(d_energy_LUT);
    return 0;
}