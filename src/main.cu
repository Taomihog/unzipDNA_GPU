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

#ifndef NOTEST
#define JJ 100
#define ZZ 800
// f should be 16.3295 at 1600
#define EXT1 1600
// f should be 13.8639 at 4402
#define EXT2 4402
#define EXT3 2140
#endif
// the x dimension of the res_arr
#define SZ 6

using Clock = std::chrono::high_resolution_clock;

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
// Partition function calculation
__global__ void kernel(int seq_len, const float * seq_energy, int LUT_j_dim, int LUT_z_dim, const float * force_LUT, const float * energy_LUT, float * res, int zmax) {
    // z (= extension) is on x-axis only
    // In the future, we may expand along y-axis the calculation of difference sequences, like:
    // int seq_index = threadIdx.y + blockIdx.y * blockDim.y;
    // int offset = z + seq_index * blockDim.x * gridDim.x;
    int z = threadIdx.x + blockIdx.x * blockDim.x;
    if (z >= LUT_z_dim || z >= zmax) {
        return;
    }

    float  temp_e, temp_f;
    float min_e = 1.0e20;
    // remember that the LUT_j_dim is larger than j's range of this sequence, i.e., seq_len.
    for (int j = 0; j < seq_len; ++j) {
        temp_e = energy_LUT[j + z * LUT_j_dim] + seq_energy[j];
        if (min_e > temp_e) {
            min_e = temp_e;
            
            if (z == EXT3-1) { 
                printf("===z = %d, LUT_z_dim = %d, zmax = %d, j = %d, min_e = %f, energy_LUT[j + z * LUT_j_dim] = %f, seq_energy[j] = %f====\n", 
                            z, LUT_z_dim, zmax, j, min_e, energy_LUT[j + z * LUT_j_dim], seq_energy[j]);
            }
        }
    }

    if (z == EXT3) { 
        printf("\n===z = %d, LUT_z_dim = %d, zmax = %d, min_e = %f====\n\n", 
                    z, LUT_z_dim, zmax, min_e);
    }

    float prob = 0;
    float Fprob = 0;
    float FFprob = 0;
    float Jprob = 0;
    float JJprob = 0;
    float temp_p,temp;
    for (int j = 0; j < seq_len; ++j) {
        temp_e = energy_LUT[j + z * LUT_j_dim] + seq_energy[j] - min_e;
        temp_f = force_LUT[j + z * LUT_j_dim];

        temp_p = temp_e > ENERGY_THRESHOLD ? 0.0f : exp(-temp_e);

        prob += temp_p;

        temp = temp_f * temp_p;
        Fprob += temp;
        FFprob += temp * temp_f;

        temp = j * temp_p;
        Jprob += temp;
        JJprob += j * temp;
    }

    res[z * SZ    ] = z; // extension_total
    float f_avg = Fprob/prob; // force_average
    res[z * SZ + 1] = z - f_avg/PILLARSTIFFNESS; // extension_DNA
    res[z * SZ + 2] = f_avg; // force_average
    res[z * SZ + 3] = sqrt(FFprob/prob - f_avg * f_avg); // force_SD
    float j_avg = Jprob/prob; // j_unzipped_average
    res[z * SZ + 4] = j_avg; // j_unzipped_average
    res[z * SZ + 5] = sqrt(JJprob/prob - j_avg * j_avg); // junzipped_SD
    // printf("z = %d, min_e = %f, prob = %f, Fprob = % f, F_avg = %f.\n", z, min_e, prob, Fprob, Fprob/prob);
#ifndef NOTEST
    if (z == EXT1) { printf("At z = %d, f_avg = %f.\n", z, f_avg); }
    if (z == EXT2) { printf("At z = %d, f_avg = %f.\n", z, f_avg); }
    if (z == EXT3) { printf("At z = %d, f_avg = %f.\n", z, f_avg); }
#endif
}
struct data {
    int zmax;
    float * res_arr;
};
struct LUT {
    int LUT_j_dim;
    int LUT_z_dim;
    float * d_force_LUT; 
    float * d_energy_LUT;
};
std::mutex mtx;
std::condition_variable cv;
std::queue<data> dataQueue;
bool finished = false;
// producer 
void calculate(LUT lut) {

    std::string path {"../test_data/NEB_H5alpha_Accessory_colonization_factor_AcfD.txt"};

    for (int data_idx = 0; data_idx < 1; ++data_idx) {
        std::lock_guard<std::mutex> lock(mtx);
        auto start = Clock::now();

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
        printf("zmax = %d.\n", zmax);
        cudaMalloc((void**) &d_res_arr, sizeof(float) * SZ * zmax);

        // run kernel then push the result to queue
        int max_thread = 128;
        int max_block = (zmax+ max_thread -1) / max_thread;
        kernel<<<max_block, max_thread>>>(seq_energy.size() - 1, d_seq_energy_arr, lut.LUT_j_dim, lut.LUT_z_dim, lut.d_force_LUT, lut.d_energy_LUT, d_res_arr, zmax);
        cudaFree(d_seq_energy_arr);
        float *res_arr = new float [SZ * zmax];
        cudaMemcpy(res_arr, d_res_arr, sizeof(float) * SZ * zmax, cudaMemcpyDeviceToHost);
        cudaFree(d_res_arr);
        
        dataQueue.push({zmax, res_arr}); // Push data to the queue
#ifndef NOTEST
            // make some test
            printf("producer: res_arr[%d] = %f\n",EXT1, res_arr[SZ * EXT1 + 2]);
            printf("producer: res_arr[%d] = %f\n",EXT2, res_arr[SZ * EXT2 + 2]);
            printf("producer: res_arr[%d] = %f\n",EXT3, res_arr[SZ * EXT3 + 2]);
#endif

        auto end = Clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;
        
        cv.notify_one(); // Notify the consumer
        //std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    finished = true;
    cv.notify_one(); // Notify the consumer that production is finished
}
// consumer 
void save(int LUT_j_dim) {
    while (!finished) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, []{ return !dataQueue.empty() || finished; }); // Wait until there is data or production is finished
        // cv.wait(lock, []{ return finished; }); // Wait until there is data or production is finished
        if (!dataQueue.empty()) {
            data p = dataQueue.front(); // Get data from the queue
            dataQueue.pop(); // Remove data from the queue
            
#ifndef NOTEST
            // make some test
            printf("consumer: p.res_arr[%d] = %f\n",EXT1, p.res_arr[SZ * EXT1 + 2]);
            printf("consumer: p.res_arr[%d] = %f\n",EXT2, p.res_arr[SZ * EXT2 + 2]);
            printf("consumer: p.res_arr[%d] = %f\n",EXT3, p.res_arr[SZ * EXT3 + 2]);
#endif

            std::ofstream fout("res.csv");
            if (!fout) {
                std::cerr << "Error opening file. Cannot save result" << std::endl;
                break;
            }

            fout << "total_extension_nm, DNA_extension_avg_nm, force_avg_pN, force_sd_pN, average_bp_unzipped_avg, bp_unzipped_sd" << std::endl;

            char delimit = ',';
            for(int i = 0; i < p.zmax; ++i) {
                // if ((p.res_arr)[i * SZ] == 0.0f) { 
                //     break;
                // }
                fout << (p.res_arr)[i * SZ] << delimit;
                fout << (p.res_arr)[i * SZ + 1] << delimit;
                fout << (p.res_arr)[i * SZ + 2] << delimit;
                fout << (p.res_arr)[i * SZ + 3] << delimit;
                fout << (p.res_arr)[i * SZ + 4] << delimit;
                fout << (p.res_arr)[i * SZ + 5] << std::endl;
            }

            // Close the file
            fout.close();
            std::cout << "Result has been written to res.csv" << "'." <<std::endl;

            delete[] (p.res_arr);
        }
    }
}
// main
int main() {
    // LUT's size along j and z directions, and:
    //      j = threadIdx.x + blockIdx.x * blockDim.x
    //      z = threadIdx.y + blockIdx.y * blockDim.y
    constexpr int LUT_j_dim = 4096;
    constexpr int LUT_z_dim = 4096;
    
    float *force_LUT;
    float *energy_LUT;
    force_LUT = new float [LUT_j_dim * LUT_z_dim];
    energy_LUT = new float [LUT_j_dim * LUT_z_dim];

    // this scope is common, will not be a thread
    {
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

        get_LUT(LUT_j_dim, LUT_z_dim, force_LUT, energy_LUT);
        FreeLibrary(hDLL);

#ifndef NOTEST
        // make some test
        printf("force[%d][%d] = %f\n", JJ,ZZ, force_LUT[JJ + ZZ * LUT_j_dim]);
        printf("energy[%d][%d] = %f\n", JJ,ZZ, energy_LUT[JJ + ZZ * LUT_j_dim]);
#endif
    }

    float *d_force_LUT;
    float *d_energy_LUT;
    cudaMalloc((void**) &d_force_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim);
    cudaMalloc((void**) &d_energy_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim);
    cudaMemcpy(d_force_LUT, force_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim, cudaMemcpyHostToDevice);
    cudaMemcpy(d_energy_LUT, energy_LUT, sizeof(float) * LUT_j_dim * LUT_z_dim, cudaMemcpyHostToDevice);
    LUT lut {LUT_j_dim, LUT_z_dim, d_force_LUT, d_energy_LUT};
    std::thread producerThread(calculate, lut);
    std::thread consumerThread(save, LUT_j_dim);

    producerThread.join();
    consumerThread.join();

    delete[] force_LUT;
    delete[] energy_LUT;
    cudaFree(d_force_LUT);
    cudaFree(d_energy_LUT);
    return 0;
}