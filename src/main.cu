// compile cmd:

// make DLL(s):
// nvcc -c -o cuLUT.o cuLUT.cu -DDLL
// nvcc -shared -o cuLUT.dll cuLUT.o -lcudart

// (objdump: A good tool of gcc users to search for the addresses of functions)
// objdump -x cuLUT.dll

// make static libraries:
// g++ -c  parser.cpp -o parser.o
// ar rcs lib_parser.a parser.o

// make main.exe
// nvcc main.cu -o main.exe -L. -l_parser -DTEST

// Do I need this?
//#include <cuda_runtime.h>
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
#include "../include/ThreadPool.h"//https://github.com/progschj/ThreadPool
*/

#ifndef NOTEST
#define JJ 100
#define ZZ 800
#endif

using Clock = std::chrono::high_resolution_clock;

__global__ void kernel(int seq_len, const float * seq_energy, int LUT_j_DIM, int LUT_z_DIM, const float * force_LUT, const float * energy_LUT, float * res, int zmax) {
    // z (= extension) is on x-axis only
    int z = threadIdx.x + blockIdx.x * blockDim.x;
    if (z >= LUT_z_DIM || z >= zmax) {
        return;
    }
    // we may expand the calculation of difference sequences on y-axis in the future
    // int seq_index = threadIdx.y + blockIdx.y * blockDim.y;
    // int offset = z + seq_index * blockDim.x * gridDim.x;

    float * temp_e, *temp_f;
    cudaMalloc((void**)&temp_e, seq_len);
    cudaMalloc((void**)&temp_f, seq_len);

    // remember that the LUT_j_DIM may be different from j's range of this sequence, i.e., seq_len.
    // When we calculate the LUTs, we use larger j and z ranges, but here the j's range depends on the trunk sequence length

    float min_e = 1.0e20;
    for (int j = 0; j < seq_len; ++j) {
        temp_f[j] = force_LUT[j + z * LUT_j_DIM];
        temp_e[j] = energy_LUT[j + z * LUT_j_DIM] + seq_energy[j];
        
        if (min_e > temp_e[j]) {
            min_e = temp_e[j];
        }
    }
    printf("z = %d\n", z);

    float prob = 0;
    float Fprob = 0;
    float FFprob = 0;
    float Jprob = 0;
    float JJprob = 0;
    float temp1,temp2,temp3;
    for (int j = 0; j < seq_len; ++j) {

        temp1 = temp_e[j] - min_e;
        temp1 = temp1 > ENERGY_THRESHOLD ? 0.0 : exp(-temp1);
        temp2 = temp_f[j];

        prob += temp1;

        temp3 = temp2 * temp1;
        Fprob += temp3;
        FFprob += temp3 * temp2;

        temp3 = j * temp1;
        Jprob += temp3;
        JJprob += j * temp3;

    }
    // extension_total
    res[z * 6    ] = z;

    // force_average
    float f_avg = Fprob/prob;
    res[z * 6 + 1] = f_avg;
    // extension_DNA
    res[z * 6 + 2] = z - f_avg/PILLARSTIFFNESS;
    // force_SD
    res[z * 6 + 3] = sqrt(FFprob/prob - f_avg * f_avg);

    // junzipped_average
    float j_avg = Jprob/prob;
    res[z * 6 + 4] = j_avg;
    // junzipped_SD
    res[z * 6 + 5] = sqrt(JJprob/prob - j_avg * j_avg);

    if (abs(z - 100) == 0) {
        printf("At z = %d, f_avg = %f.\n", z, f_avg);
    }
    printf("At z = %d, f_avg = %f.\n", z, f_avg);

    cudaFree(temp_e);
    cudaFree(temp_f);
}

std::mutex mtx;
std::condition_variable cv;
std::queue<float*> dataQueue;
bool finished = false;

// producer 
void calculate() {
    for (int i = 0; i < 1; ++i) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            std::cout << "Producing: " << i << std::endl;
            float * p = new float[2];// put the function here
            p[0] = i;
            p[1] = i * i;
            dataQueue.push(p); // Push data to the queue
        }
        cv.notify_one(); // Notify the consumer
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simulate work
    }
    finished = true;
    cv.notify_one(); // Notify the consumer that production is finished
}

// consumer 
void save() {
    while (!finished) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, []{ return !dataQueue.empty() || finished; }); // Wait until there is data or production is finished
        // cv.wait(lock, []{ return finished; }); // Wait until there is data or production is finished
        if (!dataQueue.empty()) {
            float* p = dataQueue.front(); // Get data from the queue
            dataQueue.pop(); // Remove data from the queue
            std::cout << "Consuming: " << p[0] << ", " << p[1] << std::endl;
            delete[] p; // I should never do someting like this "++p", then delete[] p; otherwise the compiler can lose track of the array size!
        }
    }
}


int main() {
    // LUT's size along j and z directions, def:
    //      j: LUT_J_DIM, threadIdx.x + blockIdx.x * blockDim.x
    //      z: LUT_Z_DIM, threadIdx.y + blockIdx.y * blockDim.y
    constexpr int max_j_unzipped = 8000;
    // constexpr int LUT_J_DIM = ((max_j_unzipped + NTHREAD - 1) / NTHREAD) * NTHREAD;
    // constexpr int LUT_Z_DIM = (((int)((max_j_unzipped + ARMLENGTH) * L0SS * 2) + NTHREAD - 1) / NTHREAD) * NTHREAD;
    constexpr int LUT_J_DIM = 8000;
    constexpr int LUT_Z_DIM = (((int)((max_j_unzipped + ARMLENGTH) * L0SS * 2) + NTHREAD - 1) / NTHREAD) * NTHREAD;
    
    float *force_LUT;
    float *energy_LUT;
    force_LUT = new float [LUT_J_DIM * LUT_Z_DIM];
    energy_LUT = new float [LUT_J_DIM * LUT_Z_DIM];

    // this scope is common, will not be a thread
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

        get_LUT(LUT_J_DIM, LUT_Z_DIM, force_LUT, energy_LUT);
        FreeLibrary(hDLL);

#ifndef NOTEST
        // make some test
        // float *testf = new float;
        // float *teste = new float;
        // cudaMemcpy(&testf, d_force_LUT+ JJ + ZZ * LUT_J_DIM, sizeof(float), cudaMemcpyDeviceToHost);
        // cudaMemcpy(&testf, d_force_LUT+ JJ + ZZ * LUT_J_DIM, sizeof(float), cudaMemcpyDeviceToHost);
        printf("force[%d][%d] = %f\n", JJ,ZZ, force_LUT[JJ + ZZ * LUT_J_DIM]);
        printf("energy[%d][%d] = %f\n", JJ,ZZ, energy_LUT[JJ + ZZ * LUT_J_DIM]);
#endif
    }

    float *d_force_LUT;
    float *d_energy_LUT;
    cudaMalloc((void**) &d_force_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM);
    cudaMalloc((void**) &d_energy_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM);
    cudaMemcpy(d_force_LUT, force_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(d_energy_LUT, energy_LUT, sizeof(float) * LUT_J_DIM * LUT_Z_DIM, cudaMemcpyHostToDevice);

    std::thread producerThread(calculate);
    std::thread consumerThread(save);

    producerThread.join();
    consumerThread.join();

    // I should use queue;
    float ** data_ptr = new float * [10000];
    int * zmax_ptr = new int[10000];
    int data_idx = 0;
    // suppose to start a loop here
    std::string path {"../test_data/NEB_H5alpha_Accessory_colonization_factor_AcfD.txt"};
    {
        auto start = Clock::now();

        //IO and parser
        // std::string seq = parser::readtxt_firstline(path);

        std::ifstream file(path);

        if (!file.is_open()) {
            std::cerr << "Error opening the file." << std::endl;
            return;
        }
        std::string seq;
        if (getline(file, seq)) {
            std::cout << "Trunk sequence '" + path + "' read successfully." << std::endl;
        } else {
            std::cerr << "Failed to read a line from the file." << std::endl;
        }

        file.close();

        //from the DNA sequence, calculate the energy at every j_unzipped.

        std::vector<double> seq_energy;
        if (seq.size() < 2) {
            std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
            return;
        }

        double accum = 0.0;
        *std::back_inserter(seq_energy) = 0.0;
        std::transform(seq.cbegin() + 1, seq.cend(), std::back_inserter(seq_energy), 
        [&accum](const char& bp2) {
            accum += BPEnergy::lookup_bp_energy(*(&bp2 - 1), bp2);
            return accum;
            }
        );//equivalent to std::partial_sum()

        //test, the "43048.29409015185" is from my python code
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

        float *d_res_arr;
        const int zmax = 1.2 * (seq_energy.size() + ARMLENGTH);
        cudaMalloc((void**) &d_res_arr, sizeof(float) * 6 * zmax);

        int max_thread = 128;
        int max_block = (zmax+ max_thread -1) / max_thread;
        kernel<<<max_block, max_thread>>>(seq_energy.size() - 1, d_seq_energy_arr, LUT_J_DIM, LUT_Z_DIM, d_force_LUT, d_energy_LUT, d_res_arr, zmax);
        cudaFree(d_seq_energy_arr);

        float *res_arr = new float [6 * zmax];
        cudaMemcpy(res_arr, d_res_arr, sizeof(float) * 6 * zmax, cudaMemcpyDeviceToHost);
        cudaFree(d_res_arr);

        data_ptr[data_idx] = res_arr;
        zmax_ptr[data_idx] = zmax;
        auto end = Clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;
    }



    // data_idx = 0;
    // float * res_arr = data_ptr[data_idx];
    // int zmax = zmax_ptr[data_idx];
    // // this is another thread
    // {

    //     std::ofstream fout("res.csv");
    //     if (!fout) {
    //         std::cerr << "Error opening file. Cannot save result" << std::endl;
    //         return 1;
    //     }

    //     fout << "total_extension_nm, DNA_extension_avg_nm, force_avg_pN, force_sd_pN, average_bp_unzipped_avg, bp_unzipped_sd" << std::endl;

    //     char delimit = ',';
    //     for(int i = 0; i < zmax; ++i) {
    //         // if (res_arr[i * 6] == 0.0f) { 
    //         //     break;
    //         // }
    //         fout << res_arr[i * 6] << delimit;
    //         fout << res_arr[i * 6 + 1] << delimit;
    //         fout << res_arr[i * 6 + 2] << delimit;
    //         fout << res_arr[i * 6 + 3] << delimit;
    //         fout << res_arr[i * 6 + 4] << delimit;
    //         fout << res_arr[i * 6 + 5] << std::endl;
    //     }

    //     // Close the file
    //     fout.close();
    //     std::cout << "Result has been written to res.csv" << "'." <<std::endl;
    // }





    //free the queue of res_arr (each data_ptr[i])
    //then free the queue itself if needed data_ptr
    delete[] force_LUT;
    delete[] energy_LUT;
    cudaFree(d_force_LUT);
    cudaFree(d_energy_LUT);
    return 0;
}