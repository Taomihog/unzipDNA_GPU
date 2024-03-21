#include "constants.h"

// CUDA kernel
// Partition function calculation to get a single point on the unzipping curve
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
        }
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
}

// run the kernel from a wrapper
#ifdef DLL
// the output function
extern "C" __declspec(dllexport) data unzip(int seq_len, const float * d_seq_energy, LUT lut) {
#else
data unzip(int seq_len, const float * d_seq_energy, LUT lut) {
#endif
    //allocate mem for the result array
    float *d_res_arr;
    const int zmax = 1.1 * (2 * L0SS * seq_len + L0DS * ARMLENGTH);
    cudaMalloc((void**) &d_res_arr, sizeof(float) * SZ * zmax);

    // run kernel then push the result to queue
    int max_thread = NTHREAD_1D;
    int max_block = (zmax+ max_thread -1) / max_thread;
    kernel<<<max_block, max_thread>>>(seq_len - 1, d_seq_energy, lut.LUT_j_dim, lut.LUT_z_dim, lut.d_force_LUT, lut.d_energy_LUT, d_res_arr, zmax);
    
    float *res_arr = new float [SZ * zmax];
    cudaMemcpy(res_arr, d_res_arr, sizeof(float) * SZ * zmax, cudaMemcpyDeviceToHost);
    cudaFree(d_res_arr);

    return {zmax, res_arr};
}

#ifndef DLL
int main(){
    // test code 
    return 0;
}
#endif