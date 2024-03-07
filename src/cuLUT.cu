/*
cmd to compile DLL file:
nvcc -c -o cuLUT.o cuLUT.cu -DDLL
nvcc -shared -o cuLUT.dll cuLUT.o -lcudart

objdump: A good tool of gcc users to search for the addresses of functions, command:
objdump -x cuLUT.dll
*/
#include <cuda_runtime.h>
#include "../include/constants.h"

/*
   --- General Information for device 0 ---
Name:  NVIDIA GeForce RTX 3070 Laptop GPU
Compute capability:  8.6
Clock rate:  1290000
Device copy overlap:  Enabled
Kernel execution timeout :  Enabled
   --- Memory Information for device 0 ---
Total global mem:  8589410304
Total constant Mem:  65536
Max mem pitch:  2147483647
Texture Alignment:  512
   --- MP Information for device 0 ---
Multiprocessor count:  40
Shared mem per mp:  49152
Registers per mp:  65536
Threads in warp:  32
Max threads per block:  1024
Max thread dimensions:  (1024, 1024, 64)
Max grid dimensions:  (2147483647, 65535, 65535)
*/

#ifndef DLL
#include <stdio.h>
#define JJ 100
#define ZZ 800
#endif

#define double float
// Calculate Langevin function and its integration from the Taylor expansions
__device__ inline double Langevin(double x) {
    if ( x > 7.9608220) {
        // I compared Coth(x) and the np.cosh(x)/np.sinh(x), and 1, and figured out this value.
        //above this value, 1.0 is more close to Coth(x).
        //this value depends on how many terms are used of course.
        return 1.0 - 1.0 / x;
    }
    double product = x*x;//x**0/0!
    double factorial = 6;//2!
    double sum1 = 1.0/2.0 - 1.0 / 6.0;
    double sum2 = 1.0;
    for (int i = 2; i <= 10; ++i){
        sum2 += product /factorial;
        factorial *= 2.0 * i;
        sum1 += product * (1.0/factorial - 1.0 /factorial/(2.0 * i + 1.0));
        factorial *= (2.0 * i + 1.0);
        product *= x*x;
    }
    return x*sum1/sum2;
}
__device__ inline double Langevin_integ(double x) {
    // = ln(sinh(x)/x)
    double sum = 0.0;
    double factor = 1.0;
    for (double i = 1.0; i < 20.0; ++i) {
        sum += factor;
        factor *= (x * x /(i * 2.0) / (i * 2.0 + 1.0));
    }
    //return sum;
    return log(sum);
}
//================================WLC high force, Odijk 1995 macromolecules============================
__device__ inline double alpha2phi_Odijk95(double alpha, double k0_eff) { 
    // if (alpha < 0.25) {
    //     return 0.0;// undefined at alpha == 0, just give it a small value
    //     //I can do this because this is force, energy must be calculated correctly!!
    // }
    return 1.0 - 0.5 / sqrt(alpha) + alpha / k0_eff; 
}
__device__ inline double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
__device__ inline double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}
// ========================MODIFIED VERSION OF FJC, Smith 1995 macromolecules==========================
//Modified version specific for ssDNA force region, and keeps accuracy
//For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
//My homemade Langevin_integ function should be accurate enough in this region.
__device__ inline double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
    return Langevin(2.0 * alpha) + alpha / k0_eff;
}
__device__ inline double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
    return 0.5 * Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
}
__device__ inline double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
    //integ actually starts from 1, but it's OK since it is for partition function calculation
    return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
}
// ======================helper for root finding, force that total extension = ext=====================
__device__ inline double ext_at_this_force(double force, double j) {
    return  force/PILLARSTIFFNESS + 
            ARMLENGTH * L0DS * alpha2phi_Odijk95(force * LPDS / KT, KDS * LPDS / KT) + 
            2.0 * j * L0SS * alpha2phi_Smith95_m(force * LPSS / KT, KSS * LPSS / KT);//increasing function with force
}
__global__ void kernel(float * force, float * energy) {
    // map from threadIdx/BlockIdx to pixel position
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    int z = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = j + z * blockDim.x * gridDim.x;

    float f = 0.0f;//meaning that the root is not found
    float e = 0.0f;

    //j == length of unzipped trunk, ext total extension
    //simple binary search to get force so calc_z_diff(force) == 0
    //force function must be monotonic
    double f1 = VALIDMINFORCE;
    double f2 = VALIDMAXFORCE;

    double y1 = z - ext_at_this_force(f1, j);
    double y2 = z - ext_at_this_force(f2, j);

#ifdef JJ
#ifdef ZZ
    if (j == JJ && z == ZZ) {
        printf("CUDA print j = %d, z = %d\n", j, z);
        printf("CUDA print z1 = %f, z2 = %f\n",ext_at_this_force(f1, j) , ext_at_this_force(f2, j));
        printf("CUDA print z-z1 = %f, z-z2 = %f\n", y1, y2);
    }
#endif
#endif

    if (y1 * y2 >= 0) {
        if (y1 < 0){
            f = VALIDMAXFORCE * 2.0;//force is too large
        } else {
            f =  VALIDMINFORCE / 2.0;//force is too small
        }
    } else {
        int cnt = 0;//in case the root is not found
        double fm, ym;
        while (++cnt < 10000) {
            
            fm = (f1 + f2) * 0.5;
            ym = z - ext_at_this_force(fm, j);
            
            if (abs(ym) <= TOR_BINARY_SEARCH) {
                f =  fm;
                break;
            } else if (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) {
                f2 = fm;
                y2 = ym;
            } else if (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) {
                f1 = fm;
                y1 = ym;
            } else {
                f =  -1.0f;//means some weird error
                printf("error");
                break;
            }
        }
    }

    // energy normalized by kT
    e =   0.5 * f * f / PILLARSTIFFNESS / KT + 
                ARMLENGTH  * L0DS * integ_alphadphi_Odijk95(f * LPDS / KT, KDS * LPDS / KT) / LPDS + 
                2.0 * j * L0SS * integ_alphadphi_Smith95_m(f * LPSS / KT, KSS * LPSS / KT) / LPSS;

#ifdef JJ
#ifdef ZZ
    if (j == JJ && z == ZZ) { printf("CUDA print f = %f, e = %f\n", f, e); }
#endif
#endif

    force[offset] = f;
    energy[offset] = e;

    // __syncthreads();

    // TODO: the summation can also be on device, but I am not be able to tell if it can optimize the performance. 
    // It probably cannot, and generate a lot of overhead on GPU (because the calculation will be sequence-dependent then).
}
#undef double

#ifdef DLL
// the output function
extern "C" __declspec(dllexport) void get_LUT(int X_DIM, int Y_DIM, float * force_out, float * energy_out) {
#else
void get_LUT(int X_DIM, int Y_DIM, float * force_out, float * energy_out) {
#endif
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    // x_index is j, y_index is z, can
    dim3    threads(NTHREAD, NTHREAD);
    dim3    blocks(X_DIM/NTHREAD, Y_DIM/NTHREAD);

    printf("Total size (j x z):     %d x %d.\n", X_DIM, Y_DIM);
    printf("Threads dim:            %d x %d.\n", NTHREAD, NTHREAD);
    printf("Blocks dim:             %d x %d.\n", X_DIM/NTHREAD, Y_DIM/NTHREAD);

    float * d_force, * d_energy;
    cudaMalloc((void**)&d_force, X_DIM * Y_DIM * sizeof(float));
    cudaMalloc((void**)&d_energy, X_DIM * Y_DIM * sizeof(float));

    kernel<<<blocks,threads>>>(d_force, d_energy);

    cudaMemcpy(force_out, d_force, X_DIM * Y_DIM * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(energy_out, d_energy, X_DIM * Y_DIM * sizeof(float), cudaMemcpyDeviceToHost);
    cudaEventRecord(stop, 0 );
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop );
    printf("Time to generate: %3.1f s\n", elapsedTime/1000);
}

#ifndef DLL
int main() {
    int max_j_unzipped = 1000;
    int X_DIM = (max_j_unzipped + NTHREAD - 1);
    X_DIM -= X_DIM % NTHREAD;
    int Y_DIM = (int)((max_j_unzipped + ARMLENGTH) * L0SS * 2) + NTHREAD - 1;
    Y_DIM -= Y_DIM % NTHREAD;
    
    float *force = new float [X_DIM * Y_DIM];
    float *energy = new float [X_DIM * Y_DIM];

    get_LUT(X_DIM, Y_DIM, force, energy);

    printf("force[%d][%d] = %f\n", JJ,ZZ, force[JJ + ZZ * X_DIM]);
    printf("energy[%d][%d] = %f\n", JJ,ZZ, energy[JJ + ZZ * X_DIM]);

    // float temp;
    // get_force(100, 1000, &temp);
    // printf("force %f\n", temp);

    delete[] force, energy;
    return 0;
}
#endif