#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <string>
#include "constants.h"
#include "utils.h"
#include "ThreadPool.h"
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

#define double float
#define MAX_PARAM_SIZE 20
#define MAX_DNA_LEN 10000
#define VALIDMAXFORCE 1000000000.0f
#define VALIDMINFORCE 0.2f
#define TOR_BINARY_SEARCH 1.0e-10
#define VERYLARGENUMBER 1.0e20

__constant__ float d_parameters[MAX_PARAM_SIZE]; // I assume the __constant__ memory is initialized to zero
__constant__ float d_sequence_energy[MAX_DNA_LEN]; // I assume the __constant__ memory is initialized to zero



__device__ double Langevin(double x) {
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

__device__ double Langevin_integ(double x) {
    // = ln(sinh(x)/x)
    double sum = 0.0;
    double factor = 1.0;
    double product = 1.0;
    for (double i = 1.0; i < 20.0; ++i) {
        sum += factor;
        factor *= (x * x /(i * 2.0) / (i * 2.0 + 1.0));
    }
    //return sum;
    return log(sum);
}

constexpr double alpha2phi_Odijk95(double alpha, double k0_eff) { 
    // if (alpha < 0.25) {
    //     return 0.0;// undefined at alpha == 0, just give it a small value
    //     //I can do this because this is force, energy must be calculated correctly!!
    // }
    return 1.0 - 0.5 / sqrt(alpha) + alpha / k0_eff; 
}
constexpr double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
constexpr double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}
// ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
//Modified version specific for ssDNA force region, and keeps accuracy
//For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
//My homemade Langevin_integ function should be accurate enough in this region.
constexpr double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
    return Langevin(2.0 * alpha) + alpha / k0_eff;
}
constexpr double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
    return 0.5 * Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
}
constexpr double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
    //integ actually starts from 1, but it's OK since it is for partition function calculation
    return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
}

constexpr double lz_ds (double force) {//dsDNA's length per base
    return ARMLENGTH * L0DS * 
            alpha2phi_Odijk95(force * LPDS / KT, KDS * LPDS / KT);
}
constexpr double lz_ss (double force, int j) {//ssDNA's length per base
    return 2.0 * j * L0SS * 
            alpha2phi_Smith95_m(force * LPSS / KT, KSS * LPSS / KT);
}
constexpr double le_ds (double force) {//function version of dsDNA's energy per bp:
    return ARMLENGTH * KT * L0DS * 
            integ_alphadphi_Odijk95(force * LPDS / KT, KDS * LPDS / KT) / LPDS;
}
constexpr double le_ss (double force, int j) {//function version of ssDNA's energy per bp:
    return 2.0 * j * KT * L0SS * 
            integ_alphadphi_Smith95_m(force * LPSS / KT, KSS * LPSS / KT) / LPSS;
}
constexpr double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
    return ext - force/PILLARSTIFFNESS - lz_ds (force) - lz_ss (force, j);//increasing function with force
}

__device__ void get_force_and_energy(int j, int z, float * force, float * energy) {
    //j == length of unzipped trunk, ext total extension
    //simple binary search to get force so calc_z_diff(force) == 0
    //force function must be monotonic
    double f1 = VALIDMINFORCE;
    double f2 = VALIDMAXFORCE;



    float ext = 1.0f;



    double y1 = delta_ext(f1, j, ext);
    double y2 = delta_ext(f2, j, ext);

    if (y1 * y2 >= 0) {
        if (y1 < 0){
            *force = VALIDMAXFORCE + 1.0;//force is too large
            return;
        } else {
            *force =  VALIDMINFORCE - 1.0;//force is too small
            return;
        }
    }

    int cnt = 0;//in case the root is not found
    while (++cnt < 10000) {
        
        double fm = (f1 + f2) * 0.5;
        double ym = delta_ext(fm, j, ext);
        
        if (abs(ym) <= TOR_BINARY_SEARCH) {
            *force =  fm;
            return;
        }

        //(&& has higher precedence)
        if (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) {
            f2 = fm;
            y2 = ym;
        } else if (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) {
            f1 = fm;
            y1 = ym;
        } else {
            *force =  VERYLARGENUMBER;//means some weird error
            return;
        }
    }
    *force =  VERYLARGENUMBER;//meaning that the root is not found
    return;
}

__global__ void kernel(float * data) {
    
    // map from threadIdx/BlockIdx to pixel position
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    int z = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = j + z * blockDim.x * gridDim.x;
    get_force_and_energy(j,z, data + offset * 2, data +offset * 2 + 1);
    __syncthreads();

    // TODO: the summation can also be on device, but I am not be able to tell if it can optimize the performance. 
}

const double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;
dp calculate_array(int extension, const std::vector<double> & seq_energy) {
        std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        temp_f.at(j) = Force(j,extension);
        temp_e.at(j) = Energy(j,extension);
        // temp_f.at(j) = lookup(j,extension,Lut_force);
        // temp_e.at(j) = lookup(j,extension,Lut_energy);
        
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
        temp1 = temp1 > energy_threshold ? 0.0 : std::exp(-temp1);
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
    point.extension_total = extension;
    point.force_average = Fprob/prob;
    point.extension_DNA = extension - point.force_average/PILLARSTIFFNESS;
    point.force_SD = std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob));
    point.junzipped_average = Jprob/prob;
    point.junzipped_SD = std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob));
    return point;
}

extern "C" __declspec(dllexport) void unzip_curve(float *sequence_energy, int * parameters, int size, float *fecurve, ) {

    cudaMemcpyToSymbol(d_parameters, parameters, sizeof(float) * PARAM_SIZE);
    cudaMemcpyToSymbol(d_sequence_energy, sequence_energy, sizeof(float) * DNA_LEN);
                                

    // parse parameters


    int DIM = size;
    dim3    blocks(DIM/16,DIM/16);
    dim3    threads(16,16);

    // should create constant memory to store the sequence_energy and DNA parameters

    cudaMalloc(&d_sequence_energy, size * sizeof(float));
    cudaMalloc(&d_parameters, parameter_size * sizeof(float));
    cudaMemcpy(d_sequence_energy, sequence_energy, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sequence_energy, sequence_energy, size * sizeof(float), cudaMemcpyHostToDevice);

    kernel<<<blocks,threads>>>( d->dev_bitmap, ticks );


    for (int j = 0; j < seq_energy.size(); ++j) {
        temp_f.at(j) = Force(j,extension);
        temp_e.at(j) = Energy(j,extension);
        // temp_f.at(j) = lookup(j,extension,Lut_force);
        // temp_e.at(j) = lookup(j,extension,Lut_energy);
        
        temp_e.at(j) += seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    float *d_input, *d_output;
    cudaMalloc(&d_input, size * sizeof(float));
    cudaMalloc(&d_output, size * sizeof(float));
    cudaMemcpy(d_input, input, size * sizeof(float), cudaMemcpyHostToDevice);
    
    int threadsPerBlock = 256;
    int numBlocks = (size + threadsPerBlock - 1) / threadsPerBlock;
    kernel<<<numBlocks, threadsPerBlock>>>(d_input, d_output, size);
    
    cudaMemcpy(output, d_output, size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(d_input);
    cudaFree(d_output);

    //calculate the partition functions:
    {
        const std::string sequence = readtxt_firstline(argv[1]);
        const std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
        //test, the "41740.760955375" is from my python code
        //std::cout << "Sequence energy differs from python program by "<< *(seq_energy.cend() - 1) - 10140.0933068047 << std::endl;// no difference, good
        
        int numThreads = std::thread::hardware_concurrency();
        std::cout << "Number of threads: " << numThreads << std::endl;
        ThreadPool pool(numThreads);

        std::vector< std::future< dp > > results;

        for(int extension = 1; extension < static_cast<int>(1.2 * seq_energy.size()); ++extension) {
            results.emplace_back(
                pool.enqueue([extension, &seq_energy]{
                    return calculate_array(extension, seq_energy);
                })
            );
        }

        std::vector<dp> result_array;
        for(auto && result: results)
            result_array.emplace_back(result.get());

        //prepare file for output
        const std::string path_out = argc >= 3 ? argv[2] : creat_path_out(argv[1]);

        std::ofstream fout(path_out);
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
        std::cout << "Result has been written to '" + path_out << "'." <<std::endl;
        
        auto end = Clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;
    }
}