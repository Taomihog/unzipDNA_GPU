// My GPU:
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

#pragma once
#include <vector>
#include <string>
#include <cmath>

//Constants
#define JOUL 4184.0f //Joul-Calorie conversion
#define AVOGADRO 6.022E+23f 
#define BOLTZMANN 0.0138065f 
#define PNNM 1.0e21  // piconewton * nanometer

//Experiment conditions
#define TEMPERATURE 298.15 
#define KT (TEMPERATURE * BOLTZMANN) 
#define ARMLENGTH 2200 //total length of the 2 dsDNA arms, unit is base-pair.
#define PILLARSTIFFNESS 0.07406f //spring constant of the pillar/optical trap/micro-needle/etc that is used for stretching.
#define SALTCONC 100.0f  //salt concentration in mM, this the the salt concentraion or the unzip experiment 
#define EFFSALTCONC (log(SALTCONC * 0.001f) / 298.0) // 298.0 K is where the energy was measured in Huguet paper, it is not equal to Condition::temperature
#define FACTOR_PNNM (PNNM * JOUL / AVOGADRO / KT) // convert to 'pN-nm' unit;

//====================================DNA mechanical parameters======================================
//These are pseudo constants..they are salt dependent and temperature dependent
#define LPDS 51.97f //dsDNA persistence length
#define KDS  1318.0f //dsDNA elastic modulus
#define L0DS 0.338f //dsDNA contour length per bp
#define LPSS 0.765f //ssDNA persistence length
#define KSS  470.0f //ssDNA elastic modulus  
#define L0SS 0.554f //ssDNA contour length per nt

//======================================numerical calculation control========================================
#define NTHREAD 32
#define NTHREAD_1D (NTHREAD * NTHREAD)
// simulation precisions, etc.
#define VALIDMAXFORCE       1000.0f
#define VALIDMINFORCE       0.001f
#define TOR_BINARY_SEARCH   1.0e-2
#define VERYLARGENUMBER     1.0e20
//When dE > this threshold, don't calculate exp(-e/kT) and set probability to 0;
#define ENERGY_THRESHOLD    50.0f 

// ==============================================data struct================================================
// the x dimension of the res_arr
#define SZ 6
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

//====================================basepair energy measured by Huguet et al===============================
//ref: Huguet, Josep M., et al. (2010) PNAS
namespace BPEnergy {

    constexpr int LUTsize = 4;

    constexpr double LUTdH[LUTsize][LUTsize] = {
        {-7.28, -4.63, -5.21, -5.80},//aa, at, ag, ac
        {-8.31, -7.28, -8.96, -8.16},//ta, tt, tg, tc
        {-8.16, -5.80, -8.57, -10.1},//ga, gt, gg, gc
        {-8.96, -5.21, -9.66, -8.57} //ca, ct, cg, cc
    };

    constexpr double LUTdS[LUTsize][LUTsize] = {
        {-20.28, -11.62, -12.89, -14.46},//aa, at, ag, ac
        {-25.06, -20.28, -24.48, -22.46},//ta, tt, tg, tc
        {-22.46, -14.46, -22.30, -25.96},//ga, gt, gg, gc
        {-24.48, -12.89, -24.43, -22.30} //ca, ct, cg, cc
    };

    constexpr double LUTm[LUTsize][LUTsize] = {
        {0.145, 0.117, 0.070, 0.099},//aa, at, ag, ac
        {0.091, 0.145, 0.091, 0.155},//ta, tt, tg, tc
        {0.155, 0.099, 0.063, 0.079},//ga, gt, gg, gc
        {0.091, 0.070, 0.132, 0.063} //ca, ct, cg, cc
    };

    // use constexpr class to simplify indexing
    class constexpr_map_class {

        public:
        constexpr int operator[] (char key) const {
                return bp2idx (key);
        }

        private:
        constexpr int bp2idx(char base) const {
            switch (base) {
            case 'a':
            case 'A':
                return 0;

            case 't':
            case 'T':
                return 1;

            case 'g':
            case 'G':
                return 2;

            case 'c':
            case 'C':
                return 3;

            default:
                return -1;
            }
        };
    };

    constexpr constexpr_map_class bp2idx_map;

    static_assert (bp2idx_map['a'] == 0, "Error.");
    static_assert (bp2idx_map['A'] == 0, "Error.");
    static_assert (bp2idx_map['t'] == 1, "Error.");
    static_assert (bp2idx_map['T'] == 1, "Error.");
    static_assert (bp2idx_map['g'] == 2, "Error.");
    static_assert (bp2idx_map['G'] == 2, "Error.");
    static_assert (bp2idx_map['c'] == 3, "Error.");
    static_assert (bp2idx_map['C'] == 3, "Error.");
    //static_assert (bp2idx_map['C'] == 4, "Error.");//will cause error

    // cannot do this! the evaluation is at runtime
    // template <char bp1, char bp2> struct bp_energy
    // {
    //     enum {val = - factor_pNnm * (BPEnergy::LUTdH[BPEnergy::bp2idx_map[bp1]][BPEnergy::bp2idx_map[bp2]] -
    //                                  (BPEnergy::LUTm[BPEnergy::bp2idx_map[bp1]][BPEnergy::bp2idx_map[bp2]] * EffSaltConc +
    //                                   BPEnergy::LUTdS[BPEnergy::bp2idx_map[bp1]][BPEnergy::bp2idx_map[bp2]] * 0.001 ) * TEMPERATURE);}
    // };
    // static_assert (bp_energy<'a','G'>::val= 0, "Error.");

    inline double lookup_bp_energy(char bp1, char bp2) {
        int idx1 = BPEnergy::bp2idx_map[bp1];
        int idx2 = BPEnergy::bp2idx_map[bp2];
        double energy = BPEnergy::LUTdH[idx1][idx2] - (BPEnergy::LUTm[idx1][idx2] * EFFSALTCONC + BPEnergy::LUTdS[idx1][idx2] * 0.001 ) * TEMPERATURE;
        return - energy * FACTOR_PNNM;//convert to 'pN-nm' unit;
    }
}