// CPU version of LUT calculation, only single point
#include <iostream>
#include "../include/constants.h"

#define NTHREAD 16
// simulation precisions, etc.
#define VALIDMAXFORCE       100.0f
#define VALIDMINFORCE       0.01f
#define TOR_BINARY_SEARCH   0.00000001f
#define VERYLARGENUMBER     1.0e10f
//When dE > this threshold, don't calculate exp(-e/kT) and set probability to 0;
#define ENERGY_THRESHOLD    50.0f
// Taylor expansion of the Langevin function and integration of Langevin function
double Langevin(double x) {
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
double Langevin_integ(double x) {
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

//==============================WLC high force, Odijk 1995 macromolecules===========================
double alpha2phi_Odijk95(double alpha, double k0_eff) { 
    // if (alpha < 0.25) {
    //     return 0.0;// undefined at alpha == 0, just give it a small value
    //     //I can do this because this is force, energy must be calculated correctly!!
    // }
    return 1.0 - 0.5 / sqrt(alpha) + alpha / k0_eff; 
}
double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
    return alpha - sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
}
double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
    return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
}
// ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
//Modified version specific for ssDNA force region, and keeps accuracy
//For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
//My homemade Langevin_integ function should be accurate enough in this region.
double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
    return Langevin(2.0 * alpha) + alpha / k0_eff;
}
double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
    return 0.5 * Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
}
double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
    //integ actually starts from 1, but it's OK since it is for partition function calculation
    return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
}
// calculate the force and energy at extension when j pairs are unzipped
double ext_at_this_force(double force, double j) {//func used to find force so the system total extension = ext
    return  force/PILLARSTIFFNESS + 
            ARMLENGTH * L0DS * alpha2phi_Odijk95(force * LPDS / KT, KDS * LPDS / KT) + 
            2.0 * j * L0SS * alpha2phi_Smith95_m(force * LPSS / KT, KSS * LPSS / KT);//increasing function with force
}
void get_force(int j, int z, float * force) {
    //j == length of unzipped trunk, ext total extension
    //simple binary search to get force so calc_z_diff(force) == 0
    //force function must be monotonic
    double f1 = VALIDMINFORCE;
    double f2 = VALIDMAXFORCE;
    std::cout << "z:" << z << ", j:" << j << std::endl;

    double y1 = z - ext_at_this_force(f1, j);
    double y2 = z - ext_at_this_force(f2, j);
    std::cout << "y1:" << y1 << ", y2:" << y2 << std::endl;

    if (y1 * y2 >= 0) {
        if (y1 < 0){
            *force = VALIDMAXFORCE * 2.0;//force is too large
            return;
        } else {
            *force =  VALIDMINFORCE / 2.0;//force is too small
            return;
        }
    }

    int cnt = 0;//in case the root is not found
    double fm, ym;
    while (++cnt < 10000) {
        
        fm = (f1 + f2) * 0.5;
        ym = z - ext_at_this_force(fm, j);
        
        if (std::abs(ym) <= TOR_BINARY_SEARCH) {
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

int main() {
    int bp_trunk = 10000;
    int X_DIM = (bp_trunk + NTHREAD - 1);
    X_DIM -= X_DIM % NTHREAD;
    int Y_DIM = (int)((bp_trunk + ARMLENGTH) * L0SS * 2) + NTHREAD - 1;
    Y_DIM -= Y_DIM % NTHREAD;
    
    float *force = new float [X_DIM * Y_DIM];
    float *energy = new float [X_DIM * Y_DIM];

    //get_LUT(X_DIM, Y_DIM, force, energy);

    printf("force[123] = %f\n", force[123]);

    float temp;
    get_force(100, 800, &temp);
    printf("force %f\n", temp);

    delete[] force, energy;
    return 0;
}