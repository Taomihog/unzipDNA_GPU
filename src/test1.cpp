#include <iostream>
#include "constants.h"
int main() {

    // test the bp energy
    char ch1 = 'a';
    char ch2 = 't';
    printf("some result: %f.\n", BPEnergy::lookup_bp_energy(ch1,ch2));


    return 0;
}