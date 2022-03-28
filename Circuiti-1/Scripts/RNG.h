// ᓚᘏᗢ
#ifndef RNG
#define RNG

#include <stdlib.h>
#include <time.h>
#include <math.h>

extern "C" {
    #include "xoshiro256plus.h"
    #include "gauss.h"
}

/* --- Linear congruential generator --- */
class LCG {
    public:
        LCG (unsigned long seed = 1);
        LCG ();
        ~LCG ();
        unsigned long random (unsigned long seed);
        unsigned long random ();

    private:
        unsigned long _seed;
};

/* --- RNG algorithms --- */
double RNG_i (int min, int max);
double RNG_d (double min = 0, double max = 1);
inline double DoubleFromBits(uint64_t i);
double RNG_i_xor (int min, int max);
double RNG_d_xor ();
double RTAC (double f (double), double xMin, double xMax, double yMax = 1/sqrt(2*M_PI));
double RTAC_xor (double f (double), double xMin, double xMax, double yMax = 1/sqrt(2*M_PI));
double RInv (double f (double));
double RInv_xor (double f (double));
double Med (int num = 10, double min = 0, double max = 1);
double RNGExp (double lambda);
double RNGPoisson (double lambda);
double RNGGauss (double mean = 0, double std_dev = 1, int sums = 10);
double RNGGauss_ziggurat (double mean = 0, double std_dev = 1);
void hex_binary_print_double (double var);

#endif