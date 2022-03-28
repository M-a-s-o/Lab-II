//  g++ -c RNG.cc -o RNG.o
// ᓚᘏᗢ
#include "RNG.h"
#include <string.h>
#include <bitset>
#include <iostream>

/* --- Linear congruential generator --- */

LCG::LCG (unsigned long seed):
    _seed (seed)
    {}

LCG::LCG ():
    _seed (1)
    {}

LCG::~LCG () {}

unsigned long LCG::random (unsigned long seed) {
    unsigned int M = 2147483647, A = 214013, C = 2531011;
    static unsigned long x = seed;
    x = (A*x+C)%M;
    return x;
}

unsigned long LCG::random () {
    return random(_seed);
}


/* --- RNG algorithms --- */

double RNG_i (int min, int max) {
    return min + rand()%(max-min);
}

double RNG_d (double min, double max) {
    if (min == 0 && max == 1)
        return rand() / static_cast<double> (RAND_MAX);
    else
        return min + (max-min)*rand() / static_cast<double> (RAND_MAX);
}

inline double DoubleFromBits(uint64_t i) {
	return (i >> 11) * 0x1.0p-53;
}

double RNG_i_xor (int min, int max) {
    return min + XoshiroNext()%(max-min);
}

double RNG_d_xor () {
    return DoubleFromBits(XoshiroNext());
}

double RTAC (double f (double), double xMin, double xMax, double yMax) {
    // Try and catch
    double x, y;
    do {
      x = RNG_d (xMin, xMax);
      y = RNG_d (0, yMax);
    } while (y > f(x));
    return x;
}

double RTAC_xor (double f (double), double xMin, double xMax, double yMax) {
    // Try and catch
    double x, y;
    do {
        x = xMin + (xMax-xMin)*RNG_d_xor();
        y = yMax*RNG_d_xor ();
    } while (y > f(x));
    return x;
}

double RInv (double f (double)) {
    double x = RNG_d (0, 1);
    return f(x);
}

double RInv_xor (double f (double)) {
    double x = RNG_d_xor();
    return f(x);
}

double Med(int num, double min, double max) {
    double sum = 0;
    for (int i = 0; i < num; ++i) {
        sum += RNG_d(min, max);
    }
    return sum/(static_cast<double> (num));
}

double RNGExp (double lambda)   {
    if (lambda == 0.)
        return 0;
    return 0.-log(1-RNG_d(0, 1))/lambda;
}

double RNGPoisson (double lambda) {
    /* https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
    * Poisson generator based upon the inversion by sequential search
    * Devroye, Luc (1986), "Discrete Univariate Distributions" (PDF), Non-Uniform Random Variate Generation, New York, NJ, USA: Springer-Verlag, pp. 485–553, doi:10.1007/978-1-4613-8643-8_10, ISBN 978-1-4613-8645-2
    * http://luc.devroye.org/rnbookindex.html
    * Chapter 10, page 505
    */

    int x = 0;
    double s, p, u = RNG_d(0,1);
    s = p = exp(0.-lambda);
    while (u > s) {
        ++x;
        p *= (lambda/(static_cast<double> (x)));
        s += p;
    }
    return x;
    /*double t_tot = RNGExp(1) ;
    int N_evt = 0 ;
    while (t_tot < lambda)
        {
            ++N_evt ;
            t_tot += RNGExp(1) ;
        }
    return N_evt ;*/
}

double RNGGauss (double mean, double std_dev, int sums) {
    // Generazione tramite il teorema centrale del limite
    int64_t sqrt3 = 0x3ffbb67ae8584caa;
    double sqrt3_d = *(double *) &sqrt3; // Quake III docet
    double dev = sqrt3_d*sqrt(sums)*std_dev;
    return Med(sums, mean-dev, mean+dev);
}

/*
double RNGGauss (double mean, double std_dev, int sums) {
    // Marsaglia polar method
    double U = RNG_d(-1,1), V = RNG_d(-1,1);
    double S = U*U+V*V;
    while (S >= 1) {
        U = RNG_d(-1,1), V = RNG_d(-1,1);
        S = U*U+V*V;
    }
    return U*sqrt(-2*log(S)/S);
}*/

double RNGGauss_ziggurat (double mean, double std_dev) {
    return Sample(mean, std_dev);
}

void hex_binary_print_double (double var) {
    char result[sizeof(double)];

    std::cout << "Decimal: " << var << std::endl;
    uint64_t u;
    memcpy(&u, &var, sizeof(double));
    std::cout << "Hex: "<< std::hex << u << std::endl;
    std::bitset<64> x(u);
    std::cout << "Binary: "<< x << std::endl;
}