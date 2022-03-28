//	gcc -c xoshiro256plus.c -O3 -o xoshiro256plus.o
//  ᓚᘏᗢ
/*
	This algorithm has been adapted from Blackman's and Vigna's work:
	https://prng.di.unimi.it/xoshiro256plus.c
	https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf

	As suggested, this xoshiro256+ implementation directly uses Vigna's Splitmix64 generator
	to seed values. See: https://xorshift.di.unimi.it/splitmix64.c

	For further information, see https://prng.di.unimi.it/

	This algorithm is not cryptographically secure.
*/
#include "xoshiro256plus.h"
#include <time.h>

// Splitmix64 state that can be seeded with any value.
uint64_t x;

// Xoshiro internal state.
static uint64_t s[4];


/* --- Splitmix64 algorithm --- */

uint64_t SplitNext () {
	uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
	z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
	return z ^ (z >> 31);
}


/* --- xoshiro256+ algorithm --- */

static inline uint64_t rotl(const uint64_t y, int k) {
	return (y << k) | (y >> (64 - k));
}

uint64_t XoshiroNext(void) {
	const uint64_t result = s[0] + s[3];

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

void SetupXoshiro () {
	x = time(NULL);
	for (int i = 0; i < 4; ++i)
		s[i] = SplitNext();
	return;
}