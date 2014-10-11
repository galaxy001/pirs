#ifndef _MT19937_64_H
#define _MT19937_64_H

#include "inttypes.h"

#define MERSENNE_TWISTER_STATE_LEN 312

extern void init_genrand64(uint64_t seed, uint64_t mt[MERSENNE_TWISTER_STATE_LEN]);
extern uint64_t genrand64_int64(uint64_t mt[MERSENNE_TWISTER_STATE_LEN]);
extern double genrand64_real1(uint64_t mt[MERSENNE_TWISTER_STATE_LEN]);
extern double genrand64_real2(uint64_t mt[MERSENNE_TWISTER_STATE_LEN]);
extern double genrand64_real3(uint64_t mt[MERSENNE_TWISTER_STATE_LEN]);

#endif /* _MT19937_64_H */
