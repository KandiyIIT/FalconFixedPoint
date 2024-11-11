#ifndef FALCON__H
#define FALCON__H

#include "scale_inner.h"
#include "inner.h"

void
Zf(keygen)(inner_shake256_context* rng,
	int8_t* f, int8_t* g, int8_t* F, int8_t* G, uint16_t* h,
	unsigned logn, uint8_t* tmp, size_t tmp_len);

#endif // FALCON__H
