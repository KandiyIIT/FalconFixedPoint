#ifndef SIGN__H
#define SIGN__H

#include "poly.h"
#include "shake256.h"

extern const uint32_t l2bound[];

void
hash_to_point_vartime(
	inner_shake256_context* sc,
	uint16_t* x, unsigned logn);

int verify_raw(const uint16_t* c0, const int16_t* s2,
	const uint16_t* h, unsigned logn, uint8_t* tmp);

void
sign_tree(int16_t* sig, inner_shake256_context* rng,
	const fxr* restrict expanded_key,
	const uint16_t* hm, unsigned logn, uint8_t* tmp);

#endif
