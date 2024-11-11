#include "inner.h"
// Thomas Pornin. https://github.com/pornin/ntrugen
void
mp_NTT(unsigned logn, uint32_t* restrict a, const uint32_t* restrict gm,
	uint32_t p, uint32_t p0i)
{
	if (logn == 0) {
		return;
	}
	size_t t = (size_t)1 << logn;
	for (unsigned lm = 0; lm < logn; lm++) {
		size_t m = (size_t)1 << lm;
		size_t ht = t >> 1;
		size_t v0 = 0;
		for (size_t u = 0; u < m; u++) {
			uint32_t s = gm[u + m];
			for (size_t v = 0; v < ht; v++) {
				size_t k1 = v0 + v;
				size_t k2 = k1 + ht;
				uint32_t x1 = a[k1];
				uint32_t x2 = mp_montymul(a[k2], s, p, p0i);
				a[k1] = mp_add(x1, x2, p);
				a[k2] = mp_sub(x1, x2, p);
			}
			v0 += t;
		}
		t = ht;
	}

}

// Thomas Pornin. https://github.com/pornin/ntrugen
void
mp_iNTT(unsigned logn, uint32_t* restrict a, const uint32_t* restrict igm,
	uint32_t p, uint32_t p0i)
{
	if (logn == 0) {
		return;
	}

	size_t t = 1;
	for (unsigned lm = 0; lm < logn; lm++) {
		size_t hm = (size_t)1 << (logn - 1 - lm);
		size_t dt = t << 1;
		size_t v0 = 0;
		for (size_t u = 0; u < hm; u++) {
			uint32_t s = igm[u + hm];
			for (size_t v = 0; v < t; v++) {
				size_t k1 = v0 + v;
				size_t k2 = k1 + t;
				uint32_t x1 = a[k1];
				uint32_t x2 = a[k2];
				a[k1] = mp_half(mp_add(x1, x2, p), p);
				a[k2] = mp_montymul(
					mp_sub(x1, x2, p), s, p, p0i);
			}
			v0 += dt;
		}
		t = dt;
	}

}