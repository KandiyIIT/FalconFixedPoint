
#include <stdio.h>
#include <memory.h>

#include "inner.h"
#include "ntt.h"
#include "scale_inner.h"
#define RNG_CONTEXT   inner_shake256_context
// use https://github.com/pornin/ntrugen

// Thomas Pornin. https://github.com/pornin/ntrugen
// Changes
//   Functiom keygen
//		poly_is_invertible do not use; 
//		added public key calculation (function compute_public)
//		added counting of errors of different types (Debug mode)
const ntru_profile SOLVE_Falcon_512 = {
	12289,
	9, 9,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	13,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

// Thomas Pornin. https://github.com/pornin/ntrugen
const ntru_profile SOLVE_Falcon_1024 = {
	12289,
	10, 10,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	11,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

// Thomas Pornin. https://github.com/pornin/ntrugen
/* Falcon, q = 12289, n = 512 -> kmax = 17 */
const uint16_t gauss_Falcon_512[] = {
	17,
		1,     4,    11,    28,    65,   146,   308,   615,
	 1164,  2083,  3535,  5692,  8706, 12669, 17574, 23285,
	29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000,
	63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
	65531, 65534
};

// Thomas Pornin. https://github.com/pornin/ntrugen
/* Falcon, q = 12289, n = 1024 -> kmax = 12 */
const uint16_t gauss_Falcon_1024[] = {
	12,
		2,     8,    28,    94,   280,   742,  1761,  3753,
	 7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
	61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533
};

// Thomas Pornin. https://github.com/pornin/ntrugen

void
gauss_sample_poly(unsigned logn, int8_t* f,
	uint16_t* tab, inner_shake256_context* rng)
{
	size_t n = (size_t)1 << logn;
	size_t kmax = tab[0];
	
	for (;;) {
		uint32_t parity = 0;
		for (size_t j = 0; j < n; j++) {
			uint32_t v = -(uint32_t)kmax;
			//uint32_t x = prng_buffer_next_u16(&pb);
			uint32_t x = 0;
			inner_shake256_extract(rng, (uint8_t*)&x, 2);
			for (size_t k = 1; k <= (kmax << 1); k++) {
				v += ((uint32_t)tab[k] - x) >> 31;
			}
			f[j] = (int8_t) * (int32_t*)&v;
			parity ^= v;
		}
		if ((parity & 1) != 0) {
			return;
		}
	}
}

// Thomas Pornin. https://github.com/pornin/ntrugen
uint32_t
poly_sqnorm(unsigned logn, const int8_t* f)
{
	size_t n = (size_t)1 << logn;

	uint32_t s = 0;
	for (size_t u = 0; u < n; u++) {
		int32_t x = f[u];
		s += (uint32_t)(x * x);
	}
	return s;
}


// Thomas Pornin. https://github.com/pornin/ntrugen
// Added public key generation (compute_public function)
// Added counting of different types of errors
// To work with fixed point data, use the scale SCALE = 26
// calculation the orthogonalized vector and solve_NTRU function)

int 
Zf(keygen)(inner_shake256_context* rng,
	int8_t* f, int8_t* g, int8_t* F, int8_t* G, uint16_t* h,
	unsigned logn, uint8_t* tmp, size_t tmp_len)
{


	size_t n;
	
	RNG_CONTEXT* rc;
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t* tt32 = (void*)utmp2;


	n = MKN(logn);

	rc = rng;
	const ntru_profile* profs[] = { &SOLVE_Falcon_512, &SOLVE_Falcon_1024 }, * prof;
	const uint16_t* gaus_tabls [] = { gauss_Falcon_512, gauss_Falcon_1024 };
	uint16_t* gaus_tabl;
	int res = (logn != 9) && (logn != 10);

	if (res == 0)
	{
		prof = profs[logn - 9];
		gaus_tabl = gaus_tabls[logn - 9];
#ifdef _DEBUG
		all_error = 0, fg_error = 0, fg_inv_error = 0, fp_fg_error = 0, 
		gcd_error = 0, reduce_error = 0, limit_error_f = 0, limit_error_g = 0;
#endif
		for (;;)
		{
			gauss_sample_poly(logn, f, gaus_tabl, rng);
			gauss_sample_poly(logn, g, gaus_tabl, rng);


			if ((poly_sqnorm(logn, f) + poly_sqnorm(logn, g)) >= 16823) {
#ifdef _DEBUG
				++all_error;
				
				++fg_error;
				
#endif

				continue;
			}
			
			/*
			* For Falcon, we need to check that the orthogonalized
			* vector also has an acceptable norm.
			*/
			size_t n = (size_t)1 << logn;
			fxr* rt1 = (fxr*)tt32;
			fxr* rt2 = rt1 + n;
			fxr* rt3 = rt2 + n;
			vect_set(logn, rt1, f);
			vect_set(logn, rt2, g);
			vect_FFT(logn, rt1);
			vect_FFT(logn, rt2);
			vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
			vect_adj_fft(logn, rt1);
			vect_adj_fft(logn, rt2);
			vect_mul_realconst(logn, rt1, fxr_q);
			vect_mul_realconst(logn, rt2, fxr_q);
			vect_mul_autoadj_fft(logn, rt1, rt3);
			vect_mul_autoadj_fft(logn, rt2, rt3);
			vect_iFFT(logn, rt1);
			vect_iFFT(logn, rt2);
			fxr sn = fxr_zero;
			for (size_t u = 0; u < n; u++) {
				sn = fxr_add(sn,
					fxr_add(fxr_sqr(rt1[u]), fxr_sqr(rt2[u])));
			}

			if (!fxr_lt(sn, fxr_bnorm_max))
			{
#ifdef _DEBUG
				++all_error;
				++fp_fg_error;
#endif
				
				continue;
			}
			
			if (!Zf(compute_public)(h, f, g, logn, tmp)) {
#ifdef _DEBUG
				++all_error;
				++public_error;
#endif

				continue;
			}
			
			res = solve_NTRU(prof, logn, f, g, tt32);
			if (res == SOLVE_OK)
			{
				
				break;
			}
	
		}
		memcpy(F, tt32, n);
		memcpy(G, (uint8_t*)tt32 + n, n);
	}
	return res;
}
