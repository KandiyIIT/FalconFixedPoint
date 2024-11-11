#ifndef POLY__H
#define POLY__H
#include "scale_inner.h"
#define poly_big_to_fixed   Zf(poly_big_to_fixed)
void 
poly_big_to_fixed(unsigned logn, fxr* restrict d,
	const uint32_t* restrict f, size_t len, uint32_t sc);

#define poly_mulselfadj_fft   Zf(poly_mulselfadj_fft)
void 
poly_mulselfadj_fft(fxr* a, unsigned logn);

#define poly_add   Zf(poly_add)
void 
poly_add(fxr* restrict a, const fxr* restrict b, unsigned logn);


#define poly_muladj_fft   Zf(poly_muladj_fft)
void 
poly_muladj_fft(fxr* restrict a, const fxr* restrict b, unsigned logn);

#define poly_LDLmv_fft   Zf(poly_LDLmv_fft)
void 
poly_LDLmv_fft(
	fxr* restrict d11, fxr* restrict l10,
	const fxr* restrict g00, const fxr* restrict g01,
	const fxr* restrict g11, unsigned logn);

#define poly_sub_scaled   Zf(poly_sub_scaled)
int poly_sub_scaled(unsigned logn,
	uint32_t* restrict F, size_t Flen,
	const uint32_t* restrict f, size_t flen,
	//const int64_t* restrict k, 
	const int32_t* restrict k,
	uint32_t sc);
#define poly_split1_fft   Zf(poly_split1_fft)
void
poly_split1_fft(
	fxr* restrict f0, fxr* restrict f1,
	const fxr* restrict f, unsigned logn);

#define poly_merge1_fft   Zf(poly_merge1_fft)
void
poly_merge1_fft(
	fxr* restrict f,
	const fxr* restrict f0, const fxr* restrict f1, unsigned logn);


#endif
