#ifndef SCALE_INNER__H
#define SCALE_INNER__H

#include <stdlib.h>
#include "inner.h"

#define SCALE	26


typedef struct {
	uint64_t v;
} fxr;

typedef struct {
	uint64_t v;
} fxr2;

#define FXR(x)   { (x) }
// Define complex data
typedef struct {
	fxr re, im;
} fxc;

extern const fxc GM_TAB[1024];
extern const fxr fxr_p2_tab[];
extern const fxr fxr_inv_sigma[];
extern const fxr fxr_sigma_min[];
extern const fxr fxr_q;
extern const fxr fxr_inverse_of_q;
extern const fxr fxr_inv_2sqrsigma0;
extern const fxr fxr_inv_log2;
extern const fxr fxr_bnorm_max;

extern const fxr fxr_invsqrt2;
extern const fxr fxr_invsqrt8;
extern const fxr fxr_log2;

extern const fxr fxr_zero ;
extern const fxr fxr_sqrt2;	


// x+y (scale = SCALE)
static inline fxr
fxr_add(fxr x, fxr y)
{
	x.v += y.v;
	return x;
}

// x - y with scale = SCALE
static inline fxr
fxr_sub(fxr x, fxr y)
{
	x.v -= y.v;
	return x;
}

// x - y with scale = SCALE * 2
static inline fxr2
fxr2_sub(fxr2 x, fxr2 y)
{
	x.v -= y.v;
	return x;
}


static inline fxr
fxr_mul(fxr fx, fxr fy)
{
	uint64_t x = fx.v, y = fy.v;

	uint64_t s1 = x >> 63;
	x = (x ^ -s1) + s1;
	uint64_t s2 = y >> 63;
	y = (y ^ -s2) + s2;
	s1 ^= s2;

	uint64_t xh = x >> SCALE;
	uint64_t yh = y >> SCALE;
	uint32_t xl = x & ((1 << SCALE) - 1);
	uint32_t yl = y & ((1 << SCALE) - 1);
	uint64_t xhyh = xh * yh;
	uint64_t xlyl = (uint64_t)xl * yl;
	uint64_t xhxl = xh + xl, yhyl = yh + yl;
	uint64_t xhxlyhyl = xhxl * yhyl - xlyl - xhyh;
	uint64_t z = (xlyl >> SCALE) + ((xlyl & ((1 << SCALE) - 1)) > (1 << (SCALE - 1))) + 
		xhxlyhyl + (xhyh << SCALE);
	
	fx.v = (z ^ -s1) + s1;
	return fx;

}


// fx * fy (scale = SCALE * 2)
static inline fxr2
fxr2_mul(fxr2 fx, fxr2 fy)
{
	uint64_t x = fx.v, y = fy.v;

	uint64_t s1 = x >> 63;
	x = (x ^ -s1) + s1;
	uint64_t s2 = y >> 63;
	y = (y ^ -s2) + s2;
	s1 ^= s2;

	uint64_t xh = x >> SCALE;
	uint64_t yh = y >> SCALE;
	uint32_t xl = x & ((1 << SCALE) - 1);
	uint32_t yl = y & ((1 << SCALE) - 1);
	uint64_t xhyh = xh * yh;
	uint64_t xlyl = (uint64_t)xl * yl;
	uint64_t xhxl = xh + xl, yhyl = yh + yl;
	uint64_t xhxlyhyl = xhxl * yhyl - xlyl - xhyh;

	uint64_t c = (xlyl + ((xhxlyhyl & ((((uint64_t)1 << SCALE)) - 1 )) << SCALE)) >
		(uint64_t)1 << (2 * SCALE - 1);

	x = (xhxlyhyl >> SCALE) + xhyh + c;

	fx.v = (x ^ -s1) + s1;
	return fx;

}

// x / y with scale = SCALE
static inline fxr
fxr_div(fxr fx, fxr fy)
{
	
	uint64_t x = fx.v, y = fy.v;
	uint64_t sx = x >> 63;
	x = (x ^ -sx) + sx;
	uint64_t sy = y >> 63;
	y = (y ^ -sy) + sy;
		
	uint64_t q = 0;
	uint64_t temp = (x << (SCALE + 2));
	uint64_t num = x >> (62 - SCALE);

	for (int i = 62; i >= 0; i--)
	{
		uint64_t b = 1 - ((num - y) >> 63);
		q = (q << 1) | b;
		num = ((num - (y & -b)) << 1) | (temp >> 63);
		temp <<= 1;
	}

	uint64_t b = 1 - ((num - y) >> 63);
	q += b;
	q &= 0x7FFFFFFFFFFFFFFF;

	sx ^= sy;
	q = (q ^ -sx) + sx;

	fx.v = q;
	return fx;
}

static inline fxr

// x * x with scale = SCALE
fxr_sqr(fxr fx)
{
	uint64_t x = fx.v;

	uint64_t s1 = x >> 63;
	x = (x ^ -s1) + s1;


	uint64_t xh = x >> SCALE;

	uint32_t xl = x & ((1 << SCALE) - 1);

	uint64_t xhxh = xh * xh;
	uint64_t xlxl = (uint64_t)xl * xl;
	uint64_t xhxl = (xh * xl) << 1;

	uint64_t z = (xlxl >> SCALE) + ((xlxl & ((1 << SCALE) - 1)) > (1 << (SCALE - 1))) +
		xhxl + (xhxh << SCALE);

	fx.v = z;
	return fx;

}

// sqer (x) with scale = SCALE
static inline fxr
fxr_sqrt(fxr fn)
{
	uint64_t x = fn.v, n = x, c = 0, d = (uint64_t)1 << 62;
	while (d > n)
	{
		d >>= 2;
	}
	while (d)
	{

		if (x >= c + d)
		{
			x -= (c + d);
			c = (c >> 1) + d;
		}
		else
		{
			c >>= 1;
		}
		d >>= 2;

	}
	fn.v = c << (SCALE / 2);
	return fn;

}

static inline fxr
fxr_sqrt1(fxr fn)
{
    #define MAX_STEP 62
    uint64_t x = fn.v, n = x, c = 0, d = (uint64_t)1ULL << MAX_STEP;

    uint64_t comp_d, comp_c_cd;
    uint64_t cd;
    for (int step = MAX_STEP; step >= 0; step -= 2)
    {
        comp_d = -(uint64_t)(n >= d);

        cd = comp_d & (c + d);

        comp_c_cd = comp_d & -(uint64_t)(x >= cd);

        x -= (comp_c_cd & cd);

        c = comp_d & (c >> 1);
        c += (comp_c_cd & d);

        d >>= 2;
    }

    fn.v = c << (SCALE / 2);
    return fn;
}


// j with scale = SCALE
static inline fxr
fxr_of(int32_t j)
{
	fxr x;

	x.v = (uint64_t)j << SCALE;
	return x;
}


// x * 2 with scale = SCALE
static inline fxr
fxr_double(fxr x)
{
	x.v <<= 1;
	return x;
}

//-x with scale = SCALE
static inline fxr
fxr_neg(fxr x)
{
	x.v = -x.v;
	return x;
}

// abs (x) with scale
static inline fxr
fxr_abs(fxr x)
{
	x.v -= (x.v << 1) & (uint64_t)((int64_t)x.v >> 63);
	return x;
}

///////////////////////////////////////////
static inline int32_t
fxr_round(fxr x)
{

	x.v += 0x02000000ul;
	return (int32_t)((int64_t)x.v >> SCALE);
}

///////////////////////////////////////
// x / (2 ^ e) with (with scale = SCALE)
static inline fxr
fxr_div2e(fxr x, unsigned n)
{
	int64_t v;


	v = (int64_t)x.v;
	x.v = (uint64_t)((v + (((int64_t)1 << n) >> 1)) >> n);
	return x;
}

// x / 2 (with scale = SCALE)
static inline fxr
fxr_half(fxr x)
{
	x = fxr_div2e(x, 1);
	return x;
}

static inline fxr
fxr_mul2e(fxr x, unsigned n)
{
	x.v <<= n;
	return x;
}

// return  x < y ? 1 : 0
static inline int
fxr_lt(fxr x, fxr y)
{
	return *(int64_t*)&x.v < *(int64_t*)&y.v;
}

static inline int64_t
fxr_trunc(fxr x)
{
	x.v = (uint64_t)(((int64_t)x.v) >> SCALE);
	return (int64_t)x.v;
}

static inline int64_t
fxr_floor(fxr x)
{
	int64_t t = fxr_trunc(x);
	int64_t d = x.v - t;

	t += (d >> 63);

	return t;
}

// fxr -> int64_t
static inline int64_t
fxr_rint(fxr x)
{

	uint64_t v = x.v;
	uint64_t s1 = v >> 63;

	v = (v ^ -s1) + s1;
	v += (1 << (SCALE - 1));
	v >>= SCALE;
	v = (v ^ -s1) + s1;
	x.v = v;

	return v;
}


// Complex type
#define FXC(re, im)   { FXR(re), FXR(im) }

// x.re + y.re, x.im + y.im (with scale = SCALE)
static inline fxc
fxc_add(fxc x, fxc y)
{
	x.re = fxr_add(x.re, y.re);
	x.im = fxr_add(x.im, y.im);
	return x;
}

// x.re - y.re, x.im - y.im (with scale = SCALE)
static inline fxc
fxc_sub(fxc x, fxc y)
{
	x.re = fxr_sub(x.re, y.re);
	x.im = fxr_sub(x.im, y.im);
	return x;
}


// x.re / 2, y.re / 2 (with scale = SCALE)
static inline fxc
fxc_half(fxc x)
{
	x.re = fxr_div2e(x.re, 1);
	x.im = fxr_div2e(x.im, 1);
	return x;
}


// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// return x * y (with scale = SCALE), x, y - complex
static inline fxc
fxc_mul(fxc x, fxc y)
{
	/*
	 * We are computing r = (a + i*b)*(c + i*d) with:
	 *   z0 = a*c
	 *   z1 = b*d
	 *   z2 = (a + b)*(c + d)
	 *   r = (z0 - z1) + i*(z2 - (z0 + z1))
	 * Since the intermediate values are truncated to our precision,
	 * the imaginary value of r _may_ be slightly different from
	 * a*d + b*c (if we had calculated it directly). For full
	 * reproducibility, all implementations should use the formulas
	 * above.
	 */
	fxr z0 = fxr_mul(x.re, y.re);
	fxr z1 = fxr_mul(x.im, y.im);
	fxr z2 = fxr_mul(fxr_add(x.re, x.im), fxr_add(y.re, y.im));
	fxc z;
	z.re = fxr_sub(z0, z1);
	z.im = fxr_sub(z2, fxr_add(z0, z1));
	return z;
}

static inline fxc
fxc_div1(fxc x, fxc y)
{

	fxc temp1 = { (int64_t)(x.re.v) >> 1, (int64_t)x.im.v >> 1 }, 
		temp2 = { (int64_t)(y.re.v) >> 1, (-(int64_t)(y.im.v)) >> 1 },
		temp3 = { 0,0 };
	
	fxr fxr_temp1, fxr_temp2;
	fxr_temp1 = fxr_sqr(temp2.im);
	fxr_temp2 = fxr_sqr(temp2.re);
	fxr_temp1 = fxr_add(fxr_temp1, fxr_temp2);
	temp1.re = fxr_div(temp1.re, fxr_temp1);
	temp1.im = fxr_div(temp1.im, fxr_temp1);
	temp3 = fxc_mul(temp1, temp2);
	return temp3;
}

static inline fxc
fxc_conj(fxc x)
{
	x.im = fxr_neg(x.im);
	return x;
}



// functions with fxr
void
expand_privkey(fxr* restrict expanded_key,
	const int8_t* f, const int8_t* g,
	const int8_t* F, const int8_t* G,
	unsigned logn, uint8_t* restrict tmp);

fxr2 fxr2_exp_(fxr2 x);
// vect with fxr

#define vect_set   Zf(vect_set)
void vect_set(unsigned logn, fxr* d, const int8_t* f);

#define vect_neg   Zf(vect_neg)
void vect_neg(unsigned logn, fxr* d);

#define vect_norm_fft   Zf(vect_norm_fft)
void
vect_norm_fft(unsigned logn, fxr* restrict d,
	const fxr* restrict a, const fxr* restrict b);

#define vect_FFT   Zf(vect_FFT)
void 
vect_FFT(unsigned logn, fxr* f);

#define vect_iFFT   Zf(vect_iFFT)
void 
vect_iFFT(unsigned logn, fxr* f);

#define vect_invnorm_fft   Zf(vect_invnorm_fft)
void 
vect_invnorm_fft(unsigned logn, fxr* restrict d,
	const fxr* restrict a, const fxr* restrict b, unsigned e);

#define vect_adj_fft   Zf(vect_adj_fft)
void 
vect_adj_fft(unsigned logn, fxr* a);

#define vect_mul_realconst   Zf(vect_mul_realconst)
void 
vect_mul_realconst(unsigned logn, fxr* a, fxr c);

#define vect_mul_autoadj_fft   Zf(vect_mul_autoadj_fft)
void 
vect_mul_autoadj_fft(unsigned logn, fxr* restrict a, const fxr* restrict b);

#define vect_div_autoadj_fft   Zf(vect_div_autoadj_fft)
void
vect_div_autoadj_fft(unsigned logn, fxr* restrict a, const fxr* restrict b);

void
vect_mul2e(unsigned logn, fxr* a, unsigned e);
void
vect_mul_fft(unsigned logn, fxr* restrict a, const fxr* restrict b);
void
vect_add(unsigned logn, fxr* restrict a, const fxr* restrict b);
void
vect_mulconst(unsigned logn, fxr* a, fxr x);
void
vect_sub(
	unsigned logn, fxr* restrict a, const fxr* restrict b);

#endif

