#ifndef INNER__H
#define INNER__H

#include <memory.h>
#include <inttypes.h>

#if defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4146 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 4334 )

#ifndef restrict
#define restrict   __restrict
#endif
#endif

#define	MKN(logn)	(1 << (logn))
#ifndef FALCON_PREFIX
#define FALCON_PREFIX   falcon_inner
#endif
#define Zf(name)             Zf_(FALCON_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name  

extern int all_error, fg_error , fp_fg_error , fg_inv_error, public_error, 
gcd_error , reduce_error, limit_error_f, limit_error_g;
extern int _all_error, _fg_error, _fp_fg_error, _fg_inv_error, _public_error,
_gcd_error, _reduce_error, _limit_error_f, _limit_error_g;

/* Error code: no error (so far) */
#define SOLVE_OK           0

/* Error code: GCD(Res(f,X^n+1), Res(g,X^n+1)) != 1 */
#define SOLVE_ERR_GCD      -1

/* Error code: reduction error (NTRU equation no longer fulfilled) */
#define SOLVE_ERR_REDUCE   -2

/* Error code: output (F,G) coefficients are off-limits */
#define SOLVE_ERR_LIMIT    -3

#define FALCON_EXPANDEDKEY_SIZE(logn) \
	(((8u * (logn) + 40) << (logn)) + 8)

typedef struct {
	union {
		uint64_t A[25];
		uint8_t dbuf[200];
	} st;
	uint64_t dptr;
} inner_shake256_context;

typedef struct {
	union {
		uint8_t d[512]; /* MUST be 512, exactly */
		uint64_t dummy_u64;
	} buf;
	size_t ptr;
	union {
		uint8_t d[256];
		uint64_t dummy_u64;
	} state;
	int type;
} prng;

typedef struct {
	uint32_t q;
	unsigned min_logn, max_logn;
	uint16_t max_bl_small[11];
	uint16_t max_bl_large[10];
	uint16_t word_win[10];
	uint32_t reduce_bits;
	uint8_t coeff_FG_limit[11];
	uint16_t min_save_fg[11];
} ntru_profile;


#define inner_shake256_init      Zf(i_shake256_init)
#define inner_shake256_inject    Zf(i_shake256_inject)
#define inner_shake256_flip      Zf(i_shake256_flip)
#define inner_shake256_extract   Zf(i_shake256_extract)

void Zf(i_shake256_init)(
	inner_shake256_context* sc);
void Zf(i_shake256_inject)(
	inner_shake256_context* sc, const uint8_t* in, size_t len);
void Zf(i_shake256_flip)(
	inner_shake256_context* sc);
void Zf(i_shake256_extract)(
	inner_shake256_context* sc, uint8_t* out, size_t len);
void Zf(prng_refill)(prng* p);
int Zf(compute_public)(uint16_t* h,
	const int8_t* f, const int8_t* g, unsigned logn, uint8_t* tmp);

#define poly_is_invertible   Zf(poly_is_invertible)
int
poly_is_invertible(unsigned logn, const int8_t* restrict f,
	uint32_t p, uint32_t p0i, uint32_t s,
	uint32_t r, uint32_t rm, unsigned rs, uint32_t* restrict tmp);

#define poly_mp_set_small   Zf(poly_mp_set_small)
void poly_mp_set_small(unsigned logn, uint32_t* restrict d,
	const int8_t* restrict f, uint32_t p);
#define mp_mkigm   Zf(mp_mkigm)
void mp_mkigm(unsigned logn, uint32_t* restrict igm,
	uint32_t ig, uint32_t p, uint32_t p0i);




#define mp_mkgm   Zf(mp_mkgm)
void mp_mkgm(unsigned logn, uint32_t* restrict gm,
	uint32_t g, uint32_t p, uint32_t p0i);

#define mp_NTT   Zf(mp_NTT)
void mp_NTT(unsigned logn, uint32_t* restrict a, const uint32_t* restrict gm,
	uint32_t p, uint32_t p0i);

#define mp_iNTT   Zf(mp_iNTT)
void mp_iNTT(unsigned logn, uint32_t* restrict a, const uint32_t* restrict igm,
	uint32_t p, uint32_t p0i);

#define zint_mod_small_unsigned   Zf(zint_mod_small_unsigned)
uint32_t zint_mod_small_unsigned(const uint32_t* d, size_t len, size_t stride,
	uint32_t p, uint32_t p0i, uint32_t R2);
#define zint_mul_small   Zf(zint_mul_small)
uint32_t zint_mul_small(uint32_t* m, size_t len, uint32_t x);

#define zint_norm_zero   Zf(zint_norm_zero)
void zint_norm_zero(uint32_t* restrict x, size_t len, size_t xstride,
	const uint32_t* restrict p);




#define zint_rebuild_CRT   Zf(rebuild_CRT)
void zint_rebuild_CRT(uint32_t* restrict xx, size_t xlen, size_t n,
	size_t num_sets, int normalize_signed, uint32_t* restrict tmp);

#define zint_bezout   Zf(zint_bezout)
int zint_bezout(uint32_t* restrict u, uint32_t* restrict v,
	const uint32_t* restrict x, const uint32_t* restrict y,
	size_t len, uint32_t* restrict tmp);
#define zint_negate   Zf(zint_negate)
void zint_negate(uint32_t* a, size_t len, uint32_t ctl);

#define poly_max_bitlength   Zf(poly_max_bitlength)
uint32_t poly_max_bitlength(unsigned logn, const uint32_t* f, size_t flen);

#define mp_mkgmigm   Zf(mp_mkgmigm)
void mp_mkgmigm(unsigned logn, uint32_t* restrict gm, uint32_t* restrict igm,
	uint32_t g, uint32_t ig, uint32_t p, uint32_t p0i);

#define DIVREM31(q, r, x)  { \
		uint32_t divrem31_q, divrem31_x; \
		divrem31_x = (x); \
		divrem31_q = (uint32_t)(divrem31_x * (uint32_t)67651) >> 21; \
		(q) = divrem31_q; \
		(r) = divrem31_x - 31 * divrem31_q; \
	} while (0)


#define poly_sub_kfg_scaled_depth1   Zf(poly_sub_kfg_scaled_depth1)
void poly_sub_kfg_scaled_depth1(unsigned logn_top,
	uint32_t* restrict F, uint32_t* restrict G, size_t FGlen,
	int32_t* restrict k,
	uint32_t sc,
	const int8_t* restrict f, const int8_t* restrict g,
	uint32_t* restrict tmp);

#define poly_sub_scaled_ntt   Zf(poly_sub_scaled_ntt)
void
poly_sub_scaled_ntt(unsigned logn, uint32_t * restrict F, size_t Flen,
	const uint32_t * restrict f, size_t flen,
	const int32_t * restrict k,
	uint32_t sc, uint32_t * restrict tmp);

#define zint_sub_scaled   Zf(zint_sub_scaled)
void zint_sub_scaled(uint32_t * restrict x, size_t xlen,
	const uint32_t * restrict y, size_t ylen, size_t stride,
	uint32_t sch, uint32_t scl);
#define poly_mp_set   Zf(poly_mp_set)
void poly_mp_set(unsigned logn, uint32_t * f, uint32_t p);

/*
 * Add k*(2^sc)*y to x. The result is assumed to fit in the array of
 * size xlen (truncation is applied if necessary).
 * Scale factor sc is provided as sch and scl, such that:
 *    sch = sc / 31
 *    scl = sc % 31  (in the 0..30 range)
 * xlen MUST NOT be lower than ylen; however, it is allowed that
 * xlen is greater than ylen.
 *
 * x[] and y[] are both signed integers, using two's complement for
 * negative values. They both use the same stride ('stride' parameter).
 */
#define zint_add_scaled_mul_small   Zf(zint_add_scaled_mul_small)
void zint_add_scaled_mul_small(uint32_t * restrict x, size_t xlen,
	const uint32_t * restrict y, size_t ylen, size_t stride,
	int32_t k, uint32_t sch, uint32_t scl);
#define poly_mp_norm   Zf(poly_mp_norm)
void poly_mp_norm(unsigned logn, uint32_t * f, uint32_t p);

#define poly_big_to_small   Zf(poly_big_to_small)
int poly_big_to_small(unsigned logn, int8_t * restrict d,
	const uint32_t * restrict s, int lim);

#define lzcnt_nonzero   lzcnt

static inline uint32_t
tbmask(uint32_t x)
{
	return (uint32_t)(*(int32_t*)&x >> 31);
}

static inline unsigned
lzcnt(uint32_t x)
{
#if NTRUGEN_AVX2
	/*
	 * All AVX2-capable CPUs have lzcnt.
	 */
	return _lzcnt_u32(x);
#else // NTRUGEN_AVX2
	uint32_t m = tbmask((x >> 16) - 1);
	uint32_t s = m & 16;
	x = (x >> 16) ^ (m & (x ^ (x >> 16)));
	m = tbmask((x >> 8) - 1);
	s |= m & 8;
	x = (x >> 8) ^ (m & (x ^ (x >> 8)));
	m = tbmask((x >> 4) - 1);
	s |= m & 4;
	x = (x >> 4) ^ (m & (x ^ (x >> 4)));
	m = tbmask((x >> 2) - 1);
	s |= m & 2;
	x = (x >> 2) ^ (m & (x ^ (x >> 2)));

	/*
	 * At this point, x fits on 2 bits. Number of leading zeros is
	 * then:
	 *    x = 0   -> 2
	 *    x = 1   -> 1
	 *    x = 2   -> 0
	 *    x = 3   -> 0
	 */
	return (unsigned)(s + ((2 - x) & tbmask(x - 3)));
#endif // NTRUGEN_AVX2
}


static inline uint32_t
mp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint64_t z = (uint64_t)a * (uint64_t)b;
	uint32_t w = (uint32_t)z * p0i;
	uint32_t d = (uint32_t)((z + (uint64_t)w * (uint64_t)p) >> 32) - p;
	return d + (p & tbmask(d));
}

/*
 * Compute R = 2^32 mod p.
 */
static inline uint32_t
mp_R(uint32_t p)
{
	/*
	 * Since 2*p < 2^32 < 3*p, we just subtract 2*p from 2^32.
	 */
	return -(p << 1);
}


static inline uint32_t
mp_set(int32_t v, uint32_t p)
{
	return (uint32_t)(v + (p & (v >> 31)));
}

/*
 * Addition modulo p.
 */
static inline uint32_t
mp_add(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a + b - p;
	return d + (p & tbmask(d));
}
static inline uint32_t
mp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a - b;
	return d + (p & tbmask(d));
}

///////////////////////////////

typedef struct {
	uint32_t p;
	uint32_t p0i;
	uint32_t R2;
	uint32_t g;
	uint32_t ig;
	uint32_t s;
} small_prime;
#define PRIMES   Zf(PRIMES)
extern const small_prime PRIMES[];
static inline uint32_t
mp_half(uint32_t a, uint32_t p)
{
	return (a + (p & -(a & 1))) >> 1;
}

/*
 * Compute R/2 = 2^31 mod p.
 */
static inline uint32_t
mp_hR(uint32_t p)
{
	/*
	 * Since p < 2^31 < (3/2)*p, we just subtract p from 2^31.
	 */
	return ((uint32_t)1 << 31) - p;
}
static inline uint32_t
mp_Rx31(unsigned e, uint32_t p, uint32_t p0i, uint32_t R2)
{
	/* x <- 2^63 mod p = Montgomery representation of 2^31 */
	uint32_t x = mp_half(R2, p);
	uint32_t d = 1;
	for (;;) {
		if ((e & 1) != 0) {
			d = mp_montymul(d, x, p, p0i);
		}
		e >>= 1;
		if (e == 0) {
			return d;
		}
		x = mp_montymul(x, x, p, p0i);
	}
}

static inline uint32_t
zint_mod_small_signed(const uint32_t* d, size_t len, size_t stride,
	uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx)
{
	if (len == 0) {
		return 0;
	}
	uint32_t z = zint_mod_small_unsigned(d, len, stride, p, p0i, R2);
	z = mp_sub(z, Rx & -(d[(len - 1) * stride] >> 30), p);
	return z;
}
static inline int32_t
mp_norm(uint32_t x, uint32_t p)
{
	uint32_t w = x - (p & tbmask((p >> 1) - x));
	return *(int32_t*)&w;
}

int
solve_NTRU(const ntru_profile* restrict prof, unsigned logn,
	const int8_t* restrict f, const int8_t* restrict g, uint32_t* tmp);




#endif
