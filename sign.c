//#include "scale_inner.h"
#include <stdio.h>

#include "sign.h"

#ifdef _DEBUG
extern int keys_number;
#endif
//* @author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline size_t
skoff_b00(unsigned logn)
{
	(void)logn;
	return 0;
}
//  * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline size_t
skoff_b01(unsigned logn)
{
	return MKN(logn);
}
//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline size_t
skoff_b10(unsigned logn)
{
	return 2 * MKN(logn);
}
//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline size_t
skoff_b11(unsigned logn)
{
	return 3 * MKN(logn);
}

//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline size_t
skoff_tree(unsigned logn)
{
	return 4 * MKN(logn);
}

//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
static inline unsigned
ffLDL_treesize(unsigned logn)
{
	/*
	 * For logn = 0 (polynomials are constant), the "tree" is a
	 * single element. Otherwise, the tree node has size 2^logn, and
	 * has two child trees for size logn-1 each. Thus, treesize s()
	 * must fulfill these two relations:
	 *
	 *   s(0) = 1
	 *   s(logn) = (2^logn) + 2*s(logn-1)
	 */
	return (logn + 1) << logn;
}

//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
//	to work with fixed point data use scale = 26
// poly_split_fft -> poly_split1_fft
// 
static 
void
ffLDL_fft_inner(fxr* restrict tree,
	fxr* restrict g0, fxr* restrict g1, unsigned logn, fxr* restrict tmp)
{
	size_t n, hn;

	n = MKN(logn);
	//int success = 0;


	if (n == 1) {
		tree[0] = g0[0];
		return;
	}
	hn = n >> 1;

	/*
	 * The LDL decomposition yields L (which is written in the tree)
	 * and the diagonal of D. Since d00 = g0, we just write d11
	 * into tmp.
	 */


	poly_LDLmv_fft(tmp, tree, g0, g1, g0, logn);

	{


		/*
		 * Split d00 (currently in g0) and d11 (currently in tmp). We
		 * reuse g0 and g1 as temporary storage spaces:
		 *   d00 splits into g1, g1+hn
		 *   d11 splits into g0, g0+hn
		 */
		
		poly_split1_fft(g1, g1 + hn, g0, logn);
		
		poly_split1_fft(g0, g0 + hn, tmp, logn);
		


		/*
		 * Each split result is the first row of a new auto-adjoint
		 * quasicyclic matrix for the next recursive step.
		 */

		
		ffLDL_fft_inner(tree + n,
			g1, g1 + hn, logn - 1, tmp);

		
		{
		
			ffLDL_fft_inner(tree + n + ffLDL_treesize(logn - 1),
				g0, g0 + hn, logn - 1, tmp);
			
		}
	}
	
}


//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
//	to work with fixed point data use scale = 26
// poly_split_fft -> poly_split1_fft

static 
void
ffLDL_fft(fxr* restrict tree, const fxr* restrict g00,
	const fxr* restrict g01, const fxr* restrict g11,
	unsigned logn, fxr* restrict tmp)
{
	size_t n, hn;
	
	fxr* d00, * d11;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = g00[0];
	}
	hn = n >> 1;
	d00 = tmp;
	d11 = tmp + n;
	tmp += n << 1;

	memcpy(d00, g00, n * sizeof * g00);

	poly_LDLmv_fft(d11, tree, g00, g01, g11, logn);

	{
		poly_split1_fft(tmp, tmp + hn, d00, logn);

		poly_split1_fft(d00, d00 + hn, d11, logn);

		memcpy(d11, tmp, n * sizeof * tmp);

			ffLDL_fft_inner(tree + n,
			d11, d11 + hn, logn - 1, tmp);

		{
				ffLDL_fft_inner(tree + n + ffLDL_treesize(logn - 1),
				d00, d00 + hn, logn - 1, tmp);
		}
	}
}


//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
//	to work with fixed point data use scale = 26

/*
 * Normalize an ffLDL tree: each leaf of value x is replaced with
 * sigma / sqrt(x).
 */

static 
void
ffLDL_binary_normalize(fxr* tree, unsigned orig_logn, unsigned logn)
{
	size_t n;

	n = MKN(logn);
	if (n == 1) {
		/*
		 * We actually store in the tree leaf the inverse of
		 * the value mandated by the specification: this
		 * saves a division both here and in the sampler.
		 */

		
		tree[0] = fxr_mul(fxr_sqrt(tree[0]), fxr_inv_sigma[orig_logn]);
		
		
	}
	else {
		ffLDL_binary_normalize(tree + n, orig_logn, logn - 1);
		
		ffLDL_binary_normalize(tree + n + ffLDL_treesize(logn - 1),
				orig_logn, logn - 1);
		
	}
	

}


extern int keys_number;

//*@author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
//	to work with fixed point data use scale = 26

void
expand_privkey(fxr* restrict expanded_key,
	const int8_t* f, const int8_t* g,
	const int8_t* F, const int8_t* G,
	unsigned logn, uint8_t* restrict tmp)

{
	size_t n;
	fxr* rf, * rg, * rF, * rG;
	fxr* b00, * b01, * b10, * b11;
	fxr* g00, * g01, * g11, * gxx;
	fxr* tree;

	n = MKN(logn);
	b00 = expanded_key + skoff_b00(logn);
	b01 = expanded_key + skoff_b01(logn);
	b10 = expanded_key + skoff_b10(logn);
	b11 = expanded_key + skoff_b11(logn);
	tree = expanded_key + skoff_tree(logn);

	/*
	 * We load the private key elements directly into the B0 matrix,
	 * since B0 = [[g, -f], [G, -F]].
	 */
	rf = b01;
	rg = b00;
	rF = b11;
	rG = b10;

	vect_set(logn, rf, f);

	vect_set(logn, rg, g);

	vect_set(logn, rF, F);

	vect_set(logn, rG, G);

	/*
	 * Compute the FFT for the key elements, and negate f and F.
	 */

	vect_FFT(logn, rf);

	vect_FFT(logn, rg);

	vect_FFT(logn, rF);

	vect_FFT(logn, rG);

	vect_neg(logn, rf);

	vect_neg(logn, rF);

	/*
	 * The Gram matrix is G = B·B*. Formulas are:
	 *   g00 = b00*adj(b00) + b01*adj(b01)
	 *   g01 = b00*adj(b10) + b01*adj(b11)
	 *   g10 = b10*adj(b00) + b11*adj(b01)
	 *   g11 = b10*adj(b10) + b11*adj(b11)
	 *
	 * For historical reasons, this implementation uses
	 * g00, g01 and g11 (upper triangle).
	 */
	g00 = (fxr*)tmp;
	g01 = g00 + n;
	g11 = g01 + n;
	gxx = g11 + n;

	memcpy(g00, b00, n * sizeof * b00);

	poly_mulselfadj_fft(g00, logn);


	memcpy(gxx, b01, n * sizeof * b01);


	poly_mulselfadj_fft(gxx, logn);


	poly_add(g00, gxx, logn);

	memcpy(g01, b00, n * sizeof * b00);
	poly_muladj_fft(g01, b10, logn);

	memcpy(gxx, b01, n * sizeof * b01);
	poly_muladj_fft(gxx, b11, logn);

	poly_add(g01, gxx, logn);


	memcpy(g11, b10, n * sizeof * b10);
	poly_mulselfadj_fft(g11, logn);

	memcpy(gxx, b11, n * sizeof * b11);
	poly_mulselfadj_fft(gxx, logn);

	poly_add(g11, gxx, logn);

	/*
	 * Compute the Falcon tree.
	 */

	ffLDL_fft(tree, g00, g01, g11, logn, gxx);
	
	/*
	 * Normalize tree.
	 */

	
	
		ffLDL_binary_normalize(tree, logn, logn);
	
}

//@author   Thomas Pornin <thomas.pornin@nccgroup.com>
void
hash_to_point_vartime(
	inner_shake256_context* sc,
	uint16_t* x, unsigned logn)
{
	/*
	 * This is the straightforward per-the-spec implementation. It
	 * is not constant-time, thus it might reveal information on the
	 * plaintext (at least, enough to check the plaintext against a
	 * list of potential plaintexts) in a scenario where the
	 * attacker does not have access to the signature value or to
	 * the public key, but knows the nonce (without knowledge of the
	 * nonce, the hashed output cannot be matched against potential
	 * plaintexts).
	 */
	size_t n;

	n = (size_t)1 << logn;
	while (n > 0) {
		uint8_t buf[2];
		uint32_t w;
		//sign_count += 2;
		//sign_count1 += 2;
		inner_shake256_extract(sc, (void*)buf, sizeof buf);
		w = ((unsigned)buf[0] << 8) | (unsigned)buf[1];
		if (w < 61445) {
			while (w >= 12289) {
				w -= 12289;
			}
			*x++ = (uint16_t)w;
			n--;
		}
	}
}


//@author   Thomas Pornin <thomas.pornin@nccgroup.com>
typedef struct {
	prng p;
	fxr sigma_min;
} sampler_context;

typedef int (*samplerZ)(void* ctx, fxr mu, fxr sigma);

int
gaussian0_sampler(prng* p)
{


	static const uint32_t dist[] = {
		10745844u,  3068844u,  3741698u,
		 5559083u,  1580863u,  8248194u,
		 2260429u, 13669192u,  2736639u,
		  708981u,  4421575u, 10046180u,
		  169348u,  7122675u,  4136815u,
		   30538u, 13063405u,  7650655u,
			4132u, 14505003u,  7826148u,
			 417u, 16768101u, 11363290u,
			  31u,  8444042u,  8086568u,
			   1u, 12844466u,   265321u,
			   0u,  1232676u, 13644283u,
			   0u,    38047u,  9111839u,
			   0u,      870u,  6138264u,
			   0u,       14u, 12545723u,
			   0u,        0u,  3104126u,
			   0u,        0u,    28824u,
			   0u,        0u,      198u,
			   0u,        0u,        1u
	};

	uint32_t v0, v1, v2, hi;
	uint64_t lo;
	size_t u;
	int z;

	/*
	 * Get a random 72-bit value, into three 24-bit limbs v0..v2.
	 */
	lo = prng_get_u64(p);
	hi = prng_get_u8(p);
	//sign_count2 += 9;
	//sign_count += 9;
	v0 = (uint32_t)lo & 0xFFFFFF;
	v1 = (uint32_t)(lo >> 24) & 0xFFFFFF;
	v2 = (uint32_t)(lo >> 48) | (hi << 16);

	/*
	 * Sampled value is z, such that v0..v2 is lower than the first
	 * z elements of the table.
	 */
	z = 0;
	for (u = 0; u < (sizeof dist) / sizeof(dist[0]); u += 3) {
		uint32_t w0, w1, w2, cc;

		w0 = dist[u + 2];
		w1 = dist[u + 1];
		w2 = dist[u + 0];
		cc = (v0 - w0) >> 31;
		cc = (v1 - w1 - cc) >> 31;
		cc = (v2 - w2 - cc) >> 31;
		z += (int)cc;
	}
	return z;


}


//// Calculation exp(x) for 0 <=x < 1
//// changes:
//// scale = 2 * SCALE
//// Constants for scale = 2 * SCALE and function 
//static fxr2 fxr2_exp_c[] = {
//	{9339440},
//	{113938847},
//	{1241225186},
//	{12410057660},
//	{111696327149},
//	{893571538674},
//	{6254999505761},
//	{37529996869837},
//	{187649984471265},
//	{750599937896511},
//	{2251799813685334},
//	{4503599627370473},
//	{4503599627370496}
//};
//
//fxr2 fxr2_exp_(fxr2 x)
//{
//	//x.v = 0.055750638246536255 * ((uint64_t)1 << 2 * SCALE);
//	fxr2 y = fxr2_exp_c[0], d = x;
//	//y = 0.000000025299506379442070029551 - y * d;
//	y = fxr2_sub(fxr2_exp_c[1], fxr2_mul(y, d));
//	//y = 0.000000275607356160477811864927 - y * d;
//	y = fxr2_sub(fxr2_exp_c[2], fxr2_mul(y, d));
//	//y = 0.000002755586350219122514855659 - y * d;
//	y = fxr2_sub(fxr2_exp_c[3], fxr2_mul(y, d));
//	//y = 0.000024801566833585381209939524 - y * d;
//	y = fxr2_sub(fxr2_exp_c[4], fxr2_mul(y, d));
//	//y = 0.000198412739277311890541063977 - y * d;
//	y = fxr2_sub(fxr2_exp_c[5], fxr2_mul(y, d));
//	//y = 0.001388888894063186997887560103 - y * d;
//	y = fxr2_sub(fxr2_exp_c[6], fxr2_mul(y, d));
//
//	//y = 0.008333333327800835146903501993 - y * d;
//	y = fxr2_sub(fxr2_exp_c[7], fxr2_mul(y, d));
//	//y = 0.041666666666110491190622155955 - y * d;
//	y = fxr2_sub(fxr2_exp_c[8], fxr2_mul(y, d));
//	//y = 0.166666666666984014666397229121 - y * d;
//	y = fxr2_sub(fxr2_exp_c[9], fxr2_mul(y, d));
//	//y = 0.500000000000019206858326015208 - y * d;
//	y = fxr2_sub(fxr2_exp_c[10], fxr2_mul(y, d));
//	//y = 0.999999999999994892974086724280 - y * d;
//	y = fxr2_sub(fxr2_exp_c[11], fxr2_mul(y, d));
//	//y = 1.000000000000000000000000000000 - y * d;
//	y = fxr2_sub(fxr2_exp_c[12], fxr2_mul(y, d));
//	return y;
//}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// Function name fpr_expm_p63 changed to name fxr_expm_p63_
// scale = SCALE

uint64_t fxr_expm_p63_(fxr x, fxr ccs)
{
	//fxr r;
	//double d, y;
	/*
	ccs.v	0.73344958102467073	double
		x.v	0.14020644588325798	double

	*/
	fxr2 fxr2_x, fxr2_ccs;
	/*x.v = 0.14020644588325798 * (1 << SCALE);
	ccs.v = 0.73344958102467073 * (1 << SCALE);*/
	fxr2_x.v = x.v << (SCALE);
	fxr2_ccs.v = ccs.v << (SCALE);
	fxr2_x = fxr2_exp_(fxr2_x);
	fxr2_x = fxr2_mul(fxr2_x, fxr2_ccs);
	return fxr2_x.v << (63 - 2 * SCALE);

}


// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// To work with fixed point data, use a scale = SCALE.

static int
BerExp(prng* p, fxr x, fxr ccs)
{
	int s, i;
	fxr r;
	uint32_t sw, w;
	uint64_t z;

	/*
	 * Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	 * and 0 <= r < log(2). Since x >= 0, we can use fxr_trunc().
	 */
	s = (int)fxr_trunc(fxr_mul(x, fxr_inv_log2));
	r = fxr_sub(x, fxr_mul(fxr_of(s), fxr_log2));

	/*
	 * It may happen (quite rarely) that s >= 64; if sigma = 1.2
	 * (the minimum value for sigma), r = 0 and b = 1, then we get
	 * s >= 64 if the half-Gaussian produced a z >= 13, which happens
	 * with probability about 0.000000000230383991, which is
	 * approximatively equal to 2^(-32). In any case, if s >= 64,
	 * then BerExp will be non-zero with probability less than
	 * 2^(-64), so we can simply saturate s at 63.
	 */
	sw = (uint32_t)s;
	sw ^= (sw ^ 63) & -((63 - sw) >> 31);
	s = (int)sw;

	/*
	 * Compute exp(-r); we know that 0 <= r < log(2) at this point, so
	 * we can use fxr_expm_p63(), which yields a result scaled to 2^63.
	 * We scale it up to 2^64, then right-shift it by s bits because
	 * we really want exp(-x) = 2^(-s)*exp(-r).
	 *
	 * The "-1" operation makes sure that the value fits on 64 bits
	 * (i.e. if r = 0, we may get 2^64, and we prefer 2^64-1 in that
	 * case). The bias is negligible since fxr_expm_p63() only computes
	 * with 51 bits of precision or so.
	 */
	if (s < 0)
		z = ((fxr_expm_p63_(r, ccs) << 1) - 1) << (-s);
	else
		z = ((fxr_expm_p63_(r, ccs) << 1) - 1) >> s;

	/*
	 * Sample a bit with probability exp(-x). Since x = s*log(2) + r,
	 * exp(-x) = 2^-s * exp(-r), we compare lazily exp(-x) with the
	 * PRNG output to limit its consumption, the sign of the difference
	 * yields the expected result.
	 */
	i = 64;
	do {
		i -= 8;
		w = prng_get_u8(p) - ((uint32_t)(z >> i) & 0xFF);
		//sign_count2 += 1;
		//sign_count += 1;
	} while (!w && i > 0);
	return (int)(w >> 31);
}


// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// To work with fixed point data, use a scale = SCALE.
int
sampler(void* ctx, fxr mu, fxr isigma)
{
	sampler_context* spc;
	int s;
	fxr r, dss, ccs;

	spc = ctx;

	/*
	 * Center is mu. We compute mu = s + r where s is an integer
	 * and 0 <= r < 1.
	 */
	s = (int)fxr_floor(mu);
	r = fxr_sub(mu, fxr_of(s));

	/*
	 * dss = 1/(2*sigma^2) = 0.5*(isigma^2).
	 */
	dss = fxr_half(fxr_sqr(isigma));

	/*
	 * ccs = sigma_min / sigma = sigma_min * isigma.
	 */
	ccs = fxr_mul(isigma, spc->sigma_min);

	/*
	 * We now need to sample on center r.
	 */
	for (;;) {
		int z0, z, b;
		fxr x;

		/*
		 * Sample z for a Gaussian distribution. Then get a
		 * random bit b to turn the sampling into a bimodal
		 * distribution: if b = 1, we use z+1, otherwise we
		 * use -z. We thus have two situations:
		 *
		 *  - b = 1: z >= 1 and sampled against a Gaussian
		 *    centered on 1.
		 *  - b = 0: z <= 0 and sampled against a Gaussian
		 *    centered on 0.
		 */
		z0 = gaussian0_sampler(&spc->p);
		b = (int)prng_get_u8(&spc->p) & 1;
		//sign_count2 += 1;
		//sign_count += 1;
		z = b + ((b << 1) - 1) * z0;

		/*
		 * Rejection sampling. We want a Gaussian centered on r;
		 * but we sampled against a Gaussian centered on b (0 or
		 * 1). But we know that z is always in the range where
		 * our sampling distribution is greater than the Gaussian
		 * distribution, so rejection works.
		 *
		 * We got z with distribution:
		 *    G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		 * We target distribution:
		 *    S(z) = exp(-((z-r)^2)/(2*sigma^2))
		 * Rejection sampling works by keeping the value z with
		 * probability S(z)/G(z), and starting again otherwise.
		 * This requires S(z) <= G(z), which is the case here.
		 * Thus, we simply need to keep our z with probability:
		 *    P = exp(-x)
		 * where:
		 *    x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		 *
		 * Here, we scale up the Bernouilli distribution, which
		 * makes rejection more probable, but makes rejection
		 * rate sufficiently decorrelated from the Gaussian
		 * center and standard deviation that the whole sampler
		 * can be said to be constant-time.
		 */
		x = fxr_mul(fxr_sqr(fxr_sub(fxr_of(z), r)), dss);
		//x = fxr_sub(x, fxr_mul(fxr_of(z0 * z0), ng_inv_2sqrsigma0));
		fxr t;
		t = fxr_of(z0 * z0);
		t = fxr_mul(t, fxr_inv_2sqrsigma0);
		x = fxr_sub(x, t);
		if (BerExp(&spc->p, x, ccs)) {
			/*
			 * Rejection sampling was centered on r, but the
			 * actual center is mu = s + r.
			 */
			return s + z;
		}
	}
}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// To work with fixed point data, use a scale = SCALE.
// Instead of the poly_split_fft function, use the poly_split1_fft function
// Instead of the poly_merge_fft function, use the poly_merge1_fft function

static void
ffSampling_fft(samplerZ samp, void* samp_ctx,
	fxr* restrict z0, fxr* restrict z1,
	const fxr* restrict tree,
	const fxr* restrict t0, const fxr* restrict t1, unsigned logn,
	fxr* restrict tmp)
{
	size_t n, hn;
	const fxr* tree0, * tree1;
	/*
	 * When logn == 2, we inline the last two recursion levels.
	 */
	if (logn == 2) {
		static int count = 0;
		++count;

		fxr x0, x1, y0, y1, w0, w1, w2, w3, sigma;
		fxr a_re, a_im, b_re, b_im, c_re, c_im;
		

		tree0 = tree + 4;
		tree1 = tree + 8;

		/*
		 * We split t1 into w*, then do the recursive invocation,
		 * with output in w*. We finally merge back into z1.
		 */
		a_re = t1[0];
		a_im = t1[2];
		b_re = t1[1];
		b_im = t1[3];
		c_re = fxr_add(a_re, b_re);
		c_im = fxr_add(a_im, b_im);
		w0 = fxr_half(c_re);
		w1 = fxr_half(c_im);
		c_re = fxr_sub(a_re, b_re);
		c_im = fxr_sub(a_im, b_im);
		w2 = fxr_mul(fxr_add(c_re, c_im), fxr_invsqrt8);
		w3 = fxr_mul(fxr_sub(c_im, c_re), fxr_invsqrt8);

		x0 = w2;
		x1 = w3;
		sigma = tree1[3];
		w2 = fxr_of(samp(samp_ctx, x0, sigma));
		w3 = fxr_of(samp(samp_ctx, x1, sigma));

		a_re = fxr_sub(x0, w2);
		a_im = fxr_sub(x1, w3);
		b_re = tree1[0];
		b_im = tree1[1];
		c_re = fxr_sub(fxr_mul(a_re, b_re), fxr_mul(a_im, b_im));
		c_im = fxr_add(fxr_mul(a_re, b_im), fxr_mul(a_im, b_re));
		x0 = fxr_add(c_re, w0);
		x1 = fxr_add(c_im, w1);
		sigma = tree1[2];

		w0 = fxr_of(samp(samp_ctx, x0, sigma));

		w1 = fxr_of(samp(samp_ctx, x1, sigma));

		a_re = w0;
		a_im = w1;
		b_re = w2;
		b_im = w3;
		c_re = fxr_mul(fxr_sub(b_re, b_im), fxr_invsqrt2);
		c_im = fxr_mul(fxr_add(b_re, b_im), fxr_invsqrt2);
		z1[0] = w0 = fxr_add(a_re, c_re);
		z1[2] = w2 = fxr_add(a_im, c_im);
		z1[1] = w1 = fxr_sub(a_re, c_re);
		z1[3] = w3 = fxr_sub(a_im, c_im);

		/*
		 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in w*.
		 */
		w0 = fxr_sub(t1[0], w0);
		w1 = fxr_sub(t1[1], w1);
		w2 = fxr_sub(t1[2], w2);
		w3 = fxr_sub(t1[3], w3);

		a_re = w0;
		a_im = w2;
		b_re = tree[0];
		b_im = tree[2];
		w0 = fxr_sub(fxr_mul(a_re, b_re), fxr_mul(a_im, b_im));
		w2 = fxr_add(fxr_mul(a_re, b_im), fxr_mul(a_im, b_re));
		a_re = w1;
		a_im = w3;
		b_re = tree[1];
		b_im = tree[3];
		w1 = fxr_sub(fxr_mul(a_re, b_re), fxr_mul(a_im, b_im));
		w3 = fxr_add(fxr_mul(a_re, b_im), fxr_mul(a_im, b_re));

		w0 = fxr_add(w0, t0[0]);
		w1 = fxr_add(w1, t0[1]);
		w2 = fxr_add(w2, t0[2]);
		w3 = fxr_add(w3, t0[3]);

		/*
		 * Second recursive invocation.
		 */
		a_re = w0;
		a_im = w2;
		b_re = w1;
		b_im = w3;
		c_re = fxr_add(a_re, b_re);
		c_im = fxr_add(a_im, b_im);
		w0 = fxr_half(c_re);
		w1 = fxr_half(c_im);
		c_re = fxr_sub(a_re, b_re);
		c_im = fxr_sub(a_im, b_im);
		w2 = fxr_mul(fxr_add(c_re, c_im), fxr_invsqrt8);
		w3 = fxr_mul(fxr_sub(c_im, c_re), fxr_invsqrt8);

		x0 = w2;
		x1 = w3;
		sigma = tree0[3];

		w2 = y0 = fxr_of(samp(samp_ctx, x0, sigma));
		w3 = y1 = fxr_of(samp(samp_ctx, x1, sigma));

		a_re = fxr_sub(x0, y0);
		a_im = fxr_sub(x1, y1);
		b_re = tree0[0];
		b_im = tree0[1];
		c_re = fxr_sub(fxr_mul(a_re, b_re), fxr_mul(a_im, b_im));
		c_im = fxr_add(fxr_mul(a_re, b_im), fxr_mul(a_im, b_re));
		x0 = fxr_add(c_re, w0);
		x1 = fxr_add(c_im, w1);
		sigma = tree0[2];
		w0 = fxr_of(samp(samp_ctx, x0, sigma));
		w1 = fxr_of(samp(samp_ctx, x1, sigma));

		a_re = w0;
		a_im = w1;
		b_re = w2;
		b_im = w3;
		c_re = fxr_mul(fxr_sub(b_re, b_im), fxr_invsqrt2);
		c_im = fxr_mul(fxr_add(b_re, b_im), fxr_invsqrt2);
		z0[0] = fxr_add(a_re, c_re);
		z0[2] = fxr_add(a_im, c_im);
		z0[1] = fxr_sub(a_re, c_re);
		z0[3] = fxr_sub(a_im, c_im);

		return;

	}

	/*
	 * Case logn == 1 is reachable only when using Falcon-2 (the
	 * smallest size for which Falcon is mathematically defined, but
	 * of course way too insecure to be of any use).
	 */
	if (logn == 1) {
		fxr x0, x1, y0, y1, sigma;
		fxr a_re, a_im, b_re, b_im, c_re, c_im;

		x0 = t1[0];
		x1 = t1[1];
		sigma = tree[3];
		z1[0] = y0 = fxr_of(samp(samp_ctx, x0, sigma));
		z1[1] = y1 = fxr_of(samp(samp_ctx, x1, sigma));
		a_re = fxr_sub(x0, y0);
		a_im = fxr_sub(x1, y1);
		b_re = tree[0];
		b_im = tree[1];
		c_re = fxr_sub(fxr_mul(a_re, b_re), fxr_mul(a_im, b_im));
		c_im = fxr_add(fxr_mul(a_re, b_im), fxr_mul(a_im, b_re));
		x0 = fxr_add(c_re, t0[0]);
		x1 = fxr_add(c_im, t0[1]);
		sigma = tree[2];
		z0[0] = fxr_of(samp(samp_ctx, x0, sigma));
		z0[1] = fxr_of(samp(samp_ctx, x1, sigma));

		return;
	}

	/*
	 * Normal end of recursion is for logn == 0. Since the last
	 * steps of the recursions were inlined in the blocks above
	 * (when logn == 1 or 2), this case is not reachable, and is
	 * retained here only for documentation purposes.

	if (logn == 0) {
		fxr x0, x1, sigma;

		x0 = t0[0];
		x1 = t1[0];
		sigma = tree[0];
		z0[0] = fxr_of(samp(samp_ctx, x0, sigma));
		z1[0] = fxr_of(samp(samp_ctx, x1, sigma));
		return;
	}

	 */

	 /*
	  * General recursive case (logn >= 3).
	  */

	n = (size_t)1 << logn;
	hn = n >> 1;
	tree0 = tree + n;
	tree1 = tree + n + ffLDL_treesize(logn - 1);

	/*
	 * We split t1 into z1 (reused as temporary storage), then do
	 * the recursive invocation, with output in tmp. We finally
	 * merge back into z1.
	 */
	 //poly_split_fft(z1, z1 + hn, t1, logn);
	poly_split1_fft(z1, z1 + hn, t1, logn);
	ffSampling_fft(samp, samp_ctx, tmp, tmp + hn,
		tree1, z1, z1 + hn, logn - 1, tmp + n);
	//poly_merge_fft(z1, tmp, tmp + hn, logn);
	poly_merge1_fft(z1, tmp, tmp + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in tmp[].
	 */
	memcpy(tmp, t1, n * sizeof * t1);
	//poly_sub(tmp, z1, logn);
	vect_sub(logn, tmp, z1);
	vect_mul_fft(logn, tmp, tree);
	vect_add(logn, tmp, t0);

	/*
	 * Second recursive invocation.
	 */
	poly_split1_fft(z0, z0 + hn, tmp, logn);
	ffSampling_fft(samp, samp_ctx, tmp, tmp + hn,
		tree0, z0, z0 + hn, logn - 1, tmp + n);
	poly_merge1_fft(z0, tmp, tmp + hn, logn);
}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
int
is_short_half(
	uint32_t sqn, const int16_t* s2, unsigned logn)
{
	size_t n, u;
	uint32_t ng;

	n = (size_t)1 << logn;
	ng = -(sqn >> 31);
	for (u = 0; u < n; u++) {
		int32_t z;

		z = s2[u];
		sqn += (uint32_t)(z * z);
		ng |= sqn;
	}
	sqn |= -(ng >> 31);

	return sqn <= l2bound[logn];
}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// To work with fixed point data, use a scale = SCALE.
static int
do_sign_tree(samplerZ samp, void* samp_ctx, int16_t* s2,
	const fxr* restrict expanded_key,
	const uint16_t* hm,
	unsigned logn, fxr* restrict tmp)
{
	size_t n, u;
	fxr* t0, * t1, * tx, * ty;
	const fxr* b00, * b01, * b10, * b11, * tree;
	fxr ni;
	uint32_t sqn, ng;
	int16_t* s1tmp, * s2tmp;
	
	n = MKN(logn);
	t0 = tmp;
	t1 = t0 + n;
	b00 = expanded_key + skoff_b00(logn);
	b01 = expanded_key + skoff_b01(logn);
	b10 = expanded_key + skoff_b10(logn);
	b11 = expanded_key + skoff_b11(logn);
	tree = expanded_key + skoff_tree(logn);

	/*
	 * Set the target vector to [hm, 0] (hm is the hashed message).
	 */
	for (u = 0; u < n; u++) {
		t0[u] = fxr_of(hm[u]);
		/* This is implicit.
		t1[u] = fxr_zero;
		*/
	}

	/*
	 * Apply the lattice basis to obtain the real target
	 * vector (after normalization with regards to modulus).
	 */

	 //FFT(t0, logn);
	vect_FFT(logn, t0);
	ni = fxr_inverse_of_q;
	memcpy(t1, t0, n * sizeof * t0);
	//poly_mul_fft(t1, b01, logn);
	vect_mul_fft(logn, t1, b01);
	//vect_mulconst(t1, fxr_neg(ni), logn);
	vect_mulconst(logn, t1, fxr_neg(ni));
	//poly_mul_fft(t0, b11, logn);
	vect_mul_fft(logn, t0, b11);
	//vect_mulconst(t0, ni, logn);
	vect_mulconst(logn, t0, ni);

	tx = t1 + n;
	ty = tx + n;
	
	/*
	 * Apply sampling. Output is written back in [tx, ty].
	 */
	ffSampling_fft(samp, samp_ctx, tx, ty, tree, t0, t1, logn, ty + n);
	
	/*
	 * Get the lattice point corresponding to that tiny vector.
	 */
	memcpy(t0, tx, n * sizeof * tx);
	memcpy(t1, ty, n * sizeof * ty);
	// vect_mul_fft(unsigned logn, fxr* restrict a, const fxr* restrict b);
	vect_mul_fft(logn, tx, b00);
	//poly_mul_fft(ty, b10, logn);
	vect_mul_fft(logn, ty, b10);
	poly_add(tx, ty, logn);
	memcpy(ty, t0, n * sizeof * t0);
	//poly_mul_fft(ty, b01, logn);
	vect_mul_fft(logn, ty, b01);

	memcpy(t0, tx, n * sizeof * tx);
	//poly_mul_fft(t1, b11, logn);
	vect_mul_fft(logn, t1, b11);
	poly_add(t1, ty, logn);

	//iFFT(t0, logn);
	vect_iFFT(logn, t0);
	//iFFT(t1, logn);
	//FILE* file;
	//file = fopen("new_t1.bin", "wb");
	//fwrite(t1, 1024, 8, file);
	vect_iFFT(logn, t1);
	/*fwrite(t1, 1024, 8, file);
	fclose(file);*/

	/*
	 * Compute the signature.
	 */
	s1tmp = (int16_t*)tx;
	sqn = 0;
	ng = 0;
	for (u = 0; u < n; u++) {
		int32_t z;

		z = (int32_t)hm[u] - (int32_t)fxr_rint(t0[u]);
		sqn += (uint32_t)(z * z);
		ng |= sqn;
		s1tmp[u] = (int16_t)z;
	}
	/*{
		FILE* file8067 = fopen("S1_8067_2.bin", "wb");
		fwrite(s1tmp, sizeof(s1tmp[0]), n, file8067);
		fclose(file8067);
	}*/
	sqn |= -(ng >> 31);

	/*
	 * With "normal" degrees (e.g. 512 or 1024), it is very
	 * improbable that the computed vector is not short enough;
	 * however, it may happen in practice for the very reduced
	 * versions (e.g. degree 16 or below). In that case, the caller
	 * will loop, and we must not write anything into s2[] because
	 * s2[] may overlap with the hashed message hm[] and we need
	 * hm[] for the next iteration.
	 */
	s2tmp = (int16_t*)tmp;
	for (u = 0; u < n; u++) {
		s2tmp[u] = (int16_t)-fxr_rint(t1[u]);
	}
	/*{
		FILE* file8067 = fopen("S2_8067_2.bin", "wb");
		fwrite(s2tmp, sizeof(s2tmp[0]), n, file8067);
		fclose(file8067);
	}*/
	if (is_short_half(sqn, s2tmp, logn)) {
		memcpy(s2, s2tmp, n * sizeof * s2);
		memcpy(tmp, s1tmp, n * sizeof * s1tmp);
		return 1;
	}
	return 0;
}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// changes:
// To work with fixed point data, use a scale = SCALE.
void
sign_tree(int16_t* sig, inner_shake256_context* rng,
	const fxr* restrict expanded_key,
	const uint16_t* hm, unsigned logn, uint8_t* tmp)
{
	fxr* ftmp;

	ftmp = (fxr*)tmp;
	for (;;) {
		/*
		 * Signature produces short vectors s1 and s2. The
		 * signature is acceptable only if the aggregate vector
		 * s1,s2 is short; we must use the same bound as the
		 * verifier.
		 *
		 * If the signature is acceptable, then we return only s2
		 * (the verifier recomputes s1 from s2, the hashed message,
		 * and the public key).
		 */
		sampler_context spc;
		samplerZ samp;
		void* samp_ctx;

		/*
		 * Normal sampling. We use a fast PRNG seeded from our
		 * SHAKE context ('rng').
		 */
		spc.sigma_min = fxr_sigma_min[logn];
		prng_init(&spc.p, rng);
		samp = sampler;
		samp_ctx = &spc;

		/*
		 * Do the actual signature.
		 */
		if (do_sign_tree(samp, samp_ctx, sig,
			expanded_key, hm, logn, ftmp))
		{
			break;
		}
	}
}

