#include "poly.h"

// Thomas Pornin. https://github.com/pornin/ntrugen
void
poly_mp_set_small(unsigned logn, uint32_t* restrict d,
	const int8_t* restrict f, uint32_t p)
{
	size_t n = (size_t)1 << logn;

	for (size_t u = 0; u < n; u++) {
		d[u] = mp_set(f[u], p);
	}
}

// Thomas Pornin. https://github.com/pornin/ntrugen
uint32_t
poly_max_bitlength(unsigned logn, const uint32_t* f, size_t flen)
{
	if (flen == 0) {
		return 0;
	}

	size_t n = (size_t)1 << logn;
	uint32_t t = 0;
	uint32_t tk = 0;
	for (size_t u = 0; u < n; u++, f++) {
		/* Extend sign bit into a 31-bit mask. */
		uint32_t m = -(f[(flen - 1) << logn] >> 30) & 0x7FFFFFFF;

		/* Get top non-zero sign-adjusted word, with index. */
		uint32_t c = 0;
		uint32_t ck = 0;
		for (size_t v = 0; v < flen; v++) {
			uint32_t w = f[v << logn] ^ m; /* sign-adjusted word */
			uint32_t nz = ((w - 1) >> 31) - 1;
			c ^= nz & (c ^ w);
			ck ^= nz & (ck ^ (uint32_t)v);
		}


		uint32_t rr = tbmask((tk - ck) | (((tk ^ ck) - 1) & (t - c)));
		t ^= rr & (t ^ c);
		tk ^= rr & (tk ^ ck);
	}

	/*
	 * Get bit length of the top word (which has been sign-adjusted)
	 * and return the result.
	 */
	return 31 * tk + 32 - lzcnt(t);
}

// Thomas Pornin. https://github.com/pornin/ntrugen
void
poly_sub_kfg_scaled_depth1(unsigned logn_top,
	uint32_t* restrict F, uint32_t* restrict G, size_t FGlen,
	//uint64_t* restrict k, 
	int32_t* restrict k32,
	uint32_t sc,
	const int8_t* restrict f, const int8_t* restrict g,
	uint32_t* restrict tmp)


{



	unsigned logn = logn_top - 1;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	uint32_t* gm = tmp;
	uint32_t* t1 = gm + n;
	uint32_t* t2 = t1 + n;

	/*uint32_t k32[1024];
	for (int i = 0; i < n; ++i)
		k32[i] = (uint32_t)k[i];*/


	/*
	 * Step 1: convert F and G to RNS. Since FGlen is equal to 1 or 2,
	 * we do it with some specialized code. We assume that the RNS
	 * representation does not lose information (i.e. each signed
	 * coefficient is lower than (p0*p1)/2, with FGlen = 2 and the two
	 * prime moduli are p0 and p1).
	 */
	if (FGlen == 1) {
		uint32_t p = PRIMES[0].p;
		for (size_t u = 0; u < n; u++) {
			uint32_t xf = F[u];
			uint32_t xg = G[u];
			xf |= (xf & 0x40000000) << 1;
			xg |= (xg & 0x40000000) << 1;
			F[u] = mp_set(*(int32_t*)&xf, p);
			G[u] = mp_set(*(int32_t*)&xg, p);
		}
	}
	else {
		uint32_t p0 = PRIMES[0].p;
		uint32_t p0_0i = PRIMES[0].p0i;
		uint32_t z0 = mp_half(PRIMES[0].R2, p0);
		uint32_t p1 = PRIMES[1].p;
		uint32_t p1_0i = PRIMES[1].p0i;
		uint32_t z1 = mp_half(PRIMES[1].R2, p1);
		for (size_t u = 0; u < n; u++) {
			uint32_t xl, xh, yl0, yh0, r0, yl1, yh1, r1;

			xl = F[u];
			xh = F[u + n] | ((F[u + n] & 0x40000000) << 1);
			yl0 = xl - (p0 & ~tbmask(xl - p0));
			yh0 = mp_set(*(int32_t*)&xh, p0);
			r0 = mp_add(yl0, mp_montymul(yh0, z0, p0, p0_0i), p0);
			yl1 = xl - (p1 & ~tbmask(xl - p1));
			yh1 = mp_set(*(int32_t*)&xh, p1);
			r1 = mp_add(yl1, mp_montymul(yh1, z1, p1, p1_0i), p1);
			F[u] = r0;
			F[u + n] = r1;

			xl = G[u];
			xh = G[u + n] | ((G[u + n] & 0x40000000) << 1);
			yl0 = xl - (p0 & ~tbmask(xl - p0));
			yh0 = mp_set(*(int32_t*)&xh, p0);
			r0 = mp_add(yl0, mp_montymul(yh0, z0, p0, p0_0i), p0);
			yl1 = xl - (p1 & ~tbmask(xl - p1));
			yh1 = mp_set(*(int32_t*)&xh, p1);
			r1 = mp_add(yl1, mp_montymul(yh1, z1, p1, p1_0i), p1);
			G[u] = r0;
			G[u + n] = r1;
		}
	}

	/*
	 * Step 2: for FGlen small primes, convert F and G to RNS+NTT,
	 * and subtract (2^sc)*(ft,gt). The (ft,gt) polynomials are computed
	 * in RNS+NTT dynamically.
	 */
	for (size_t u = 0; u < FGlen; u++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		uint32_t R3 = mp_montymul(R2, R2, p, p0i);
		mp_mkgm(logn, gm, PRIMES[u].g, p, p0i);

		/*
		 * k <- (2^sc)*k (and into NTT).
		 */
		uint32_t scv = mp_montymul(
			(uint32_t)1 << (sc & 31), R2, p, p0i);
		for (uint32_t m = sc >> 5; m > 0; m--) {
			scv = mp_montymul(scv, R2, p, p0i);
		}
		for (size_t v = 0; v < n; v++) {
			uint32_t x = mp_set(*(int32_t*)&k32[v], p);
			k32[v] = mp_montymul(scv, x, p, p0i);
		}
		mp_NTT(logn, k32, gm, p, p0i);

		/*
		 * Convert F and G to NTT.
		 */
		uint32_t* Fu = F + (u << logn);
		uint32_t* Gu = G + (u << logn);
		mp_NTT(logn, Fu, gm, p, p0i);
		mp_NTT(logn, Gu, gm, p, p0i);

		/*
		 * Given the top-level f, we obtain ft = N(f) (the f at
		 * depth 1) with:
		 *    f = f_e(X^2) + X*f_o(X^2)
		 * with f_e and f_o being modulo X^n+1. Then:
		 *    N(f) = f_e^2 - X*f_o^2
		 * The NTT representation of X is obtained from the gm[] tab:
		 *    NTT(X)[2*j + 0] = gm[j + n/2]
		 *    NTT(X)[2*j + 1] = -NTT(X)[2*j + 0]
		 * Note that the values in gm[] are in Montgomery
		 * representation.
		 */
		for (size_t v = 0; v < n; v++) {
			t1[v] = mp_set(f[(v << 1) + 0], p);
			t2[v] = mp_set(f[(v << 1) + 1], p);
		}
		mp_NTT(logn, t1, gm, p, p0i);
		mp_NTT(logn, t2, gm, p, p0i);
		for (size_t v = 0; v < hn; v++) {
			uint32_t xe0 = t1[(v << 1) + 0];
			uint32_t xe1 = t1[(v << 1) + 1];
			uint32_t xo0 = t2[(v << 1) + 0];
			uint32_t xo1 = t2[(v << 1) + 1];
			uint32_t xv0 = gm[hn + v];
			uint32_t xv1 = p - xv0;
			xe0 = mp_montymul(xe0, xe0, p, p0i);
			xe1 = mp_montymul(xe1, xe1, p, p0i);
			xo0 = mp_montymul(xo0, xo0, p, p0i);
			xo1 = mp_montymul(xo1, xo1, p, p0i);
			uint32_t xf0 = mp_sub(xe0,
				mp_montymul(xo0, xv0, p, p0i), p);
			uint32_t xf1 = mp_sub(xe1,
				mp_montymul(xo1, xv1, p, p0i), p);

			uint32_t xkf0 = mp_montymul(
				mp_montymul(xf0, k32[(v << 1) + 0], p, p0i),
				R3, p, p0i);
			uint32_t xkf1 = mp_montymul(
				mp_montymul(xf1, k32[(v << 1) + 1], p, p0i),
				R3, p, p0i);
			Fu[(v << 1) + 0] = mp_sub(Fu[(v << 1) + 0], xkf0, p);
			Fu[(v << 1) + 1] = mp_sub(Fu[(v << 1) + 1], xkf1, p);
		}

		/*
		 * Same treatment for G and gt.
		 */
		for (size_t v = 0; v < n; v++) {
			t1[v] = mp_set(g[(v << 1) + 0], p);
			t2[v] = mp_set(g[(v << 1) + 1], p);
		}
		mp_NTT(logn, t1, gm, p, p0i);
		mp_NTT(logn, t2, gm, p, p0i);
		for (size_t v = 0; v < hn; v++) {
			uint32_t xe0 = t1[(v << 1) + 0];
			uint32_t xe1 = t1[(v << 1) + 1];
			uint32_t xo0 = t2[(v << 1) + 0];
			uint32_t xo1 = t2[(v << 1) + 1];
			uint32_t xv0 = gm[hn + v];
			uint32_t xv1 = p - xv0;
			xe0 = mp_montymul(xe0, xe0, p, p0i);
			xe1 = mp_montymul(xe1, xe1, p, p0i);
			xo0 = mp_montymul(xo0, xo0, p, p0i);
			xo1 = mp_montymul(xo1, xo1, p, p0i);
			uint32_t xg0 = mp_sub(xe0,
				mp_montymul(xo0, xv0, p, p0i), p);
			uint32_t xg1 = mp_sub(xe1,
				mp_montymul(xo1, xv1, p, p0i), p);

			uint32_t xkg0 = mp_montymul(
				mp_montymul(xg0, k32[(v << 1) + 0], p, p0i),
				R3, p, p0i);
			uint32_t xkg1 = mp_montymul(
				mp_montymul(xg1, k32[(v << 1) + 1], p, p0i),
				R3, p, p0i);
			Gu[(v << 1) + 0] = mp_sub(Gu[(v << 1) + 0], xkg0, p);
			Gu[(v << 1) + 1] = mp_sub(Gu[(v << 1) + 1], xkg1, p);
		}

		/*
		 * Convert back F and G to RNS.
		 */
		mp_mkigm(logn, t1, PRIMES[u].ig, p, p0i);
		mp_iNTT(logn, Fu, t1, p, p0i);
		mp_iNTT(logn, Gu, t1, p, p0i);

		/*
		 * We replaced k (plain 32-bit) with (2^sc)*k (NTT). We must
		 * put it back to its initial value if there should be another
		 * iteration.
		 */
		if ((u + 1) < FGlen) {
			mp_iNTT(logn, k32, t1, p, p0i);
			scv = (uint32_t)1 << (-sc & 31);
			for (uint32_t m = sc >> 5; m > 0; m--) {
				scv = mp_montymul(scv, 1, p, p0i);
			}
			for (size_t v = 0; v < n; v++) {
				k32[v] = (uint32_t)mp_norm(
					mp_montymul(scv, k32[v], p, p0i), p);
			}
		}
	}

	/*
	 * Output F and G are in RNS (non-NTT), but we want plain integers.
	 */
	if (FGlen == 1) {
		uint32_t p = PRIMES[0].p;
		for (size_t u = 0; u < n; u++) {
			F[u] = (uint32_t)mp_norm(F[u], p) & 0x7FFFFFFF;
			G[u] = (uint32_t)mp_norm(G[u], p) & 0x7FFFFFFF;
		}
	}
	else {
		uint32_t p0 = PRIMES[0].p;
		uint32_t p1 = PRIMES[1].p;
		uint32_t p1_0i = PRIMES[1].p0i;
		uint32_t s = PRIMES[1].s;
		uint64_t pp = (uint64_t)p0 * (uint64_t)p1;
		uint64_t hpp = pp >> 1;
		for (size_t u = 0; u < n; u++) {
			/*
			 * Apply CRT with two primes on the coefficient of F.
			 */
			uint32_t x0 = F[u];      /* mod p0 */
			uint32_t x1 = F[u + n];  /* mod p1 */
			uint32_t x0m1 = x0 - (p1 & ~tbmask(x0 - p1));
			uint32_t y = mp_montymul(
				mp_sub(x1, x0m1, p1), s, p1, p1_0i);
			uint64_t z = (uint64_t)x0 + (uint64_t)p0 * (uint64_t)y;
			z -= pp & -((hpp - z) >> 63);
			F[u] = (uint32_t)z & 0x7FFFFFFF;
			F[u + n] = (uint32_t)(z >> 31) & 0x7FFFFFFF;
		}
		for (size_t u = 0; u < n; u++) {
			/*
			 * Apply CRT with two primes on the coefficient of G.
			 */
			uint32_t x0 = G[u];      /* mod p0 */
			uint32_t x1 = G[u + n];  /* mod p1 */
			uint32_t x0m1 = x0 - (p1 & ~tbmask(x0 - p1));
			uint32_t y = mp_montymul(
				mp_sub(x1, x0m1, p1), s, p1, p1_0i);
			uint64_t z = (uint64_t)x0 + (uint64_t)p0 * (uint64_t)y;
			z -= pp & -((hpp - z) >> 63);
			G[u] = (uint32_t)z & 0x7FFFFFFF;
			G[u + n] = (uint32_t)(z >> 31) & 0x7FFFFFFF;
		}
	}
}


// Thomas Pornin. https://github.com/pornin/ntrugen
void
poly_sub_scaled_ntt(unsigned logn, uint32_t* restrict F, size_t Flen,
	const uint32_t* restrict f, size_t flen,
	//const int64_t* restrict k, 
	const int32_t* restrict k,
	uint32_t sc, uint32_t* restrict tmp)

{


	size_t n = (size_t)1 << logn;
	size_t tlen = flen + 1;
	uint32_t* gm = tmp;
	uint32_t* igm = gm + n;
	uint32_t* fk = igm + n;
	uint32_t* t1 = fk + (tlen << logn);
	uint32_t sch, scl;
	DIVREM31(sch, scl, sc);

	/*
	 * Compute k*f in fk[], in RNS notation.
	 * f is assumed to be already in RNS+NTT over flen+1 words.
	 */
	for (size_t u = 0; u < tlen; u++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		mp_mkgmigm(logn, gm, igm, PRIMES[u].g, PRIMES[u].ig, p, p0i);

		for (size_t v = 0; v < n; v++) {
			t1[v] = mp_set(k[v], p);
		}

		mp_NTT(logn, t1, gm, p, p0i);

		const uint32_t* fs = f + (u << logn);
		uint32_t* ff = fk + (u << logn);

		for (size_t v = 0; v < n; v++) {
			ff[v] = mp_montymul(
				mp_montymul(t1[v], fs[v], p, p0i), R2, p, p0i);
		}

		mp_iNTT(logn, ff, igm, p, p0i);
	}

	/*
	 * Rebuild k*f.
	 */
	zint_rebuild_CRT(fk, tlen, n, 1, 1, t1);

	/*
	 * Subtract k*f, scaled, from F.
	 */
	for (size_t u = 0; u < n; u++) {
		zint_sub_scaled(F + u, Flen, fk + u, tlen, n, sch, scl);
	}
}

// Thomas Pornin. https://github.com/pornin/ntrugen
void
poly_mp_set(unsigned logn, uint32_t* f, uint32_t p)
{
	size_t n = (size_t)1 << logn;

	for (size_t u = 0; u < n; u++) {
		uint32_t x = f[u];
		x |= (x & 0x40000000) << 1;
		f[u] = mp_set(*(int32_t*)&x, p);
	}
}
// Thomas Pornin. https://github.com/pornin/ntrugen
void
poly_mp_norm(unsigned logn, uint32_t* f, uint32_t p)
{
	size_t n = (size_t)1 << logn;

	for (size_t u = 0; u < n; u++) {
		f[u] = (uint32_t)mp_norm(f[u], p) & 0x7FFFFFFF;
	}
}

// Thomas Pornin. https://github.com/pornin/ntrugen
int
poly_big_to_small(unsigned logn, int8_t* restrict d,
	const uint32_t* restrict s, int lim)
{
	size_t n = (size_t)1 << logn;
	for (size_t u = 0; u < n; u++) {
		uint32_t x = s[u];
		x |= (x & 0x40000000) << 1;
		int32_t z = *(int32_t*)&x;
		if (z < -lim || z > lim) {
			return 0;
		}
		d[u] = (int8_t)z;
	}
	return 1;
}

// Thomas Pornin. https://github.com/pornin/ntrugen
// Changes:
//	to work with fixed point data use scale = 26

void
poly_mulselfadj_fft(fxr* a, unsigned logn)
{
	/*
	 * Since each coefficient is multiplied with its own conjugate,
	 * the result contains only real values.
	 */
	size_t n, hn, u;

	n = (size_t)1 << logn;
	hn = n >> 1;

	for (u = 0; u < hn; u++) {
		fxr a_re, a_im;

		a_re = a[u];
		a_im = a[u + hn];
		a[u] = fxr_add(fxr_sqr(a_re), fxr_sqr(a_im));
		a[u + hn] = fxr_zero;
	}

}

// Thomas Pornin. https://github.com/pornin/ntrugen
// Changes:
//	to work with fixed point data use scale = 26

void
poly_add(
	fxr* restrict a, const fxr* restrict b, unsigned logn)
{
	size_t n, u;

	n = (size_t)1 << logn;

	for (u = 0; u < n; u++) {
		a[u] = fxr_add(a[u], b[u]);
	}

}
// Thomas Pornin. https://github.com/pornin/ntrugen
// Changes:
//	to work with fixed point data use scale = 26
// Multiply complex number on adj next number
void
poly_muladj_fft(
	fxr* restrict a, const fxr* restrict b, unsigned logn)
{
	size_t n, hn, u;

	n = (size_t)1 << logn;
	hn = n >> 1;

	for (u = 0; u < hn; u++) {
		fxc ac, bc;
		ac.re = a[u];
		ac.im = a[u + hn];
		bc.re = b[u];
		bc.im.v = -b[u + hn].v;
		
		ac = fxc_mul(ac, bc);
		a[u] = ac.re;
		a[u + hn] = ac.im;
	}

}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
//	to work with fixed point data use scale = 26
void
poly_LDLmv_fft(
	fxr* restrict d11, fxr* restrict l10,
	const fxr* restrict g00, const fxr* restrict g01,
	const fxr* restrict g11, unsigned logn)
{

	size_t n, hn, u;

	//int success = 0;
	n = (size_t)1 << logn;
	hn = n >> 1;

	// 	for (u = 0; u < hn; u++) 
	for (u = 0; u < hn; u++)
	{

		fxc g00_, g01_, g11_, temp_;
		//fxr mu_re, mu_im;
		fxc mu_;

		g00_.re = g00[u];
		g00_.im = g00[u + hn];
		g01_.re = g01[u];
		g01_.im = g01[u + hn];
		temp_.re = g01_.re;
		temp_.im.v = -g01_.im.v;


		g11_.re = g11[u];
		g11_.im = g11[u + hn];
		

		
		//mu_ = fxc_div1_with_success(g01_, g00_, &success);
		mu_ = fxc_div1/*_with_success*/(g01_, g00_);
		
		//if (success == 0)
		{
			g01_ = fxc_mul(mu_, temp_);
			// FPC_SUB(d11[u], d11[u + hn], g11_re, g11_im, g01_re, g01_im);
			temp_ = fxc_sub(g11_, g01_);
			d11[u].v = temp_.re.v;
			d11[u + hn].v = temp_.im.v;
			l10[u] = mu_.re;
			l10[u + hn] = fxr_neg(mu_.im);
		}
	}
	//return success;

}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
// Function name poly_split_fft changed to poly_split1_fft
//	to work with fixed point data use scale = 26
//   Added logn = 2 handling 
void
poly_split1_fft(
	fxr* restrict f0, fxr* restrict f1,
	const fxr* restrict f, unsigned logn)
{

	/*
	 * The FFT representation we use is in bit-reversed order
	 * (element i contains f(w^(rev(i))), where rev() is the
	 * bit-reversal function over the ring degree. This changes
	 * indexes with regards to the Falcon specification.
	 */
	 
	if (logn > 2)
	{
		size_t n, hn, qn, u;

		n = (size_t)1 << logn;
		hn = n >> 1;
		qn = hn >> 1;


		/*
		 * We process complex values by pairs. For logn = 1, there is only
		 * one complex value (the other one is the implicit conjugate),
		 * so we add the two lines below because the loop will be
		 * skipped.
		 */
	

		 // for (u = 0; u < qn; u++) 
		for (u = 0; u < qn; u++)
		{
	

			fxc a, b;
			fxc t1, t2;

	
			a.re = f[(u << 1) + 0];
			a.im = f[(u << 1) + 0 + hn];
			b.re = f[(u << 1) + 1];
			b.im = f[(u << 1) + 1 + hn];
			// FPC_ADD(t_re, t_im, a_re, a_im, b_re, b_im);
			t1 = fxc_add(a, b);

	
			t1 = fxc_half(t1);
			f0[u].v = t1.re.v;
	
			f0[u + qn].v = t1.im.v;

	

			//FPC_SUB(t_re, t_im, a_re, a_im, b_re, b_im);
			t1 = fxc_sub(a, b);
	

	
			t2.re = GM_TAB[u + hn].re;
			t2.im.v = -GM_TAB[u + hn].im.v;
	
			t1 = fxc_mul(t1, t2);



	
			f1[u] = fxr_div2e(t1.re, 1);
			f1[u + qn] = fxr_div2e(t1.im, 1);
		}
	}
	else
	{
		if (logn == 2)
		{
			fxr a_re, a_im, b_re, b_im;
			fxr t_re, t_im;
			a_re = f[0];
			a_im = f[2];
			b_re = f[1];
			b_im = f[3];

	
			t_re = fxr_add(a_re, b_re);
			t_im = fxr_add(a_im, b_im);

			f0[0] = fxr_half(t_re);
			f0[1] = fxr_half(t_im);

	

			t_re = fxr_sub(a_re, b_re);
			t_im = fxr_sub(a_im, b_im);
			f1[0] = fxr_mul(fxr_add(t_re, t_im), fxr_invsqrt8);
			f1[1] = fxr_mul(fxr_sub(t_im, t_re), fxr_invsqrt8);
		}
		else
		{
			f0[0] = f[0];
			f1[0] = f[1];
		}
	}


}

// @author   Thomas Pornin <thomas.pornin@nccgroup.com>
// Changes:
// Function name poly_merge_fft changed to poly_merge1_fft
//	to work with fixed point data use scale = 26
//   Added logn = 2 handling 
void
poly_merge1_fft(
	fxr* restrict f,
	const fxr* restrict f0, const fxr* restrict f1, unsigned logn)
{
	
	if (logn > 2)
	{
		size_t n, hn, qn, u;
		fxc t1, t2;
		n = (size_t)1 << logn;
		hn = n >> 1;
		qn = hn >> 1;


		/*
		 * An extra copy to handle the special case logn = 1.
		 */

		 // for (u = 0; u < qn; u++) 
		for (u = 0; u < qn; u++)
		{
	
			fxc a, b;
			
			a.re = f0[u];
			a.im = f0[u + qn];

			t1.re = GM_TAB[((u + hn))].re;
			t1.im = GM_TAB[((u + hn))].im;
			t2.re = f1[u];
			t2.im = f1[u + qn];
			
			b = fxc_mul(t2, t1);



			
			t1 = fxc_add(a, b);
			

			f[(u << 1) + 0] = t1.re;
			f[(u << 1) + 0 + hn] = t1.im;
			
			t1 = fxc_sub(a, b);
			
			f[(u << 1) + 1] = t1.re;
			
			f[(u << 1) + 1 + hn] = t1.im;
		}
	}
	else
	{
		if (logn == 2)
		{
			fxr a_re, a_im, b_re, b_im;
			fxr t_re, t_im;

			a_re = f0[0];
			a_im = f0[1];
			b_re = f1[0];
			b_im = f1[1];

			

			t_re = fxr_mul(fxr_sub(b_re, b_im), fxr_invsqrt2);
			t_im = fxr_mul(fxr_add(b_re, b_im), fxr_invsqrt2);
			f[0] = fxr_add(a_re, t_re);
			f[2] = fxr_add(a_im, t_im);
			f[1] = fxr_sub(a_re, t_re);
			f[3] = fxr_sub(a_im, t_im);
		}
		else
		{
			f[0] = f0[0];
			f[1] = f1[0];
		}

	}

}

int poly_sub_scaled(unsigned logn,
	uint32_t* restrict F, size_t Flen,
	const uint32_t* restrict f, size_t flen,
	//const int64_t* restrict k, 
	const int32_t* restrict k,
	uint32_t sc)

{
	if (flen == 0) {
		return 0;
	}
	uint32_t sch, scl;
	DIVREM31(sch, scl, sc);
	if (sch >= Flen) {
		return 0;
	}
	F += (size_t)sch << logn;
	Flen -= sch;

	switch (logn) {
	case 1: {
		uint32_t t0 = 0;
		uint32_t t1 = 0;
		uint32_t signf0 = -(f[(flen << 1) - 2] >> 30) >> 1;
		uint32_t signf1 = -(f[(flen << 1) - 1] >> 30) >> 1;

		/*int64_t k0 = k[0];
		int64_t k1 = k[1];*/
		int32_t k0 = k[0];
		int32_t k1 = k[1];

		int64_t cc0 = 0;
		int64_t cc1 = 0;
		for (size_t u = 0; u < Flen; u++) {
			/*
			 * Next word, shifted.
			 */
			uint32_t f0, f1;
			if (u < flen) {
				f0 = f[(u << 1) + 0];
				f1 = f[(u << 1) + 1];
			}
			else {
				f0 = signf0;
				f1 = signf1;
			}
			uint32_t fs0 = ((f0 << scl) & 0x7FFFFFFF) | t0;
			uint32_t fs1 = ((f1 << scl) & 0x7FFFFFFF) | t1;
			t0 = f0 >> (31 - scl);
			t1 = f1 >> (31 - scl);

			uint32_t F0 = F[(u << 1) + 0];
			uint32_t F1 = F[(u << 1) + 1];
			int64_t z0 = (int64_t)F0 + cc0
				- (int64_t)fs0 * (int64_t)k0
				+ (int64_t)fs1 * (int64_t)k1;
			int64_t z1 = (int64_t)F1 + cc1
				- (int64_t)fs0 * (int64_t)k1
				- (int64_t)fs1 * (int64_t)k0;
			F[(u << 1) + 0] = (uint32_t)z0 & 0x7FFFFFFF;
			F[(u << 1) + 1] = (uint32_t)z1 & 0x7FFFFFFF;
			cc0 = (z0 >> 31)/* - carry00 + carry11*/;
			cc1 = (z1 >> 31)/* - carry01  + carry10*/;
		}
		return 0;
	}

	case 2: {

		uint32_t t0 = 0;
		uint32_t t1 = 0;
		uint32_t t2 = 0;
		uint32_t t3 = 0;
		uint32_t signf0 = -(f[(flen << 2) - 4] >> 30) >> 1;
		uint32_t signf1 = -(f[(flen << 2) - 3] >> 30) >> 1;
		uint32_t signf2 = -(f[(flen << 2) - 2] >> 30) >> 1;
		uint32_t signf3 = -(f[(flen << 2) - 1] >> 30) >> 1;
		int32_t k0 = k[0];
		int32_t k1 = k[1];
		int32_t k2 = k[2];
		int32_t k3 = k[3];
		int64_t cc0 = 0;
		int64_t cc1 = 0;
		int64_t cc2 = 0;
		int64_t cc3 = 0;
		for (size_t u = 0; u < Flen; u++) {
			/*
			 * Next word, shifted.
			 */
			uint32_t f0, f1, f2, f3;

			if (u < flen) {
				f0 = f[(u << 2) + 0];
				f1 = f[(u << 2) + 1];
				f2 = f[(u << 2) + 2];
				f3 = f[(u << 2) + 3];
			}
			else {
				f0 = signf0;
				f1 = signf1;
				f2 = signf2;
				f3 = signf3;
			}
			uint32_t fs0 = ((f0 << scl) & 0x7FFFFFFF) | t0;
			uint32_t fs1 = ((f1 << scl) & 0x7FFFFFFF) | t1;
			uint32_t fs2 = ((f2 << scl) & 0x7FFFFFFF) | t2;
			uint32_t fs3 = ((f3 << scl) & 0x7FFFFFFF) | t3;
			t0 = f0 >> (31 - scl);
			t1 = f1 >> (31 - scl);
			t2 = f2 >> (31 - scl);
			t3 = f3 >> (31 - scl);

			uint32_t F0 = F[(u << 2) + 0];
			uint32_t F1 = F[(u << 2) + 1];
			uint32_t F2 = F[(u << 2) + 2];
			uint32_t F3 = F[(u << 2) + 3];
			


			int64_t z0 = (int64_t)F0 + cc0
				- (int64_t)fs0 * (int64_t)k0
				//- res[0][0]
				+ (int64_t)fs1 * (int64_t)k3
				//+ res[1][3]
				+ (int64_t)fs2 * (int64_t)k2
				//+ res[2][2]
				+ (int64_t)fs3 * (int64_t)k1;
			//+ res[3][1];

			int64_t z1 = (int64_t)F1 + cc1
				- (int64_t)fs0 * (int64_t)k1
				//- res[0][1]
				- (int64_t)fs1 * (int64_t)k0
				//- res[1][0]
				+ (int64_t)fs2 * (int64_t)k3
				//+ res[2][3]
				+ (int64_t)fs3 * (int64_t)k2;
			//+ res[3][2];

			int64_t z2 = (int64_t)F2 + cc2
				- (int64_t)fs0 * (int64_t)k2
				//- res[0][2]
				- (int64_t)fs1 * (int64_t)k1
				//- res[1][1]
				- (int64_t)fs2 * (int64_t)k0
				//- res[2][0]
				+ (int64_t)fs3 * (int64_t)k3;
			//+ res[3][3];

			int64_t z3 = (int64_t)F3 + cc3
				- (int64_t)fs0 * (int64_t)k3
				//- res[0][3]
				- (int64_t)fs1 * (int64_t)k2
				//- res[1][2]
				- (int64_t)fs2 * (int64_t)k1
				//- res[2][1]
				- (int64_t)fs3 * (int64_t)k0;
			//- res[3][0];

			F[(u << 2) + 0] = (uint32_t)z0 & 0x7FFFFFFF;
			F[(u << 2) + 1] = (uint32_t)z1 & 0x7FFFFFFF;
			F[(u << 2) + 2] = (uint32_t)z2 & 0x7FFFFFFF;
			F[(u << 2) + 3] = (uint32_t)z3 & 0x7FFFFFFF;
			cc0 = (z0 >> 31) /*- carry[0][0] + carry[1][3] + carry[2][2] + carry[3][1]*/;
			cc1 = (z1 >> 31) /*- carry[0][1] - carry[1][0] + carry[2][3] + carry[3][2]*/;
			cc2 = (z2 >> 31) /*- carry[0][2] - carry[1][1] - carry[2][0] + carry[3][3]*/;
			cc3 = (z3 >> 31)/* - carry[0][3] - carry[1][2] - carry[2][1] - carry[3][0]*/;
		}

		return 0;
	}

	case 3: {

		uint32_t t0 = 0;
		uint32_t t1 = 0;
		uint32_t t2 = 0;
		uint32_t t3 = 0;
		uint32_t t4 = 0;
		uint32_t t5 = 0;
		uint32_t t6 = 0;
		uint32_t t7 = 0;
		uint32_t signf0 = -(f[(flen << 3) - 8] >> 30) >> 1;
		uint32_t signf1 = -(f[(flen << 3) - 7] >> 30) >> 1;
		uint32_t signf2 = -(f[(flen << 3) - 6] >> 30) >> 1;
		uint32_t signf3 = -(f[(flen << 3) - 5] >> 30) >> 1;
		uint32_t signf4 = -(f[(flen << 3) - 4] >> 30) >> 1;
		uint32_t signf5 = -(f[(flen << 3) - 3] >> 30) >> 1;
		uint32_t signf6 = -(f[(flen << 3) - 2] >> 30) >> 1;
		uint32_t signf7 = -(f[(flen << 3) - 1] >> 30) >> 1;
		int32_t k0 = k[0];
		int32_t k1 = k[1];
		int32_t k2 = k[2];
		int32_t k3 = k[3];
		int32_t k4 = k[4];
		int32_t k5 = k[5];
		int32_t k6 = k[6];
		int32_t k7 = k[7];
		int64_t cc0 = 0;
		int64_t cc1 = 0;
		int64_t cc2 = 0;
		int64_t cc3 = 0;
		int64_t cc4 = 0;
		int64_t cc5 = 0;
		int64_t cc6 = 0;
		int64_t cc7 = 0;
		for (size_t u = 0; u < Flen; u++) {
			/*
			 * Next word, shifted.
			 */
			uint32_t f0, f1, f2, f3, f4, f5, f6, f7;
			uint32_t fs0, fs1, fs2, fs3, fs4, fs5, fs6, fs7;
			//uint32_t fs[8];
			if (u < flen) {
				f0 = f[(u << 3) + 0];
				f1 = f[(u << 3) + 1];
				f2 = f[(u << 3) + 2];
				f3 = f[(u << 3) + 3];
				f4 = f[(u << 3) + 4];
				f5 = f[(u << 3) + 5];
				f6 = f[(u << 3) + 6];
				f7 = f[(u << 3) + 7];
			}
			else {
				f0 = signf0;
				f1 = signf1;
				f2 = signf2;
				f3 = signf3;
				f4 = signf4;
				f5 = signf5;
				f6 = signf6;
				f7 = signf7;
			}

			

			fs0 = ((f0 << scl) & 0x7FFFFFFF) | t0;
			fs1 = ((f1 << scl) & 0x7FFFFFFF) | t1;
			fs2 = ((f2 << scl) & 0x7FFFFFFF) | t2;
			fs3 = ((f3 << scl) & 0x7FFFFFFF) | t3;
			fs4 = ((f4 << scl) & 0x7FFFFFFF) | t4;
			fs5 = ((f5 << scl) & 0x7FFFFFFF) | t5;
			fs6 = ((f6 << scl) & 0x7FFFFFFF) | t6;
			fs7 = ((f7 << scl) & 0x7FFFFFFF) | t7;

			t0 = f0 >> (31 - scl);
			t1 = f1 >> (31 - scl);
			t2 = f2 >> (31 - scl);
			t3 = f3 >> (31 - scl);
			t4 = f4 >> (31 - scl);
			t5 = f5 >> (31 - scl);
			t6 = f6 >> (31 - scl);
			t7 = f7 >> (31 - scl);
			

			uint32_t F0 = F[(u << 3) + 0];
			uint32_t F1 = F[(u << 3) + 1];
			uint32_t F2 = F[(u << 3) + 2];
			uint32_t F3 = F[(u << 3) + 3];
			uint32_t F4 = F[(u << 3) + 4];
			uint32_t F5 = F[(u << 3) + 5];
			uint32_t F6 = F[(u << 3) + 6];
			uint32_t F7 = F[(u << 3) + 7];
			int64_t z0 = (int64_t)F0 + cc0
				- (int64_t)fs0 * (int64_t)k0
				//- res[0][0]
				+ (int64_t)fs1 * (int64_t)k7
				//+ res[1][7]
				+ (int64_t)fs2 * (int64_t)k6
				//+ res[2][6]
				+ (int64_t)fs3 * (int64_t)k5
				//+ res[3][5]
				+ (int64_t)fs4 * (int64_t)k4
				//+ res[4][4]
				+ (int64_t)fs5 * (int64_t)k3
				//+ res[5][3]
				+ (int64_t)fs6 * (int64_t)k2
				//+ res[6][2]

				+ (int64_t)fs7 * (int64_t)k1;
			//+ res[7][1];
			int64_t z1 = (int64_t)F1 + cc1
				- (int64_t)fs0 * (int64_t)k1
				//- res[0][1]
				- (int64_t)fs1 * (int64_t)k0
				//- res[1][0]
				+ (int64_t)fs2 * (int64_t)k7
				//+ res[2][7]
				+ (int64_t)fs3 * (int64_t)k6
				//+ res[3][6]
				+ (int64_t)fs4 * (int64_t)k5
				//+ res[4][5]
				+ (int64_t)fs5 * (int64_t)k4
				//+ res[5][4]
				+ (int64_t)fs6 * (int64_t)k3
				//+ res[6][3]
				+ (int64_t)fs7 * (int64_t)k2;
			//+ res[7][2];
			int64_t z2 = (int64_t)F2 + cc2
				- (int64_t)fs0 * (int64_t)k2
				//- res[0][2]
				- (int64_t)fs1 * (int64_t)k1
				//- res[1][1]
				- (int64_t)fs2 * (int64_t)k0
				//- res[2][0]
				+ (int64_t)fs3 * (int64_t)k7
				//+ res[3][7]
				+ (int64_t)fs4 * (int64_t)k6
				//+ res[4][6]
				+ (int64_t)fs5 * (int64_t)k5
				//+ res[5][5]
				+ (int64_t)fs6 * (int64_t)k4
				//+ res[6][4]

				+ (int64_t)fs7 * (int64_t)k3;
			//+ res[7][3];

			int64_t z3 = (int64_t)F3 + cc3
				- (int64_t)fs0 * (int64_t)k3
				//res[0][3]
				- (int64_t)fs1 * (int64_t)k2
				//res[1][2]
				- (int64_t)fs2 * (int64_t)k1
				//- res[2][1]
				- (int64_t)fs3 * (int64_t)k0
				//- res[3][0]
				+ (int64_t)fs4 * (int64_t)k7
				//+ res[4][7]
				+ (int64_t)fs5 * (int64_t)k6
				//+ res[5][6]
				+ (int64_t)fs6 * (int64_t)k5
				//+ res[6][5]
				+ (int64_t)fs7 * (int64_t)k4;
			//+ res[7][4];

			int64_t z4 = (int64_t)F4 + cc4
				- (int64_t)fs0 * (int64_t)k4
				//- res[0][4]
				- (int64_t)fs1 * (int64_t)k3
				//- res[1][3]
				- (int64_t)fs2 * (int64_t)k2
				//- res[2][2]
				- (int64_t)fs3 * (int64_t)k1
				//- res[3][1]
				- (int64_t)fs4 * (int64_t)k0
				//- res[4][0]
				+ (int64_t)fs5 * (int64_t)k7
				//+ res[5][7]
				+ (int64_t)fs6 * (int64_t)k6
				//+ res[6][6]
				+ (int64_t)fs7 * (int64_t)k5;
			//+ res[7][5];

			int64_t z5 = (int64_t)F5 + cc5
				- (int64_t)fs0 * (int64_t)k5
				//- res[0][5]
				- (int64_t)fs1 * (int64_t)k4
				//- res[1][4]
				- (int64_t)fs2 * (int64_t)k3
				//- res[2][3]
				- (int64_t)fs3 * (int64_t)k2
				//- res[3][2]
				- (int64_t)fs4 * (int64_t)k1
				//- res[4][1]
				- (int64_t)fs5 * (int64_t)k0
				//- res[5][0]
				+ (int64_t)fs6 * (int64_t)k7
				//+ res[6][7]
				+ (int64_t)fs7 * (int64_t)k6;
			//+ res[7][6];

			int64_t z6 = (int64_t)F6 + cc6
				- (int64_t)fs0 * (int64_t)k6
				//- res[0][6]
				- (int64_t)fs1 * (int64_t)k5
				//- res[1][5]
				- (int64_t)fs2 * (int64_t)k4
				//- res[2][4]
				- (int64_t)fs3 * (int64_t)k3
				//- res[3][3]
				- (int64_t)fs4 * (int64_t)k2
				//- res[4][2]
				- (int64_t)fs5 * (int64_t)k1
				//- res[5][1]
				- (int64_t)fs6 * (int64_t)k0
				//- res[6][0]
				+ (int64_t)fs7 * (int64_t)k7;
			//+ res[7][7];

			int64_t z7 = (int64_t)F7 + cc7
				- (int64_t)fs0 * (int64_t)k7
				//- res[0][7]
				- (int64_t)fs1 * (int64_t)k6
				//- res[1][6]
				- (int64_t)fs2 * (int64_t)k5
				//- res[2][5]
				- (int64_t)fs3 * (int64_t)k4
				//- res[3][4]
				- (int64_t)fs4 * (int64_t)k3
				//- res[4][3]
				- (int64_t)fs5 * (int64_t)k2
				//- res[5][2]
				- (int64_t)fs6 * (int64_t)k1
				//- res[6][1]
				- (int64_t)fs7 * (int64_t)k0;
			//- res[7][0];

			F[(u << 3) + 0] = (uint32_t)z0 & 0x7FFFFFFF;
			F[(u << 3) + 1] = (uint32_t)z1 & 0x7FFFFFFF;
			F[(u << 3) + 2] = (uint32_t)z2 & 0x7FFFFFFF;
			F[(u << 3) + 3] = (uint32_t)z3 & 0x7FFFFFFF;
			F[(u << 3) + 4] = (uint32_t)z4 & 0x7FFFFFFF;
			F[(u << 3) + 5] = (uint32_t)z5 & 0x7FFFFFFF;
			F[(u << 3) + 6] = (uint32_t)z6 & 0x7FFFFFFF;
			F[(u << 3) + 7] = (uint32_t)z7 & 0x7FFFFFFF;
			cc0 = (z0 >> 31) /*- carry[0][0] + carry [1][7] + carry[2][6] + carry[3][5] +
				carry [4][4] + carry[5][3] + carry[6][2] + carry[7][1]*/;
			cc1 = (z1 >> 31) /*- carry[0][1] - carry[1][0] + carry[2][7] + carry[3][6] +
				carry[4][5] + carry[5][4] + carry[6][3] + carry[7][2]*/;
			cc2 = (z2 >> 31)/* - carry[0][2] - carry[1][1] - carry[2][0] + carry[3][7] +
				carry[4][6] + carry[5][5] + carry[6][4] + carry[7][3]*/;
			cc3 = (z3 >> 31) /*- carry[0][3] - carry[1][2] - carry[2][1] + carry[3][0] +
				carry[4][7] + carry[5][6] + carry[6][5] + carry[7][4]*/;
			cc4 = (z4 >> 31) /*- carry[0][4] - carry[1][3] - carry[2][2] - carry[3][1] -
				carry[4][0] + carry[5][7] + carry[6][6] + carry[7][5]*/;
			cc5 = (z5 >> 31)/* - carry[0][5] - carry[1][4] - carry[2][3] - carry[3][2] -
				carry[4][1] - carry[5][0] + carry[6][7] + carry[7][6]*/;
			cc6 = (z6 >> 31) /*- carry[0][6] - carry[1][5] - carry[2][4] - carry[3][3] -
				carry[4][2] - carry[5][1] - carry[6][0] + carry[7][7]*/;
			cc7 = (z7 >> 31) /*- carry[0][7] - carry[1][6] - carry[2][5] - carry[3][4] -
				carry[4][3] - carry[5][2] - carry[6][1] - carry[7][0]*/;
		}

		return 0;
	}
	}
	size_t n = (size_t)1 << logn;
	for (size_t u = 0; u < n; u++) {
		int32_t kf = -k[u];
		uint32_t* x = F + u;
		for (size_t v = 0; v < n; v++) {
			zint_add_scaled_mul_small(
				x, Flen, f + v, flen, n, kf, 0, scl);
			if (u + v == n - 1) {
				x = F;
				kf = -kf;
			}
			else {
				x++;
			}
		}
	}
	return 0;
}

void
poly_big_to_fixed(unsigned logn, fxr* restrict d, const uint32_t* restrict f,
	size_t len, uint32_t sc)
{
	size_t n = (size_t)1 << logn;
	if (len == 0) {
		memset(d, 0, n * sizeof * d);
		return;
	}

	/*
	 * We split the bit length into sch and scl such that:
	 *   sc = 31*sch + scl
	 * We also want scl in the 1..31 range, not 0..30. It may happen
	 * that sch becomes -1, which will "wrap around" (harmlessly).
	 *
	 * For each coefficient, we need three words, each with a given
	 * left shift (negative for a right shift):
	 *    sch-1   1 - scl
	 *    sch     32 - scl
	 *    sch+1   63 - scl
	 */
	uint32_t sch, scl;


	DIVREM31(sch, scl, sc);


	uint32_t z = (scl - 1) >> 31;
	sch -= z;
	scl |= 31 & -z;

	uint32_t t0 = (uint32_t)(sch - 1) & 0xFFFFFF;
	uint32_t t1 = sch & 0xFFFFFF;
	uint32_t t2 = (uint32_t)(sch + 1) & 0xFFFFFF;

	for (size_t u = 0; u < n; u++, f++) {
		uint32_t w0, w1, w2, ws, xl, xh;

		w0 = 0;
		w1 = 0;
		w2 = 0;
		uint32_t t, w;
		for (size_t v = 0; v < len; v++) {


			w = f[v << logn];
			t = (uint32_t)v & 0xFFFFFF;
			w0 |= w & -((uint32_t)((t ^ t0) - 1) >> 31);
			w1 |= w & -((uint32_t)((t ^ t1) - 1) >> 31);
			w2 |= w & -((uint32_t)((t ^ t2) - 1) >> 31);
		}


		/*
		 * If there were not enough words for the requested
		 * scaling, then we must supply copies with the proper
		 * sign.
		 */
		ws = -(f[(len - 1) << logn] >> 30) >> 1;
		w0 |= ws & -((uint32_t)((uint32_t)len - sch) >> 31);
		w1 |= ws & -((uint32_t)((uint32_t)len - sch - 1) >> 31);
		w2 |= ws & -((uint32_t)((uint32_t)len - sch - 2) >> 31);


		/*
		 * Assemble the 64-bit value with the shifts. We assume
		 * that shifts on 32-bit values are constant-time with
		 * regard to the shift count (this should be true on all
		 * modern architectures; the last notable arch on which
		 * shift timing depended on the count was the Pentium IV).
		 *
		 * Since the shift count (scl) is guaranteed to be in 1..31,
		 * we do not have special cases to handle.
		 *
		 * We must sign-extend w2 to ensure the sign bit is properly
		 * set in the fnr value.
		 */
		w2 |= (uint32_t)(w2 & 0x40000000) << 1;
		xl = (w0 >> (scl - 1)) | (w1 << (32 - scl));
		xh = (w1 >> scl) | (w2 << (31 - scl));



		int64_t temp = (int64_t)(int32_t)xh;
		d[u].v = ((uint64_t)xl >> 6) | ((int64_t)temp << (32 - 6));


	}

}


