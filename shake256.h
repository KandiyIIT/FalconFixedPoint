#ifndef SHAKE256__H
#define SHAKE256__H

#include "inner.h"

typedef struct {
	uint64_t opaque_contents[26];
}shake256_context;

//typedef struct {
//	union {
//		uint64_t A[25];
//		uint8_t dbuf[200];
//	} st;
//	uint64_t dptr;
//} inner_shake256_context;

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

//int
//shake256_init_prng_from_system(shake256_context* sc);

//typedef struct {
//	union {
//		uint8_t d[512]; /* MUST be 512, exactly */
//		uint64_t dummy_u64;
//	} buf;
//	size_t ptr;
//	union {
//		uint8_t d[256];
//		uint64_t dummy_u64;
//	} state;
//	int type;
//} prng;


void prng_init(prng* p, inner_shake256_context* src);
void prng_refill(prng* p);
void prng_get_bytes(prng* p, void* dst, size_t len);
static inline uint64_t
prng_get_u64(prng* p)
{
	size_t u;

	/*
	 * If there are less than 9 bytes in the buffer, we refill it.
	 * This means that we may drop the last few bytes, but this allows
	 * for faster extraction code. Also, it means that we never leave
	 * an empty buffer.
	 */
	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9) {
		prng_refill(p);
		u = 0;
	}
	p->ptr = u + 8;

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
#if FALCON_LE && FALCON_UNALIGNED  // yyyLEU+1
	return *(uint64_t*)(p->buf.d + u);
#else  // yyyLEU+0
	return (uint64_t)p->buf.d[u + 0]
		| ((uint64_t)p->buf.d[u + 1] << 8)
		| ((uint64_t)p->buf.d[u + 2] << 16)
		| ((uint64_t)p->buf.d[u + 3] << 24)
		| ((uint64_t)p->buf.d[u + 4] << 32)
		| ((uint64_t)p->buf.d[u + 5] << 40)
		| ((uint64_t)p->buf.d[u + 6] << 48)
		| ((uint64_t)p->buf.d[u + 7] << 56);
#endif  // yyyLEU-
}


static inline unsigned
prng_get_u8(prng* p)
{
	unsigned v;

	v = p->buf.d[p->ptr++];
	if (p->ptr == sizeof p->buf.d) {
		prng_refill(p);
	}
	return v;
}


#endif