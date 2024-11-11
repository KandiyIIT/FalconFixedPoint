#ifndef __API__H
#define __API__H

int
expand_privkey(fxr* restrict expanded_key,
	const int8_t* f, const int8_t* g,
	const int8_t* F, const int8_t* G,
	unsigned logn, uint8_t* restrict tmp);

#endif // __API__H
