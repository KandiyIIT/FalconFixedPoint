#ifndef NTT__H
#define NTT__H

inline uint32_t
mq_conv_small(int x);

void
mq_NTT(uint16_t* a, unsigned logn);
void
mq_iNTT(uint16_t* a, unsigned logn);

void
mq_poly_tomonty(uint16_t* f, unsigned logn);

#endif // NTT__H
