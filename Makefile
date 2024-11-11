.POSIX:
CC = gcc
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -march=native
LD = gcc 
LDFLAGS = 
LIBS = -lm

OBJ = fxp_scale.o keygen.o mp.o mp31.o mq_12289.o ntru.o ntt.o poly.o rng.o shake.o sign.o
all: test_falcon

clean:
	-rm -f $(OBJ) main.o test_falcon

test_falcon: main.o $(OBJ)
	$(LD) $(LDFLAGS) -o test_falcon main.o $(OBJ) $(LIBS)
	
fxp_scale.o:
	$(CC) $(CFLAGS) -c -o fxp_scale.o fxp_scale.c
keygen.o:
	$(CC) $(CFLAGS) -c -o keygen.o keygen.c
mp.o:
	$(CC) $(CFLAGS) -c -o mp.o mp.c
mp31.o:
	$(CC) $(CFLAGS) -c -o mp31.o mp31.c
mq_12289.o:
	$(CC) $(CFLAGS) -c -o mq_12289.o mq_12289.c
ntru.o:
	$(CC) $(CFLAGS) -c -o ntru.o ntru.c
ntt.o:
	$(CC) $(CFLAGS) -c -o ntt.o ntt.c
poly.o:
	$(CC) $(CFLAGS) -c -o poly.o poly.c
rng.o:
	$(CC) $(CFLAGS) -c -o rng.o rng.c
shake.o:
	$(CC) $(CFLAGS) -c -o shake.o shake.c
#shake256.o:
#	$(CC) $(CFLAGS) -c -o shake256.o shake256.c
sign.o:
	$(CC) $(CFLAGS) -c -o sign.o sign.c
