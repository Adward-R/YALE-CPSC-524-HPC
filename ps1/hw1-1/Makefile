IC = icc
IFLAGS1 = -g -O0 -fno-alias -std=c99
IFLAGS2 = -g -O1 -fno-alias -std=c99
IFLAGS3 = -g -O3 -no-vec -no-simd -fno-alias -std=c99
IFLAGS4 = -g -O3 -xHost -fno-alias -std=c99

pi_1 pi_2 pi_3 pi_4: pi.c timing.o
	${IC} ${IFLAGS1} -o pi_1 pi.c timing.o
	${IC} ${IFLAGS2} -o pi_2 pi.c timing.o
	${IC} ${IFLAGS3} -o pi_3 pi.c timing.o
	${IC} ${IFLAGS4} -o pi_4 pi.c timing.o
