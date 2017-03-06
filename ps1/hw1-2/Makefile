IC = icc
IFLAGS = -g -O3 -xHost -fno-alias -std=c99

triad: triad.c dummy.c timing.o
	${IC} ${IFLAGS} -o triad triad.c dummy.c timing.o
