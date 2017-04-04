IC = icc
CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

serial: serial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(IC) -o $@ $(CFLAGS) $^

fserial: fserial.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f serial fserial *.o