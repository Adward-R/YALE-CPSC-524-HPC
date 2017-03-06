CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

task4:	task4.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f task4 *.o
