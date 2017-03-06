CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

task2:	task2.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f task2 *.o
