CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

task3:	task3.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f task3 *.o
