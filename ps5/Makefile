CC = icpc

CFLAGS= -g -O3 -xHost -mkl -fno-alias -qopenmp
#CFLAGS= -g -O0 -xHost -mkl -fno-alias -qopenmp

OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt 

TARGETS = nbody2
TARGETOBJECTS = nbody2.o

NSIZE=16384

.SUFFIXES: .o .c

all: $(TARGETS)


$(TARGETS): $(TARGETOBJECTS)
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $(OPTFLAGS) -o $@ $<

clean: 
	rm -f $(TARGETOBJECTS) $(TARGETS) *.optrpt
