CC = icpc

CFLAGS= -g -O3 -xHost -mkl -fno-alias -qopenmp

OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt 

TARGETS = serial parallel1 parallel2 parallel3
TARGETOBJECTS = serial.o parallel1.o parallel2.o parallel3.o

.SUFFIXES: .o .cpp

# all: $(TARGETS)

# $(TARGETS): $(TARGETOBJECTS)
# 	$(CC) -o $@ $(CFLAGS) $^


serial: serial.o
	$(CC) -o $@ $(CFLAGS) $^
parallel1: parallel1.o
	$(CC) -o $@ $(CFLAGS) $^
parallel2: parallel2.o
	$(CC) -o $@ $(CFLAGS) $^
parallel3: parallel3.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $(OPTFLAGS) -o $@ $<

clean: 
	rm -f $(TARGETOBJECTS) $(TARGETS) *.optrpt
