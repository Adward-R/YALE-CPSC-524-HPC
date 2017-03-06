# Set the compiler command
CC = mpicc

# Set the compiler options
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

# The first target ("all") is built if you just run "make".
# "all" depends on hello, so to build "all", make will first build "hello"
all: task1 

# This target ("task1") is your program. You can build it using "make task1".
# "task1" depends on 3 object files. You have been given the last two, but you
# need to create task1.o.
task1: task1.o timing.o rwork.o

# The following command is used to create the executable "hello" once you have all 3 dependencies.
# $(CC) is replaced by the value of CC (mpicc, in this case). $@ is replaced by the target name.
# $(CFLAGS) is replaced by the value of CFLAGS. $^ is replaced by the string of dependencies.
# IMPORTANT NOTE: Indented lines like this MUST begin with a TAB, not a bunch of spaces!
	$(CC) -o $@ $(CFLAGS) $^

# This is a special type of rule that tells make how to make a .o file form a .c file. Make will use
# it to create task1.o from task1.c (assuming it can find task1.c). $< is replaced by the name of
# the .c file.
.c.o:
	$(CC) $(CFLAGS) -c $<

# This command is just another target that is often included to clean up so that you can build 
# everything from scratch. It simply deletes the file that are built by other parts of the Makefile
clean:
	rm -f task1 task1.o

