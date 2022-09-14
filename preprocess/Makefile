# Author - MPGuruajan. Last update 22.12.03

# According to make manual, every makefile is supposed to contain this
# line

SHELL = /bin/sh

#
#
#
# To change the random number seed: GSL_RNG_SEED=123 ./a.out
#
#
#
#


# Suffi(ces?)-xes of files used in this makefile 

.SUFFIXES: .c .h .o

# Compilation and linking options

#COMPOPS = -fprofile-arcs -ftest-coverage -g -Wall
COMPOPS = -g -Wall
LINKOPS = -static -L/usr/local/lib -lgsl -lgslcblas -lm -g

# List of header files

init_headers = stdio.h stdlib.h math.h init_functions.h \
nrutil.h

# List of source files

init_sources = muse_init.c nrutil.c add_noise.c random_ppt.c \
spinodal.c circular.c elliptic.c initial_configuration.c \
thinfilm.c twocirc.c

# List of object files 

init_objects = INIT.o nrutil.o initial_configuration.o add_noise.o \
random_ppt.o spinodal.o circular.o elliptic.o thinfilm.o twocirc.o

# Directory paths for the source and header files 

vpath %.c ./src
vpath %.h ./headers
#vpath %.h /usr/lib/gcc-lib/x86_64-redhat-linux/3.3.3/include/
#vpath %.h /usr/lib/gcc-lib/i386-redhat-linux/3.3.3/include/
vpath %.h /usr/lib/gcc/i686-apple-linux-darwin8/4.0.1/include/
vpath %.h /usr/include/
vpath %.h /usr/local/include/
vpath %.h /usr/local/include/gsl/
vpath %.h /usr/include/gsl/

# Compile and link

init: $(init_objects) $(init_headers)
	gcc $(init_objects) $(LINKOPS)

INIT.o: $(init_sources) $(init_headers)
	gcc -o $@ -c ./src/muse_init.c $(COMPOPS)

spinodal.o: spinodal.c
	gcc -c $(COMPOPS) $<

circular.o: circular.c math.h
	gcc -c $(COMPOPS) $<

elliptic.o: elliptic.c math.h
	gcc -c $(COMPOPS) $<

twocirc.o: twocirc.c math.h
	gcc -c $(COMPOPS) $<

random_ppt.o: random_ppt.c math.h stdlib.h gsl_rng.h gsl_math.h
	gcc -c $(COMPOPS) $<

thinfilm.o: thinfilm.c math.h
	gcc -c $(COMPOPS) $<

initial_configuration.o: initial_configuration.c random_ppt.c \
circular.c spinodal.c twocirc.c stdio.h stdlib.h math.h
	gcc -c $(COMPOPS) $<

add_noise.o: add_noise.c gsl_rng.h gsl_math.h \
stdio.h stdlib.h
	gcc -c $(COMPOPS) $<

# To remove the intermediate files

clean:
	rm -rf *.o
	rm -rf *.bbg 
	rm -rf *.bb 
	rm -rf *.da 
	rm -rf *.gcov 

# To remove all intermediate files

CLEAN:
	rm -rf a.out *.o

# End of the makefile
