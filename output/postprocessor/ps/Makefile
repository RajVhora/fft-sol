# Author - MPGuruajan. Last update 22.12.03

# According to make manual, every makefile is supposed to contain this
# line

SHELL = /bin/sh

# Suffi(ces?)-xes of files used in this makefile 

.SUFFIXES: .c .h .o

# Compilation and linking options

#COMPOPS = -fprofile-arcs -ftest-coverage -g -Wall
COMPOPS = -g -Wall
LINKOPS = -lm

# Directory paths for the source and header files 

vpath %.c ./src
vpath %.h /usr/lib/gcc-lib/x86_64-redhat-linux/3.3.3/include/
vpath %.h /usr/include/
vpath %.h /usr/local/include/

# To generate ps files 
psfiles: psmakere.o psmakerc.o psfile_generator.o stdlib.h stdio.h
	 gcc -o eta.out psmakere.o psfile_generator.o
	 gcc -o comp.out psmakerc.o psfile_generator.o
 	 
psmakerc.o: psmakerc.c psfile_generator.c stdlib.h stdio.h
	gcc -c $(COMPOPS) $<

psmakere.o: psmakere.c psfile_generator.c stdlib.h stdio.h
	gcc -c $(COMPOPS) $<

psfile_generator.o: psfile_generator.c stdio.h stdlib.h
	gcc -c $(COMPOPS) $<

# To remove the intermediate files

clean:
	rm -rf *.o
	rm -rf *.bbg 
	rm -rf *.bb 
	rm -rf *.da 
	rm -rf *.gcov 
	rm -rf ../../../postprocess/ps/*.ps

# To remove all intermediate files

CLEAN:
	rm -rf comp.out *.o
	rm -rf eta.out *.o
	rm -rf a.out *.o

# End of the makefile
