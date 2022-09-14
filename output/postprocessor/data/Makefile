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
vpath %.h ./headers
vpath %.h /usr/lib/gcc-lib/x86_64-redhat-linux/3.3.3/include/
vpath %.h /usr/include/
vpath %.h /usr/local/include/

# To generate data files 

#datafiles: datamaker.o datafile_generator.o stdlib.h stdio.h
#	gcc datamaker.o datafile_generator.o

#datamaker.o: datamaker.c datafile_generator.c stdlib.h stdio.h
#	gcc -c $(COMPOPS) $<

#datafile_generator.o: datafile_generator.c stdio.h stdlib.h
#	gcc -c $(COMPOPS) $<

datafiles: datamakere.o datamakerc.o datafile_generator.o stdlib.h stdio.h
	 gcc -o eta.out datamakere.o datafile_generator.o
	 gcc -o comp.out datamakerc.o datafile_generator.o
 	 
datamakerc.o: datamakerc.c datafile_generator.c stdlib.h stdio.h
	gcc -c $(COMPOPS) $<

datamakere.o: datamakere.c datafile_generator.c stdlib.h stdio.h
	gcc -c $(COMPOPS) $<

datafile_generator.o: datafile_generator.c stdio.h stdlib.h
	gcc -c $(COMPOPS) $<


# To remove the intermediate files

clean:
	rm -rf *.o
	rm -rf *.bbg 
	rm -rf *.bb 
	rm -rf *.da 
	rm -rf *.gcov 
	rm -rf ../../../postprocess/data/*.data
	rm -rf ../../../postprocess/data/*.dat

# To remove all intermediate files

CLEAN:
	rm -rf a.out *.o

# End of the makefile
