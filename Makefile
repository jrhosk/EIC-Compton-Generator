# Makefile to compile compton generator                                                                                                                    
#       Joshua hoskins                                                                                                                                    
#         October 2016                                                                                                                                       
#                                                                                                                                                         

ROOTLIBS   = $(shell root-config --libs ) -lSpectrum
ROOTGLIBS  = $(shell root-config --glibs)
INCLUDES   = -I$(shell root-config --incdir) -Iinclude/ 
CC         = g++ ${INCLUDES}
SRC        = src
CFLAGS     = -O -Wall ${INCLUDES}

all: generator

%.o: %.cc
	${CC} ${CFLAGS} ${EXTRAFLAGS} -c -o $@ $< 
generator : generator.o ${SRC}/Generator.o
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${ROOTLIBS} ${ROOTGLIBS}
clean:
	rm -f *.o ~* 
