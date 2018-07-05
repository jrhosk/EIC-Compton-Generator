#
#  Makefile to compile generator code
#    
#
#

CXX=g++
LD=g++
CXXFLAGS=-g -O2 -O -Wall
INCLUDE= include
SRC= src
ROOTLIBS=-L/u/apps/root/5.34.21/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -L/u/apps/root/5.34.21/root/lib -lGui
ROOTINC=/u/apps/root/5.34.21/root/include

all: generator

%.o: %.cc
	${CXX} ${CXXFLAGS} -I${INCLUDE} -I${ROOTINC} -c -o $@ $< 

generator: generator.o ${SRC}/Generator.o
	$(CXX) -I${INCLUDE} -I${ROOTINC} -o $@ ${CXXFLAGS} $^ ${ROOTLIBS}

.PHONY: clean
clean:
	rm -f generator *.o
