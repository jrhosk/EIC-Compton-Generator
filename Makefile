#
#  Makefile to compile generator code
#    
#
#

CXX=g++
LD=g++
CXXFLAGS=-g -O2
AC_CXXFLAGS= -I/u/apps/root/5.34.21/root/include
AC_LDFLAGS= -L/u/apps/root/5.34.21/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -L/u/apps/root/5.34.21/root/lib -lGui

all: generator

%.o: %.cc
	${CXX} ${CXXFLAGS} -Iinclude -I${AC_CXXFLAGS} -c -o $@ $< 

generator: generator.o src/Generator.o
	$(CXX) -Iinclude -I${AC_CXXFFLAGS} -o $@ ${CXXFLAGS} $^ ${AC_LDFLAGS}

.PHONY: clean
clean:
	rm -f generator *.o
