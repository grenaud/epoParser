CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS = -lm -O3 -lz -I${LIBGAB}  -c
LDFLAGS  = -lz


all: epoParser ParseEntireEPOBlock.o



ParseEntireEPOBlock.o:	${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlock.cpp

epoParser.o:	epoParser.cpp
	${CXX} ${CXXFLAGS} epoParser.cpp




epoParser:	epoParser.o ${LIBGAB}utils.o  gzstream/gzstream.o ParseEntireEPOBlock.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)



clean :
	rm -f epoParser.o epoParser ParseEntireEPOBlock.o

