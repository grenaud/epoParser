
CXX      = g++   
LIBGAB   = /net/node07/storage/ctools/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/  -c
LDFLAGS  = -lz


all: epoParser ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o



ParseEntireEPOBlock.o:	${LIBGAB}libgab.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlock.cpp

ParseEntireEPOBlockRelax.o:	${LIBGAB}libgab.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlockRelax.cpp

epoParser.o:	epoParser.cpp
	${CXX} ${CXXFLAGS} epoParser.cpp


epoParser:	epoParser.o ${LIBGAB}libgab.o  ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 



clean :
	rm -f epoParser.o epoParser ParseEntireEPOBlock.o  ParseEntireEPOBlockRelax.o

