
CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/VCFparser/gzstream/  -c
LDFLAGS  = -lz


all: epoParser ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o



ParseEntireEPOBlock.o:	${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlock.cpp

ParseEntireEPOBlockRelax.o:	${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlockRelax.cpp

epoParser.o:	epoParser.cpp
	${CXX} ${CXXFLAGS} epoParser.cpp


epoParser:	epoParser.o ${LIBGAB}utils.o  ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 



clean :
	rm -f epoParser.o epoParser ParseEntireEPOBlock.o  ParseEntireEPOBlockRelax.o

