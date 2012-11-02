CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS = -lm -O3 -lz -I${LIBGAB}  -c
LDFLAGS  = -lz


all: epoParser ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o



ParseEntireEPOBlock.o:	${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlock.cpp

ParseEntireEPOBlockRelax.o:	${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} ParseEntireEPOBlockRelax.cpp

epoParser.o:	epoParser.cpp
	${CXX} ${CXXFLAGS} epoParser.cpp


epoParser:	epoParser.o ${LIBGAB}utils.o  gzstream/gzstream.o ParseEntireEPOBlock.o ParseEntireEPOBlockRelax.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 



clean :
	rm -f epoParser.o epoParser ParseEntireEPOBlock.o  ParseEntireEPOBlockRelax.o

