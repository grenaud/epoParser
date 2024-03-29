/*
 * ParseEntireEPOBlockRelax
 * Date: Aug-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign goes here] gmail.com
 *
 */

#ifndef ParseEntireEPOBlockRelax_h
#define ParseEntireEPOBlockRelax_h

#include <iostream>
#include <string>
#include "gzstream/gzstream.h"
#include <vector> 
#include <stdlib.h>

#include "libgab.h"
/* #include "ParseEntireEPOBlock.h" */

using namespace std;

/* typedef struct{ */
/*     string chrName; */
/*     unsigned int chrPos; */
/*     char alleleHuman; */
/*     char alleleChimp; */
/*     char alleleHumanChimpAnc; */
/*     char alleleGorilla; */
/*     char alleleGorillaHumanChimpAnc; */
/*     char alleleOrang; */
/*     char alleleOrangGorillaHumanChimpAnc; */
/*     bool cpg; */
/*     bool blockFailed; */
/* } snpLine; */

/* typedef struct { */
/*     unsigned int start; */
/*     unsigned int end; */
/*     string       filename; */
/*     unsigned int line;  */
/* } alignmentBlock; */


pair<unsigned int,bool> parseEntireEPOBlockRelax(string inFile,unsigned int lineSought,unsigned int chrPos,unsigned int chrStart,unsigned int chrEnd,string chrName);

/* void  printsnpLine(snpLine & toPrint); */
/* void  printEmptyLine(string chrName,unsigned int chrPos); */
/* void  processLinesEPO(vector<snpLine> &  previousLines,snpLine & toPrint); */

#endif
