/*
 * ParseEntireEPOBlock
 * Date: Aug-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign goes here] gmail.com
 *
 */

#ifndef ParseEntireEPOBlock_h
#define ParseEntireEPOBlock_h

#include <iostream>
#include <string>
#include "gzstream/gzstream.h"
#include <vector> 
#include <stdlib.h>

#include "utils.h"

using namespace std;

//names prefixes for the species in EPO aligment 
static string speciesNameHuman="homo_sapiens";
static string speciesNameChimp="pan_troglodytes";
static string speciesNameGorilla="gorilla_gorilla";
static string speciesNameOrang="pongo_abelii";
static string speciesNameAnc="ancestral_sequences";

typedef struct{
    string chrName;
    unsigned int chrPos;
    char alleleHuman;
    char alleleChimp;
    char alleleHumanChimpAnc;
    char alleleGorilla;
    char alleleGorillaHumanChimpAnc;
    char alleleOrang;
    char alleleOrangGorillaHumanChimpAnc;
    bool cpg;
    bool blockFailed;
} snpLine;

typedef struct {
    unsigned int start;
    unsigned int end;
    string       filename;
    unsigned int line; 
} alignmentBlock;



unsigned int parseEntireEPOBlock(string inFile,unsigned int lineSought,unsigned int chrPos,unsigned int chrStart,unsigned int chrEnd,string chrName);

/* void  printsnpLine(snpLine & toPrint); */
/* void  printEmptyLine(string chrName,unsigned int chrPos); */
/* void  processLinesEPO(vector<snpLine> &  previousLines,snpLine & toPrint); */





/**
   This function formats and prints the information in
   "toPrint" to the stdout
	   
	 
   @param[in]     toPrint  Struct with the information from the EPO alignment.

*/
inline void  printsnpLine(snpLine & toPrint){
    cout<<toPrint.chrName<<"\t"
	<<toPrint.chrPos<<"\t"
	<<toPrint.alleleHuman<<"\t"
	<<toPrint.alleleHumanChimpAnc<<"\t"
	<<toPrint.alleleChimp<<"\t"

	<<toPrint.alleleGorillaHumanChimpAnc<<"\t"
	<<toPrint.alleleGorilla<<"\t"

	<<toPrint.alleleOrangGorillaHumanChimpAnc<<"\t"
	<<toPrint.alleleOrang<<"\t"
	<<toPrint.cpg<<"\t"
	<<toPrint.blockFailed<<"\t"

	<<endl;

}

/**
   This function prints an empty line if the block was not satifactory
	   
	 
   @param[in]     chrName Chromosome in use
   @param[in]     chrPos Position on the chrName

*/
inline void  printEmptyLine(string chrName,
			    unsigned int chrPos){
    snpLine toPrint;			
    toPrint.chrName = chrName;
    toPrint.chrPos  = chrPos;
    toPrint.alleleHuman                     = 'N';
    toPrint.alleleChimp                     = 'N';
    toPrint.alleleHumanChimpAnc             = 'N';
    toPrint.alleleGorilla                   = 'N';
    toPrint.alleleGorillaHumanChimpAnc      = 'N';
    toPrint.alleleOrang                     = 'N';
    toPrint.alleleOrangGorillaHumanChimpAnc = 'N';			    
    toPrint.cpg=false;
    toPrint.blockFailed=true;

    printsnpLine(toPrint);
}


/**
   This function detects CpGs and prints an the information provided by the EPO
   alignment . It will print the line in the vector "previousLines" and push toPrint into 
   the vector "previousLines".
	   
	 
   @param[in]     previousLines Vector of the previous lines found, this will be printed
   @param[in]     toPrint Current information, will be printed during the next loop.

*/
inline void  processLinesEPO(vector<snpLine> &  previousLines,snpLine & toPrint){
    if(previousLines.back().chrName != toPrint.chrName){
	cerr << "Wrong chr in EPO"<<endl;
	exit(1);       
    }

    //checking for CpG 
    //check only ancestors,
    if( (previousLines.back().chrPos+1) == toPrint.chrPos){
	if( (previousLines.back().alleleHuman         == 'C' &&
	     toPrint.alleleHuman                      == 'G' ) 
	    ||
	    (previousLines.back().alleleHumanChimpAnc == 'C' &&
	     toPrint.alleleHumanChimpAnc              == 'G' ) 
	    ||	    
	    // (previousLines.back().alleleChimp         == 'C' &&
	    //  toPrint.alleleChimp                      == 'G' )
	    // ||	   	    
	    (previousLines.back().alleleGorillaHumanChimpAnc == 'C' &&
	     toPrint.alleleGorillaHumanChimpAnc              == 'G' )
	    ||
	    // (previousLines.back().alleleGorilla       == 'C' &&
	    //  toPrint.alleleGorilla                    == 'G' )
	    (previousLines.back().alleleOrangGorillaHumanChimpAnc == 'C' &&
	     toPrint.alleleOrangGorillaHumanChimpAnc              == 'G' )
	    ){
	    previousLines.back().cpg = true;
	    toPrint.cpg              = true;
	}															       
    }

    printsnpLine(previousLines.back());

    previousLines.clear();
    previousLines.push_back(toPrint);
}



#endif
