/*
 * epoParser
 * Date: Aug-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign goes here] gmail.com
 *
 */



#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string> 
#include <vector> 
#include <map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "gzstream/gzstream.h"
#include "ParseEntireEPOBlock.h"

using namespace std;



bool cmdAlignBlocks (alignmentBlock i,alignmentBlock j) { return (i.start<j.start); }




int main (int argc, char *argv[]) {

    string usage=string(""+string(argv[0])+" options [hsa_emf.index] [fasta index] [name chr]"+
			"\nThis program an EPO alignment index and fasta index\n"+
			"for the genome and produces a tab delimited output\n");
			
    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    string epoIndex      =argv[argc-3];

    string epoIndexPath;
    size_t found=epoIndex.find_last_of("/");

    if(found == string::npos){
	cerr << "Must specify full path for hsa_emf.index"<<endl;
	return 1;       
    }else{
	epoIndexPath=epoIndex.substr(0,found);
    }

    string fastaIndex=argv[argc-2];

    string chrName   =argv[argc-1];

    ifstream myFile;
    string line;

    map<string,unsigned int> name2Length;
    myFile.open(fastaIndex.c_str(), ios::in);
    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    vector<string> fields = allTokens(line,"\t");
	    //cout<<fields[0]<<endl;
	    name2Length[fields[0]]=string2uint(fields[1]);
	}
    }
    myFile.close();

    if(name2Length.find(chrName) == name2Length.end()){ //chr not found
	cerr << "Chromosome name not found in fasta index file "<<usage<<endl;
	return 1;     
    }


    vector<alignmentBlock> alignmentBlocks;

    myFile.open(epoIndex.c_str(), ios::in);
    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    vector<string> fields = allTokens(line,"\t");
	    if( fields[0] == chrName ){
		alignmentBlock toadd;
		toadd.start    = string2uint(fields[1])+1;
		toadd.end      = string2uint(fields[2]);
		toadd.filename = fields[3];
		toadd.line     = string2uint(fields[5]);
		if(toadd.end<toadd.start){
		    cerr << "Start greater than end in line  "<<line<<endl;
		    return 1;     
		}		   
		alignmentBlocks.push_back(toadd);
	    }
	}
    }
    myFile.close();

    sort (alignmentBlocks.begin(), alignmentBlocks.end(), cmdAlignBlocks); 
    for(int i=0;i<alignmentBlocks.size();i++){
    	cerr<<i<<alignmentBlocks[i].start<<"\t"
    	    <<alignmentBlocks[i].end<<"\t"
    	    <<alignmentBlocks[i].filename<<"\t"
    	    <<alignmentBlocks[i].line<<"\n";
    }

    for(int i=0;i<(alignmentBlocks.size()-1);i++){
	if(alignmentBlocks[i+1].start<alignmentBlocks[i].end){
	    cerr << "Start ("<<alignmentBlocks[i+1].start<<") greater than end ("<<alignmentBlocks[i].end<<") of next  "<<endl;	   
	    return 1;       
	}       	   
    }

    //MAIN LOOP
    int indexAlignVector=0;
    bool noMoreBlocks=false;
    for(unsigned int chrIndex=1;
	chrIndex<name2Length[chrName];
	chrIndex++){

	if(!noMoreBlocks && chrIndex>=alignmentBlocks[indexAlignVector].start){

	    cerr<<endl;
	    cerr<<"block "<<indexAlignVector<<" of "<<alignmentBlocks.size()<<endl;
	    cerr<<"filename  "<<alignmentBlocks[indexAlignVector].filename<<endl;
	    cerr<<"lineBlock "<<alignmentBlocks[indexAlignVector].line<<endl;

	    cerr<<"start  "<<alignmentBlocks[indexAlignVector].start<<endl;
	    cerr<<"end    "<<alignmentBlocks[indexAlignVector].end<<endl;

	    cerr<<"chrindex  "<<chrIndex<<endl;
	    cerr<<"chrname   "<<chrName<<endl;
	    cerr<<"indexAlignVector   "<<indexAlignVector<<endl;

	    chrIndex=parseEntireEPOBlock(epoIndexPath+"/"+alignmentBlocks[indexAlignVector].filename,
					 alignmentBlocks[indexAlignVector].line,					 
					 chrIndex,
					 alignmentBlocks[indexAlignVector].start,
					 alignmentBlocks[indexAlignVector].end,
					 chrName);
	    //indexAlignVector=min(indexAlignVector+1,int(alignmentBlocks.size()-1));
	    if(chrIndex != alignmentBlocks[indexAlignVector].end){
		cerr<<"The EPO does not contain the right number of nucleotides"<<endl;
		return 1;
	    }
	    indexAlignVector++;
	    if(indexAlignVector == int(alignmentBlocks.size())){ //encountered last block
		noMoreBlocks=true;
	    }
	    
	}else{
	    printEmptyLine(chrName,chrIndex);				
	}
    }


    
    // parseEntireEPOBlock("/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v64/Compara.6_primates_EPO.chrX_1.emf.gz",
    // 			1000022, //line to begin
    // 			5039413, //
    // 			chrName);



    cerr<<"Terminated gracefully"<<endl;
    
    return 0;
}

