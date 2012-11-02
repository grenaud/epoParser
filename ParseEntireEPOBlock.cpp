/*
 * ParseEntireEPOBlock
 * Date: Aug-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign goes here] gmail.com
 *
 */

#include "ParseEntireEPOBlock.h"




/**
   This function parses a given EPO block in a file and prints the contents 
   on the screen (cout).
	   
	 
   @param[in]     inFile  full path to the gzipped file containing the block
   @param[in]     lineSought at which line does the first SEQ start ?
   @param[in]     chrPos Chromosomal position to print
   @param[in]     chrStart Chromosomal position where the alignment begins (on human)
   @param[in]     chrEnd At the end, chrPos should be equal to chrPos
   @param[in]     chrName Name of chr for this alignment (again on human)

   @return Last chromosome position printed on the screen
*/

unsigned int parseEntireEPOBlock(string inFile,unsigned int lineSought,unsigned int chrPos,unsigned int chrStart,unsigned int chrEnd,string chrName){
     // ;//= "/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v64/Compara.6_primates_EPO.chrX_1.emf.gz";
     // ; //=1000022;
     // ; //=5039413;

    unsigned int lineCounter=0;
    igzstream igzStream;
    bool header=true;
    bool data  =false;
    string line;

    bool validBlock=true;

    vector<string> headerLines;
    vector<snpLine> previousLines;
    int indexHuman=-1;
    int indexChimp=-1;
    int indexHumanChimpAnc=-1;
    int indexGorilla=-1;
    int indexGorillaHumanChimpAnc=-1;
    int indexOrang=-1;
    int indexOrangGorillaHumanChimpAnc=-1;

    igzStream.open(inFile.c_str());
    if( igzStream.fail()){
	cerr << "Cannot open file "<<inFile<<endl;
	exit(1);      
    }
    while(getline(igzStream, line)){


	if(lineCounter>=lineSought){
	   
	    if(header){

		if(line.substr(0,3) == "SEQ"){
		    headerLines.push_back( line.substr(4) );
		}else{
		    header=false;

		    //Here we detect the indices of the species in the header
		    for(int i=0;i<(headerLines.size());i++){
			cerr<<i<<"\t"<<headerLines[i]<<"\t"<<validBlock<<endl;

			//detect human
			if( indexHuman!=-1 && starsWith(headerLines[i],speciesNameHuman)){//we do not allow two humans
			    validBlock=false;
			}

			if( indexHuman==-1 && starsWith(headerLines[i],speciesNameHuman) ){
			    vector<string> fields = allTokens(headerLines[i]," ");
			    

			    if(             fields[1] == chrName  &&
					    (string2uint(fields[2]) == chrStart) &&
					    (string2uint(fields[3]) >= chrEnd) ){
				if(fields[4] != "1"){
				    cerr<<"Wrong strand"<<endl;
				    validBlock=false;
				    //exit(1);
				}
				indexHuman=i;
			    }
			}			
			

			//Detect chimp 		       
			if( indexChimp!=-1 && starsWith(headerLines[i],speciesNameChimp)){			
			    validBlock=false;
			}
			if( indexChimp==-1 && starsWith(headerLines[i],speciesNameChimp)){			
			    if( starsWith(headerLines[i-1],speciesNameAnc) ){
				indexChimp=i;
				indexHumanChimpAnc=i-1;
			    }else{
				validBlock=false;
			    }
			}

			//Detect Gorilla
			if(  indexGorilla!=-1 && starsWith(headerLines[i],speciesNameGorilla)){			
			    validBlock=false;
			}
			if(  indexGorilla==-1 && starsWith(headerLines[i],speciesNameGorilla)){			
			    if( starsWith(headerLines[i-1],speciesNameAnc) ){
				indexGorilla=i;
				indexGorillaHumanChimpAnc=i-1;				
			    }else{
				validBlock=false;
			    }
			}
		
			//Detect orang
			if( indexOrang!=-1 && starsWith(headerLines[i],speciesNameOrang) ){
			    validBlock=false;
			}
			if( indexOrang==-1 && starsWith(headerLines[i],speciesNameOrang) ){
			    if( starsWith(headerLines[i-1],speciesNameAnc)){
				indexOrang=i;
				indexOrangGorillaHumanChimpAnc=i-1;
			    }else{
				validBlock=false;
			    }
			}
			

		    } 

		}//end header lines

		//reached the end of the species
		if(line.substr(0,4) == "TREE"){
		    if(indexHuman != 0 ){ //the human should always be the first block
			validBlock=false;
		    }

		    // cout<<line<<endl;
		    
		    cerr<<"index human  "<<indexHuman<<endl;
		    cerr<<"index chimp  "<<indexChimp<<endl;
		    cerr<<"index CH_anc "<<indexHumanChimpAnc<<endl;
		    cerr<<"index grila  "<<indexGorilla<<endl;
		    cerr<<"index gr_anc "<<indexGorillaHumanChimpAnc<<endl;
		    cerr<<"index orang  "<<indexOrang<<endl;
		    cerr<<"index or_anc "<<indexOrangGorillaHumanChimpAnc<<endl;

		    if(validBlock){
			cerr<<"Block accepted"<<endl;
		    }else{
			cerr<<"Block rejected"<<endl;
		    }
		}
	
	    }else{//end if header



		if(data){
		    if(line == "//"){ //reached the end of the data
			if(!previousLines.empty()){
			    printsnpLine(previousLines.back());
			}
			chrPos-=1;
			break; //exit main loop
		    }

		    if(validBlock){
			snpLine toPrint;
			
			toPrint.chrName = chrName;
			toPrint.chrPos  = chrPos;
			toPrint.alleleHuman         = upper(line[indexHuman]);
			toPrint.blockFailed=!validBlock;

			if(indexChimp != -1){
			    toPrint.alleleChimp         = upper(line[indexChimp]);
			    toPrint.alleleHumanChimpAnc = upper(line[indexHumanChimpAnc]);
			}else{
			    toPrint.alleleChimp         = 'N';
			    toPrint.alleleHumanChimpAnc = 'N';		       
			}

			if(indexGorillaHumanChimpAnc != -1){
			    toPrint.alleleGorilla              = upper(line[indexGorilla]);
			    toPrint.alleleGorillaHumanChimpAnc = upper(line[indexGorillaHumanChimpAnc]);
			}else{
			    toPrint.alleleGorilla                   = 'N';
			    toPrint.alleleGorillaHumanChimpAnc      = 'N';
			}

			if(indexOrangGorillaHumanChimpAnc != -1){
			    toPrint.alleleOrangGorillaHumanChimpAnc = upper(line[indexOrangGorillaHumanChimpAnc]);
			    toPrint.alleleOrang                     = upper(line[indexOrang]);
			}else{
			    toPrint.alleleOrangGorillaHumanChimpAnc = 'N';
			    toPrint.alleleOrang                     = 'N';
			}

			// if( (toPrint.alleleChimp                     == '-' && //those, they have inserts in the other species
			//      toPrint.alleleHumanChimpAnc             == '-' &&
			//      toPrint.alleleGorilla                   == '-' &&
			//      toPrint.alleleGorillaHumanChimpAnc      == '-' &&
			//      toPrint.alleleOrang                     == '-' &&
			//      toPrint.alleleOrangGorillaHumanChimpAnc == '-' )){
			//     //skip 
			// }else{

			if(toPrint.alleleHuman  != '-' ){ //if there is an insert in human, skip it, no point in printing those
			    toPrint.cpg=false;
			    if(previousLines.empty()){
				previousLines.push_back(toPrint);
			    }else{ //not empty				
				processLinesEPO(previousLines,toPrint);
			    }
			}
			    
			    //			}
			
			if(isValidDNA(upper(line[indexHuman]))){ //we only increase the counter for A,C,G,T,N
			    chrPos++;
			}
			    
		    }else{ //not valid block
			if(isValidDNA(upper(line[indexHuman]))){ //even if the strand is wrong, the number of bases should be the same
			    printEmptyLine(chrName,chrPos);
			    chrPos++;
			}
		    }
		    
		}else{ //haven't reached data yet

		    if(line.substr(0,4) == "DATA"){
			data=true;
		    }
		}

	    } //end not header

	    
	} // not reached line sought
	lineCounter++;
    } //end of line
    igzStream.close();

    if(chrPos!= chrEnd){
	cerr<<"The end was not equal to the iterator over the chromsome length"<<endl;
	exit(1);
    }
    return chrPos;
}
