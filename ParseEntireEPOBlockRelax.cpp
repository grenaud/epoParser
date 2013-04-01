/*
 * ParseEntireEPOBlockRelax
 * Date: Aug-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign goes here] gmail.com
 *
 */

#include "ParseEntireEPOBlockRelax.h"
#include "ParseEntireEPOBlock.h"



//names prefixes for the species in EPO aligment 
// string speciesNameHuman="homo_sapiens";
// string speciesNameChimp="pan_troglodytes";
// string speciesNameGorilla="gorilla_gorilla";
// string speciesNameOrang="pongo_abelii";
// string speciesNameAnc="ancestral_sequences";

// /**
//    This function formats and prints the information in
//    "toPrint" to the stdout
	   
	 
//    @param[in]     toPrint  Struct with the information from the EPO alignment.

// */
// inline void  printsnpLine(snpLine & toPrint){
//     cout<<toPrint.chrName<<"\t"
// 	<<toPrint.chrPos<<"\t"
// 	<<toPrint.alleleHuman<<"\t"
// 	<<toPrint.alleleHumanChimpAnc<<"\t"
// 	<<toPrint.alleleChimp<<"\t"

// 	<<toPrint.alleleGorillaHumanChimpAnc<<"\t"
// 	<<toPrint.alleleGorilla<<"\t"

// 	<<toPrint.alleleOrangGorillaHumanChimpAnc<<"\t"
// 	<<toPrint.alleleOrang<<"\t"
// 	<<toPrint.cpg<<"\t"
// 	<<toPrint.blockFailed<<"\t"

// 	<<endl;

// }

// /**
//    This function prints an empty line if the block was not satifactory
	   
	 
//    @param[in]     chrName Chromosome in use
//    @param[in]     chrPos Position on the chrName

// */
// void  printEmptyLine(string chrName,
// 			    unsigned int chrPos){
//     snpLine toPrint;			
//     toPrint.chrName = chrName;
//     toPrint.chrPos  = chrPos;
//     toPrint.alleleHuman                     = 'N';
//     toPrint.alleleChimp                     = 'N';
//     toPrint.alleleHumanChimpAnc             = 'N';
//     toPrint.alleleGorilla                   = 'N';
//     toPrint.alleleGorillaHumanChimpAnc      = 'N';
//     toPrint.alleleOrang                     = 'N';
//     toPrint.alleleOrangGorillaHumanChimpAnc = 'N';			    
//     toPrint.cpg=false;
//     toPrint.blockFailed=true;

//     printsnpLine(toPrint);
// }


// /**
//    This function detects CpGs and prints an the information provided by the EPO
//    alignment . It will print the line in the vector "previousLines" and push toPrint into 
//    the vector "previousLines".
	   
	 
//    @param[in]     previousLines Vector of the previous lines found, this will be printed
//    @param[in]     toPrint Current information, will be printed during the next loop.

// */
// inline void  processLinesEPO(vector<snpLine> &  previousLines,snpLine & toPrint){
//     if(previousLines.back().chrName != toPrint.chrName){
// 	cerr << "Wrong chr in EPO"<<endl;
// 	exit(1);       
//     }

//     //checking for CpG 
//     //check only ancestors,
//     if( (previousLines.back().chrPos+1) == toPrint.chrPos){
// 	if( (previousLines.back().alleleHuman         == 'C' &&
// 	     toPrint.alleleHuman                      == 'G' ) 
// 	    ||
// 	    (previousLines.back().alleleHumanChimpAnc == 'C' &&
// 	     toPrint.alleleHumanChimpAnc              == 'G' ) 
// 	    ||	    
// 	    // (previousLines.back().alleleChimp         == 'C' &&
// 	    //  toPrint.alleleChimp                      == 'G' )
// 	    // ||	   	    
// 	    (previousLines.back().alleleGorillaHumanChimpAnc == 'C' &&
// 	     toPrint.alleleGorillaHumanChimpAnc              == 'G' )
// 	    ||
// 	    // (previousLines.back().alleleGorilla       == 'C' &&
// 	    //  toPrint.alleleGorilla                    == 'G' )
// 	    (previousLines.back().alleleOrangGorillaHumanChimpAnc == 'C' &&
// 	     toPrint.alleleOrangGorillaHumanChimpAnc              == 'G' )
// 	    ){
// 	    previousLines.back().cpg = true;
// 	    toPrint.cpg              = true;
// 	}															       
//     }

//     printsnpLine(previousLines.back());

//     previousLines.clear();
//     previousLines.push_back(toPrint);
// }



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
pair<unsigned int,bool>  parseEntireEPOBlockRelax(string inFile,unsigned int lineSought,unsigned int chrPos,unsigned int chrStart,unsigned int chrEnd,string chrName){
     // ;//= "/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v64/Compara.6_primates_EPO.chrX_1.emf.gz";
     // ; //=1000022;
     // ; //=5039413;

    unsigned int lineCounter=0;
    igzstream igzStream;
    bool header=true;
    bool data  =false;
    string line;

    bool validBlock=true;

    cerr<<"Relaxed version"<<endl;

    vector<string> headerLines;
    vector<snpLine> previousLines;
    int indexHuman=-1;
    vector<int> indexChimps;

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
			// if( indexHuman!=-1 && starsWith(headerLines[i],speciesNameHuman)){//we do not allow two humans
			//     validBlock=false;
			// }

			//if( indexHuman==-1 && starsWith(headerLines[i],speciesNameHuman) ){
			if( starsWith(headerLines[i],speciesNameHuman) ){

			    vector<string> fields = allTokens(headerLines[i],' ');
			    

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
			
			if( starsWith(headerLines[i],speciesNameChimp)){
			    indexChimps.push_back(i);			   
			}

			// //Detect chimp
			// if( indexChimp!=-1 && starsWith(headerLines[i],speciesNameChimp)){			
			//     validBlock=false;
			// }
			// if( indexChimp==-1 && starsWith(headerLines[i],speciesNameChimp)){			
			//     if( starsWith(headerLines[i-1],speciesNameAnc) ){
			// 	indexChimp=i;
			// 	indexHumanChimpAnc=i-1;
			//     }else{
			// 	validBlock=false;
			//     }
			// }

			// //Detect Gorilla
			// if(  indexGorilla!=-1 && starsWith(headerLines[i],speciesNameGorilla)){			
			//     validBlock=false;
			// }
			// if(  indexGorilla==-1 && starsWith(headerLines[i],speciesNameGorilla)){			
			//     if( starsWith(headerLines[i-1],speciesNameAnc) ){
			// 	indexGorilla=i;
			// 	indexGorillaHumanChimpAnc=i-1;				
			//     }else{
			// 	validBlock=false;
			//     }
			// }
		
			// //Detect orang
			// if( indexOrang!=-1 && starsWith(headerLines[i],speciesNameOrang) ){
			//     validBlock=false;
			// }
			// if( indexOrang==-1 && starsWith(headerLines[i],speciesNameOrang) ){
			//     if( starsWith(headerLines[i-1],speciesNameAnc)){
			// 	indexOrang=i;
			// 	indexOrangGorillaHumanChimpAnc=i-1;
			//     }else{
			// 	validBlock=false;
			//     }
			// }
			
			

		    } 

		}//end header lines

		//reached the end of the species
		//reached the end of the species
		if(line.substr(0,4) == "TREE"){
		    // if(indexHuman != 0 ){ //the human should always be the first block
		    // 	validBlock=false;
		    // }

		    // cout<<line<<endl;
		    if(indexChimps.empty()){
			validBlock=false;
		    }

		    cerr<<"index   human  "<<indexHuman<<endl;
		    cerr<<"indices chimp  "<<vectorToString(indexChimps)<<endl;
		    
		    // cerr<<"index chimp  "<<indexChimp<<endl;
		    // cerr<<"index CH_anc "<<indexHumanChimpAnc<<endl;
		    // cerr<<"index grila  "<<indexGorilla<<endl;
		    // cerr<<"index gr_anc "<<indexGorillaHumanChimpAnc<<endl;
		    // cerr<<"index orang  "<<indexOrang<<endl;
		    // cerr<<"index or_anc "<<indexOrangGorillaHumanChimpAnc<<endl;

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

			// if(indexChimp != -1){
			//     toPrint.alleleChimp         = upper(line[indexChimp]);
			//     //toPrint.alleleHumanChimpAnc = upper(line[indexHumanChimpAnc]);
			// }else{
			toPrint.alleleChimp         = 'N';
			toPrint.alleleHumanChimpAnc = 'N';		       
			// }
			//toPrint.alleleChimp         = 'N';
			for(unsigned int indexChimp=0;indexChimp<indexChimps.size();indexChimp++){
			    if( indexChimp == 0){
				toPrint.alleleChimp = upper(line[indexChimp]);
			    }else{
				if(toPrint.alleleChimp == 'N'){
				    //leave it
				}else{
				    if( toPrint.alleleChimp != upper(line[indexChimp]) ){ //if the chimp alleles differ, we assign an 'N'
					toPrint.alleleChimp         = 'N';
				    }
				}
			    }
			}
			
			// if(indexGorillaHumanChimpAnc != -1){
			//     toPrint.alleleGorilla              = upper(line[indexGorilla]);
			//     toPrint.alleleGorillaHumanChimpAnc = upper(line[indexGorillaHumanChimpAnc]);
			// }else{
			toPrint.alleleGorilla                   = 'N';
			toPrint.alleleGorillaHumanChimpAnc      = 'N';
			// }

			// if(indexOrangGorillaHumanChimpAnc != -1){
			//     toPrint.alleleOrangGorillaHumanChimpAnc = upper(line[indexOrangGorillaHumanChimpAnc]);
			//     toPrint.alleleOrang                     = upper(line[indexOrang]);
			// }else{
			toPrint.alleleOrangGorillaHumanChimpAnc = 'N';
			toPrint.alleleOrang                     = 'N';
			// }

			if(toPrint.alleleHuman  != '-' ){ //if there is an insert in human, skip it, no point in printing those
			    toPrint.cpg=false;
			    if(previousLines.empty()){
				previousLines.push_back(toPrint);
			    }else{ //not empty				
				processLinesEPO(previousLines,toPrint);
			    }
			}

			if(isValidDNA(upper(line[indexHuman]))){
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

    pair<unsigned int,bool> toReturn  (chrPos,validBlock);
    
    return toReturn;
    //    return chrPos;
}
