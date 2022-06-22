/* 
  CalculateEnergy.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 4/2004, Justin Bois 1/2007
*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "pfuncUtilsHeader.h" //contains functions and structures
#include "DNAExternals.h"
//***********************************************************

DBL_TYPE naEnergy( char *prefix, int seq[]) {
  // Give energy of DNA strand with all other params set to defaults
  return naEnergyFull( prefix, seq, DNA, 1, 37, 1.0, 0.0, 0);
}



DBL_TYPE naEnergyFull( char prefix[], int inputSeq[], int naType, 
		       int dangles, DBL_TYPE temperature,
		       DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {
  // Give energy with possibly symmetry set to 1 and gets secondary structure from prefix.fold

  return naEnergyFullWithSym( prefix, inputSeq, naType, dangles, temperature,
			      1, sodiumconc, magnesiumconc, uselongsalt);
}


DBL_TYPE naEnergyFullWithSym( char prefix[], int inputSeq[], int naType, 
			      int dangles, DBL_TYPE temperature, int possibleSymmetry,
			      DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {

  fold thefold;
  DBL_TYPE energy;
  char *foldFile;
  int *seq;
  int size;
  int nicks[ MAXSTRANDS];
  int seqlength, nStrands;

  seqlength = getSequenceLengthInt( inputSeq, &nStrands);
  seq = (int*) malloc( (seqlength+1)*sizeof( int) ); 
  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);

  TEMP_K = temperature + ZERO_C_IN_KELVIN;
  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;

  size = strlen(prefix) + 6;
  foldFile = (char*) malloc( size*sizeof( char));
  strcpy( foldFile, prefix); 
  strcat( foldFile, ".fold");

  //LoadEnergies();  Energies are loaded in GetEnergy

  LoadFold( &thefold, foldFile); // get input file
  thefold.seq = seq;

  energy = GetEnergy( &thefold); //Calculates Energy  

  energy += LOG_FUNC( (DBL_TYPE) checkSymmetry( thefold.pairs, seqlength, nicks, 
                                          possibleSymmetry, nStrands))*kB*TEMP_K;
  energy += (BIMOLECULAR + SALT_CORRECTION)*(nStrands-1);


  free( foldFile);
  free( thefold.pairs);
  free( thefold.pknots);
  free( thefold.fixedBases);
  free( thefold.isNicked);
  free( seq);

  return energy;
}


/* ******************* */
DBL_TYPE naEnergyPairsOrParens( int *thepairs, char *parens, int seq[]) {
  // Returns energy of structure for DNA with all other parameters set to defaults

  return naEnergyPairsOrParensFull( thepairs, parens, seq, DNA, 1, 37, 1.0, 0.0, 0);
}

/* ******************* */
DBL_TYPE naEnergyPairsOrParensFull( int *thepairs, char *parens, 
                                   int inputSeq[], int naType,
				    int dangles, DBL_TYPE temperature,
				    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, 
				    int uselongsalt) {

  // Give energy with possibly symmetry set to 1
 
  return naEnergyPairsOrParensFullWithSym( thepairs, parens,
                                         inputSeq, naType, dangles, 
					   temperature, 1, sodiumconc, magnesiumconc,
					   uselongsalt);
}

/* ******************* */
DBL_TYPE naEnergyPairsOrParensFullWithSym( int *thepairs, char *parens, 
                                          int inputSeq[], int naType,
                                          int dangles, DBL_TYPE temperature,
					   int possibleSymmetry,
					   DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, 
					   int uselongsalt) {
                                            
  fold thefold;
  DBL_TYPE energy;
  int nStrands;
  int *seq; //without the strand breaks
  int seqlength;
  int nicks[ MAXSTRANDS];
  int i;
  
  for( i = 0; i < MAXSTRANDS; i++) { //initialize nicks array
    nicks[i] = -1;
  }
  
  seqlength = getSequenceLengthInt( inputSeq, &nStrands);
  seq = (int*) malloc( (seqlength+1)*sizeof( int) ); 
  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);
  
  TEMP_K = temperature + ZERO_C_IN_KELVIN;
  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;
  
  MakeFold( &thefold, seqlength, seq, parens, thepairs);
  thefold.seq = seq;
  
  thefold.isNicked = (int *) calloc( seqlength, sizeof( int));
  i = 0;
  while( nicks[i] != -1) {
    thefold.isNicked[ nicks[i++]] = 1;
  }
  
  energy = GetEnergy( &thefold); //Calculates Energy  
  energy += LOG_FUNC( (DBL_TYPE) checkSymmetry( thefold.pairs, seqlength, nicks, 
                                           possibleSymmetry, nStrands))*kB*TEMP_K;
  energy += (BIMOLECULAR + SALT_CORRECTION)*(nStrands-1);
  
  free( thefold.pairs);
  free( thefold.pknots);
  free( thefold.isNicked);
  free( seq);
  
  return energy;
}



/* ***************************************************** */
void MakeFold( fold *thefold, int seqlength, int seq[], char *parens, int *thepairs) {
  
  
  int init, i; // loop indices for initializations
  
  int pairsFromParens[ MAXSEQLENGTH];
  
  
  thefold->seqlength = seqlength;
  thefold->seq = seq;
  
  pairsFromParens[0] = -5; 
  if( parens != NULL) 
    getStructureFromParens( parens, pairsFromParens,  
                           thefold->seqlength);
  
  thefold->pairs = 
    (int*) calloc( thefold->seqlength+1, sizeof(int));
  if( thefold->pairs == NULL) {
    printf("Unable to allocate fold file!\n");
    exit(1);
  }
  
  thefold->pknots = 
    (int*) calloc( thefold->seqlength+1, sizeof(int));
  if( thefold->pknots == NULL) {
    printf("Unable to allocate fold file!\n");
    exit(1);
  }
  
  //enumerate the actual pknot tracker
  thefold->actualPknots = 
    (int*) calloc( thefold->seqlength+1, sizeof(int));
  if( thefold->actualPknots == NULL) {
    printf("Unable to allocate fold file for new pknot stuff!\n");
    exit(1);
  }

  for( init = 0; init <= thefold->seqlength; init++) {
    thefold->pairs[init] = -1;
    thefold->pknots[init] = -1;
  }
  

  //this populates the pairs lisin fold
  for( init = 0; init <= thefold->seqlength - 1; init++) {
    if( parens != NULL) 
      thefold->pairs[init] = pairsFromParens[init];
    else
      thefold->pairs[init] = thepairs[init];
  }
  
  
  //new pknot routine

  bool isPair = FALSE;
z
  //initialize to TRUE and first nuc will actually change it
  bool noPair_NextNuc = FALSE;
  
  //nucs are 0 indexed so start nuc is index 0
  int searchNuc_index = 0;

  //search nuc comlempentary pair index
  int searchNuc_compPair = -1;

  //range [0] = i-1, [1] = start of this run. This is in 1 based index
  //SL = Smal to large nuc order, LS = Large to small nuc order
  //Small to large
  int nucsLenght= thefold->seqlength;
  int *nucs_Behind_SL = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Behind_SL_Count =0;
  int *nucs_Front_SL = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Front_SL_Count =nucsLenght;

  //large to small
  int *nucs_Behind_LS = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Behind_LS_Count =0;
  int *nucs_Front_LS = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Front_LS_Count =nucsLenght;

  int *gapNucs = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int gapNucs_Count =0;

  int *pairsTracker = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int pairsTracker_Count =0;

  int indexToAdd = -1;


  int nucsLenght= thefold->seqlength;
  int *nucs_Behind_SL_Test = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Behind_SL_Test_Count =0;
  int *nucs_Front_Test_SL = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Front_Test_SL_Count =nucsLenght;

  //large to small
  int *nucs_Behind_Test_LS = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Behind_Test_LS_Count =0;
  int *nucs_Front_Test_LS = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int nucs_Front_Test_LS_Count =nucsLenght;

  int *gapNucs_Test = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int gapNucs_Test_Count =0;

  int *pairsTracker_Test = (int*) calloc( thefold->seqlength+1, sizeof(int));
  int pairsTracker_Test_Count =0;

  int indexToAdd = -1;
  //initialize nucs front list
  for (int index = 1; index <= nucsLenght, nucsLenght++)
  {
    nucs_Front_SL[index] = index;
    
    nucs_Front_LS[index] = index;
    nucs_Behind_SL[index] = -1;    
    nucs_Behind_LS[index] = -1;
    gapNucs[index] = -1;
    pairsTracker[index] = -1;
  }

  bool inGap = FALSE;
  bool inStack = FALSE;
  bool inStack_nextNuc = FALSE;
  bool isPknot_suspected = FALSE;
  bool isPknot_confident = FALSE;


  while (noPair_NextNuc == TRUE)
  { 
   
    //this is the entry for each nuc and stuff happens
    //now act on searchNuc_index nucleotide

    //get the  actual nuc number and th complementary pair
    int actualNuc = searchNuc_index+1;
    searchNuc_compPair = thefold->pairs[searchNuc_index];
    int actualNuc_compPair = searchNuc_compPair+1;

    //add the last nuc to the list
    //add to the list immediatly when it is hit
    indexToAdd = nucsLenght-nucs_Behind_SL_Count;
    nucs_Behind_LS[indexToAdd]=actualNuc;
    nucs_Behind_SL_Count++;

    if (inGap==TRUE)
    {
      //indexToAdd = nucsLenght-gapNucs_Count;
      //gapNucs[indexToAdd]=actualNuc;
    }

    if (inStack==TRUE)
    {

    }

    //now remove from the front
    for (int index = 1; index <= nucsLenght, nucsLenght++)
    {
      if (nucs_Front_SL[index] == actualNuc)
      {
        nucs_Front_SL[index] = -1;
        nucs_Front_SL_Count--;
      }
    }
    
    //now we know what should be behind and ahead and what has been found in its gap if applicable

    //if not -1 then it is paired if -1 then gap
    if (searchNuc_compPair != -1)
    {
      isPair=TRUE;
    }
    
    if (isPair==TRUE)
    {
      //this nuc is paired so record and poke around poke around
      indexToAdd = nucsLenght-pairsTracker_Count;
      pairsTracker[indexToAdd]=searchNuc_index;
      pairsTracker_Count++;

      indexToAdd = nucsLenght-pairsTracker_Test_Count;
      pairsTracker_Test[indexToAdd]=searchNuc_index;
      pairsTracker_Test_Count++;

      //now need to test to see if next makes sense based on trackers
      bool isValid_inFront = FALSE;
      for (int index = 1; index <=nucsLenght; index++)
      {
        if (nucs_Front_SL[index]==nextSearchNuc_compPair)
        {
          //its in the listof expected front nucs so its valid ad should be able to loop to next nuck
          isValid_inFront = TRUE;     
        }
      }


      bool testNextPair =FALSE; 
      if ( isValid_inFront==TRUE)
      {
        //we are progressing through a normal sequence and stacks
        //check next nuc pair
        testNextPair=TRUE;
      }
      else
      {
        //this is potentially in a pknot 
        isPknot_suspected=TRUE;
      }
      
      //now check if next nuc after searchNuc_compPair is a pair
      int nextSearchNuc=searchNuc_index+1;
      int nextSearchNuc_compPair = thefold->pairs[nextSearchNuc];

      if (testNextPair==TRUE)
      {
        if (nextSearchNuc_compPair == -1)
        {
          //its in a gap so record that
          inGap=TRUE;
          indexToAdd = nucsLenght-gapNucs_Count;
          gapNucs[indexToAdd]=actualNuc;
          gapNucs_Count++;

          //should be in a loop now so check nucs
          //fold should be in expected range after last pairing

          //check if nucpair_compPair isin front
          isValid_inFront = FALSE;
          for (int index = 1; index <=nucsLenght; index++)
          {
            if (nucs_Front_SL[index]==nextSearchNuc_compPair)
            {
              //its in the listof expected front nucs so its valid ad should be able to loop to next nuck
              isValid_inFront = TRUE;

              //this basicaly would be the equivialant of coming upon a pair while walking a gap and then jumping the pair
              //the next nucs is a gap and is part of expected nucs so its basically in a lcoal loop            
            }
          }

          if(isValid_inFront==TRUE)
          {
            //this basicaly would be the equivialant of coming upon a pair while walking a gap and then jumping the pair
            //the next nucs is a gap and is part of expected nucs so its basically in a lcoal loop   

            //the nucpair is valid as a standard pair as its in the front
            //set test front array to exclude everything before the next nucsearch_compPair
            for(int index = 1; index <=nextSearchNuc_compPair, index++)
            {
              nucs_Front_Test_SL[index]=-1;
            } 

          }
          else
          {
            //its not in the expected front so it is probably a pknot but need to test
            isPknot_suspected=TRUE;

            //this basicaly would be the equivialant of coming upon a pair while walking a gap and then jumping the pair
            //the next nucs is a pknot
          }
        }
        else
        {
          //there is a pair net so figure it out
          //this is the next nuc pair we are looking at
          //should be in a stack so walk it and make sure it makes sense
          //basically repeat while loop now
          inStack=TRUE;
          indexToAdd = nucsLenght-pairsTracker_Test_Count;
          pairsTracker_Test[indexToAdd]=nextSearchNuc_compPair;
          pairsTracker_Test_Count++;
          
        }

      }
      else
      {
        //no point in testing next pair but need a good reason why
      }
     
    }
  }

  //at this point the while loop


  //the following pknot finding routine should be optimized
  for( init = 0; init <= thefold->seqlength-1; init++) {
    if( thefold->pairs[init] > init) {
      for( i = 0; i < init; i++) {
        if( thefold->pairs[i] > init && 
        thefold->pairs[i] < thefold->pairs[init]) {
          if( thefold->pknots[i] == -1) {
            
            thefold->pknots[i] = thefold->pairs[init];
            thefold->pknots[ thefold->pairs[init]] = i;
          }
          break;
        }
      }
      
    }
  }
  
}




/* *************** */


