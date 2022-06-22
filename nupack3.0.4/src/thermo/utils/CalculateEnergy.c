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
  int *nucs_Behind_SL = (int*) calloc( thefold->seqlength, sizeof(int));
  int nucs_Behind_SL_Count =0;
  int *nucs_Front_SL = (int*) calloc( thefold->seqlength, sizeof(int));
  int nucs_Front_SL_Count =0;

  //large to small
  int *nucs_Behind_LS = (int*) calloc( thefold->seqlength, sizeof(int));
  int nucs_Behind_LS_Count =0;
  int *nucs_Front_LS = (int*) calloc( thefold->seqlength, sizeof(int));
  int nucs_Front_LS_Count =0;


  //initialize nucs front list

  while (noPair_NextNuc == TRUE)
  { 
    //this is the entry for each nuc and stuff happens
    //now act on searchNuc_index nucleotide

    //add the last nuc to the list
    //add to the list immediatly when it is hit
    int indexToAdd = nucsLenght-nucs_Behind_SL_Count;
    nucs_Behind_LS[indexToAdd]=searchNuc_index;
    

    //get the complementary pair
    searchNuc_compPair = thefold->pairs[searchNuc_index]


    //if not -1 then it is paired if -1 then gap
    if (searchNuc_compPair != -1)
    {
      isPair=TRUE;
    }
    
    if (isPair==TRUE)
    {
      //check if next nuc is a pair
      int nextSearchNuc=searchNuc_index+1;
      int nextSearchNuc_compPair = thefold->pairs[searchNuc_index];

      if (nextSearchNuc_compPair == -1)
      {
        //its in a gap so record that


      }
      else
      {
        //there is a pair so keep going and loging
      }


    }

  }



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


