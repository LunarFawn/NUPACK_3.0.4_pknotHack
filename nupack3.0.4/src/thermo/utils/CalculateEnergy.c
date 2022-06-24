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


typedef struct PknotDetectionData
{ 
  //main seeking trackers
  //small to large trackers
  int* small_behind_y_trackerList;
  int small_behind_y_trackerList_Count;
  int* small_front_y_trackerList;
  int small_front_y_trackerList_Count;
  
  //large to small trackers
  int* large_behind_y_trackerList;
  int large_behind_y_trackerList_Count;
  int* large_front_y_trackerList;
  int large_front_y_trackerList_Count;
  
  
  //gap and trackers 
  int* gapNucs_trackerList;
  int gapNucs_trackerList_Count;
  int* pairsNucs_trackerList;
  int pairsNucs_trackerList_Count;
  
  //utils
  int nucsLenght;
  int fullSequenceLenght;
  int thisSegmentLenght;
  int currentNuc_y;
  int currentNuc_d;
  int endNuc_y;
  bool isPaired_current_y;
  bool is_yGRTd_current;
  bool is_dGRTy_current;
  
  int nextNuc_y;
  int nextNuc_d;
  bool isPaired_next_y;
  bool is_yGRTd_next;
  bool is_dGRTy_next;
  
  //loop, bulge, stack, pknot boolean logic variables
  bool inGap;

  bool isLoop_suspected;
  bool isLoop_confident;
  bool isLoop_confirmed;

  bool isBulge_suspected;
  bool isBulge_confident;
  bool isBulge_confirmed; 
 
  bool isStack_suspected;
  bool isStack_confident;
  bool isStack_confirmed;

  bool isPknot_suspected;
  bool isPknot_confident;
  bool isPknot_confirmed;
};


void InitalizePknotStruct(fold *thefold, struct PknotDetectionData *tempPknot)
{     
  //main seeking trackers
  //small to large trackers

  *tempPknot.small_behind_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.small_behind_y_trackerList_Count = -1;
  *tempPknot.small_front_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.small_front_y_trackerList_Count = -1;
  
  //large to small trackers
  *tempPknot.large_behind_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.large_behind_y_trackerList_Count = -1;
  *tempPknot.large_front_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.large_front_y_trackerList_Count = -1;
  
  
  //gap and trackers 
  *tempPknot.gapNucs_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.gapNucs_trackerList_Count = -1;
  *tempPknot.pairsNucs_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.pairsNucs_trackerList_Count = -1;


 int indexToAdd = -1;
  //initialize nucs front list
  for (int index = 0; index < nucsLenght, nucsLenght++)
  {
    *tempPknot.small_behind_y_trackerList[index] = -1;  
    *tempPknot.small_front_y_trackerList[index] = index;
    
    //large to small trackers
    *tempPknot.large_behind_y_trackerList[index] = -1;  
    *tempPknot.large_front_y_trackerList[index] = index;  
    
    //gap and trackers 
    *tempPknot.gapNucs_trackerList[index] = -1; 
    *tempPknot.pairsNucs_trackerList[index] = -1;  
  }

  //utils
  *tempPknot.nucsLenght = -1;
  *tempPknot.fullSequenceLenght = thefold->seqlength;
  *tempPknot.thisSegmentLenght = -1;
  *tempPknot.currentNuc_y = -1;
  *tempPknot.currentNuc_d = -1;
  *tempPknot.endNuc_y = -1;
  *tempPknot.isPaired_current_y = FALSE;
  *tempPknot.is_yGRTd_current = FALSE;
  *tempPknot.is_dGRTy_current = FALSE;
  
  *tempPknot.nextNuc_y = -1;
  *tempPknot.nextNuc_d = -1;
  *tempPknot.isPaired_next_y = FALSE;
  *tempPknot.is_yGRTd_next = FALSE;
  *tempPknot.is_dGRTy_next = FALSE;
  
  //loop, bulge, stack, pknot boolean logic variables
  *tempPknot.inGap = FALSE;

  *tempPknot.isLoop_suspected = FALSE;
  *tempPknot.isLoop_confident = FALSE;
  *tempPknot.isLoop_confirmed = FALSE;

  *tempPknot.isBulge_suspected = FALSE;
  *tempPknot.isBulge_confident = FALSE;
  *tempPknot.isBulge_confirmed = FALSE; 
 
  *tempPknot.isStack_suspected = FALSE;
  *tempPknot.isStack_confident = FALSE;
  *tempPknot.isStack_confirmed = FALSE;

  *tempPknot.isPknot_suspected = FALSE;
  *tempPknot.isPknot_confident = FALSE;
  *tempPknot.isPknot_confirmed = FALSE;
};

void AddNuc_TrackerList(int y, int *trackerList, int *trackerList_Count)
{
  int indexToAdd = segmentLenght-trackerList_Count;
  *trackerList[indexToAdd] = y;
  *trackerList_Count++;  
};

void RemoveNuc_TrackerList(int y, int *trackerList, int *trackerList_Count)
{
  //now remove from the front
  for (int index = 0; index < trackerList_Count, index++)
  {
    if (trackerList[index] == y)
    {
      *trackerList[index] = -1;
      *trackerList_Count--;
    }
  } 
};

bool NucInTrackerList(int nucIndexNum, int trakerList, int trackerList_Count)
{
 bool inTrackerList = FALSE;
    for (int index = 0; index < trackerList_Count; index++)
    {
      if (trakerList[index]==nucIndexNum)
      {
        //its in the listof expected front nucs so its valid ad should be able to loop to next nuck
        inTrackerList = TRUE;     
      }
    }

  return inTrackerList;
};

void SetStructureCondidence(bool suspected, bool confident, bool confirmned,
                           bool *suspectedTracker, bool *confidentTracker, bool *confirmnedTracker)
{
  *suspectedTracker =  suspected;
  *confidentTracker = confident;
  *confirmnedTracker = confirmned;
};

bool WalkAndTest_Structure(int startNuc_y, int endNuc_y, bool doSmallToLargeNuc, struct PknotDetectionData *pknotData_mainStruct, fold *thefold)
{
  //This is the logic for walking and testing the nucc pairing in a secondary structure
  //feed it a starting y nuc and it will walk the structure eacah nuc 1 at a time and test if paired
  //this first bit initializes a bunch of stuff and tehn it checks if curreny y is paired
  bool doSmall = doSmallToLargeNuc;
  int fullSeqLength = -1;
  int segmentLength = -1;
  if (doSmallToLargeNuc == TRUE)
  {
    *pknotData_mainStruct->thisSegmentLenght = (endNuc_y - startNuc_y);
    
  }
  else
  {
    *pknotData_mainStruct->thisSegmentLenght = (startNuc_y - endNuc_y);
  }

  segmentLength = pknotData_mainStruct->thisSegmentLenght;

  *tempPknot->endNuc_y = endNuc_y;
  *pknotData_mainStruct->currentNuc_y = startNuc_y;
  *pknotData_mainStruct->currentNuc_d = thefold->pairs[startNuc_y];
  
  *pknotData_mainStruct->nextNuc_y = startNuc_y+1;
  *pknotData_mainStruct->nextNuc_d = thefold->pairs[startNuc_y];
  
  
  *pknotData_mainStruct->nucsLenght = thefold->seqlength;
  fullSeqLength = pknotData_mainStruct->nucsLenght;


  if(pknotData_mainStruct->currentNuc_y == -1)
  {
    *pknotData_mainStruct->isPaired_current_y= FALSE;
  }
  else
  {
    *pknotData_mainStruct->isPaired_current_y= TRUE;
    if (currentNuc_y > currentNuc_d)
    {
      *pknotData_mainStruct->is_yGRTd_current = TRUE;
      *pknotData_mainStruct->is_dGRTy_current = FALSE;
      
    }
    else
    {
      *pknotData_mainStruct->is_yGRTd_current = FALSE;
      *pknotData_mainStruct->is_dGRTy_current = TRUE;
    }
  }

  if(pknotData_mainStruct->nextNuc_y == -1)
  {
    pknotData_mainStruct->isPaired_next_y= FALSE;
  }
  else
  {
    pknotData_mainStruct->isPaired_next_y= TRUE;
     

    if (currentNuc_y > currentNuc_d)
    {
      *pknotData_mainStruct->is_yGRTd_next = TRUE;
      *pknotData_mainStruct->is_dGRTy_next = FALSE;
    }
    else
    {
      *pknotData_mainStruct->is_yGRTd_next = FALSE;
      *pknotData_mainStruct->is_dGRTy_next = TRUE;
    }
  }
  
  // everything is initialized in the structure for the walk

    //add the last nuc to the list
    //add to the list immediatly when it is hit
    
  AddNuc_TrackerList(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->small_behind_y_trackerList,
                    pknotData_mainStruct->small_behind_y_trackerList_Count);

  if (inGap==TRUE)
  {
    //indexToAdd = nucsLenght-gapNucs_Count;
    //gapNucs[indexToAdd]=actualNuc;
  }

  if (inStack==TRUE)
  {
    //placeholder for this?
  }

  //now remove from the front
  RemoveNuc_TrackerList(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->small_front_y_trackerList
                    pknotData_mainStruct->small_front_y_trackerList_Count);
  
  //now we know what should be behind and ahead and what has been found in its gap if applicable
  //make a test data structure
  struct PknotDetectionData pknotData_testStruct;
  InitalizePknotStruct(*thefold, *pknotData_testStruct);
  *pknotData_testStruct.currentNuc_y=pknotData_mainStruct->currentNuc_y;
  *pknotData_testStruct.currentNuc_d=pknotData_mainStruct->currentNuc_d;
  *pknotData_testStruct.fullSequenceLenght=pknotData_mainStruct->fullSequenceLenght;
  *pknotData_testStruct.thisSegmentLenght=0;

  if (pknotData_mainStruct.isPaired_current_y==TRUE)
  {
    //this nuc is paired so record in main structure tracker as well as the primary test tracker
    AddNuc_TrackerList(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->pairsNucs_trackerList,
                    pknotData_mainStruct->pairsNucs_trackerList_Count)
    
    AddNuc_TrackerList(pknotData_testStruct->currentNuc_y, pknotData_testStruct->pairsNucs_trackerList,
                    pknotData_testStruct->pairsNucs_trackerList_Count);
    
    //now need to see if next nuc is a gap, stack, or pknot in first part of identifying structure

    //first check if the current nucs pair is in the current list of nucs in front
    bool isValid_inFront = FALSE;
    isValid_inFront = NucInTrackerList(pknotData_mainStruct.currentNuc_d, pknotData_mainStruct->small_front_y_trackerList,
                                      pknotData_mainStruct->small_front_y_trackerList_Count);

    if (isValid_inFront==TRUE)
    {
      //we are progressing through a normal sequence and stacks and it looks normal
      //if in front is valid then its most likely not a pknot
      //check next nuc pair

      //this is potentialy a stack
      SetStructureCondidence(TRUE, FALSE, FALSE,
                            pknotData_mainStruct->isStack_suspected, 
                            pknotData_mainStruct->isStack_confident, 
                            pknotData_mainStruct->isStack_confirmed);
      
      SetStructureCondidence(FALSE, FALSE, FALSE,
                            pknotData_mainStruct->isPknot_suspected, 
                            pknotData_mainStruct->isPknot_confident, 
                            pknotData_mainStruct->isPknot_confirmed);
      
    }
    else
    {
      //this is potentially in a pknot 
      SetStructureCondidence(FALSE, FALSE, FALSE,
                            pknotData_mainStruct->isStack_suspected, 
                            pknotData_mainStruct->isStack_confident, 
                            pknotData_mainStruct->isStack_confirmed);
      
      SetStructureCondidence(TRUE, FALSE, FALSE,
                            pknotData_mainStruct->isPknot_suspected, 
                            pknotData_mainStruct->isPknot_confident, 
                            pknotData_mainStruct->isPknot_confirmed);
    }
    
    //at this point we have suspicions only about the nucs\
    //potential structures are LOOP, STACK, BULGE, PKNOT
    //potential configs of individual nuc is UNPAIRD, PAIRED
    //if in LOOP each nuc after a stack can be another stack nuc but each jump will always be in the original front list
    //if in STACK then each y jump will have a d jump that is equivalent
    //if in PKNOT each jump should give a wild number that may or may not return to orignal front list. that is why need reverse search in this case is suspected
    
    //test for LOOP, STACK, BULGE
    if (pknotData_mainStruct->isStack_suspected==TRUE)
    {
      //test for stack as it is suspected to be the pair that paw detected during walk
      TestAfterPair_Stack(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->currentNuc_d, pknotData_mainStruct, )
    }

    if (pknotData_mainStruct->isPknot_suspected==TRUE)
    {
      //test for pknot
    }
    



    //now get next nuc info as part of test routine to determine next course of action
    int nextSearchNuc=searchNuc_index+1;
    int nextSearchNuc_compPair = thefold->pairs[nextSearchNuc];

    if (isStack_suspected==TRUE)
    {
      

    }
    else
    {
      //no point in testing next pair but need a good reason why
    }
    
  }  
};

bool TestAfterPair_Loop(int start_y, int start_d, struct PknotDetectionData *pknotData_mainStruct, struct PknotDetectionData *pknotData_testStruct, fold *thefold)
{
  //is LOOP if both i_first+1 and j_first-1 are not paired it is a GAPand any number of loops and that is all I care about as a gap means it is a new strcuture domain
  bool isLoop = FALSE;
  
  //we know that start_y i spaired so assume that
  int i_second = start_y+1;
  int j_second = start_d-1;
  if (i_second == -1 && j_second == -1)
  {
    isLoop = TRUE
  }
  return isLoop;
}

bool TestAfterPair_Bulge(int start_y, int start_d, struct PknotDetectionData *pknotData_mainStruct, struct PknotDetectionData *pknotData_testStruct, fold *thefold)
{
  //is BULGE if i_second IS paired but j_second IS NOT or i_second IS NOT paired but j_second IS  
  bool isBulge = FALSE;

  //we know that start_y i spaired so assume that
  int i_second = start_y+1;
  int j_second = start_d-1;
  if ((i_second != -1 && j_second == -1) || (i_second == -1 && j_second != -1))
  {
    isBulge = TRUE
  }
  return isBulge;
}

bool TestAfterPair_Stack(int start_y, int start_d, fold *thefold)
{
  bool isStack=FALSE;
  
  //this fucntion tests if teh nuc (y) passed is part of a stack of 2 or more base pairs
  //such that the first pair is i,j and the second pair is i+1, j-1
  //startNuc_y is the first nuc in what is suspected to be a stack
  bool isPaired_first = FALSE;
  bool isPaired_second = FALSE;
  
  int i_first = start_y;
  int j_first = start_d;
  if (i_first != -1)
  {
    isPaired_first=TRUE;
  }  

  int i_second =  i_first+1;
  int j_second = thefold->pairs[i_second];
  if (i_second != -1)
  {
    isPaired_second = TRUE
  }

  //if both nucs are paired then it might be in a stack
  //if they are not both paired they can not be a stack
  if( isPaired_first==TRUE && isPaired_second==TRUE)
  {
    //might be a stack
    //make sure all nucs pairs point to the right nucs
    //then check the math
    if (thefold->pairs[j_first]==i_first && thefold->pairs[j_second]==i_second)
    {
      //both y and y+1 has a pair that matches a reverse query of the j nucs
      //now test the nuc order
      if (i_second==i_first+1 && j_second == j_first-1)
      {
        //it is a stack yay
        isStack=TRUE;
      }
    }
  }
  else
  {
    //it cant be a stack
    //do nothing as initialized to FALSE
    
  }
  return isStack;
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

  //logic for deciding if stack or pknot
  bool isStack_suspected = FALSE;
  bool isStack_confident = FALSE;
  bool isPknot_suspected = FALSE;
  bool isPknot_confident = FALSE;


  while (noPair_NextNuc == TRUE)
  { 
   
    /

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


