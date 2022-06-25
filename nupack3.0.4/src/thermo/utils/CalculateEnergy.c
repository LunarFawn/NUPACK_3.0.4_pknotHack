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

      SetStructureCondidence(TRUE, FALSE, FALSE,
                            pknotData_mainStruct->isBulge_suspected, 
                            pknotData_mainStruct->isBulge_confident, 
                            pknotData_mainStruct->isBulge_confirmed);
      
      SetStructureCondidence(TRUE, FALSE, FALSE,
                            pknotData_mainStruct->isLoop_suspected, 
                            pknotData_mainStruct->isLoop_confident, 
                            pknotData_mainStruct->isLoop_confirmed);
      
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
      
      SetStructureCondidence(FALSE, FALSE, FALSE,
                            pknotData_mainStruct->isBulge_suspected, 
                            pknotData_mainStruct->isBulge_confident, 
                            pknotData_mainStruct->isBulge_confirmed);
      
      SetStructureCondidence(FALSE, FALSE, FALSE,
                            pknotData_mainStruct->isLoop_suspected, 
                            pknotData_mainStruct->isLoop_confident, 
                            pknotData_mainStruct->isLoop_confirmed);
      
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
    bool structureFound=FALSE;
    if (pknotData_mainStruct->isStack_suspected==TRUE)
    {
      //test for stack as it is suspected to be the pair that paw detected during walk
      bool isStack = TestAfterPair_Stack(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->currentNuc_d);
      
      if (isStack==TRUE)
      {
        SetStructureCondidence(TRUE, TRUE, FALSE,
                              pknotData_mainStruct->isStack_suspected, 
                              pknotData_mainStruct->isStack_confident, 
                              pknotData_mainStruct->isStack_confirmed);
      }
      else
      {
        SetStructureCondidence(FALSE, FALSE, FALSE,
                              pknotData_mainStruct->isStack_suspected, 
                              pknotData_mainStruct->isStack_confident, 
                              pknotData_mainStruct->isStack_confirmed);
      }
    }

    if (pknotData_mainStruct->isBulge_suspected==TRUE)
    {
      //test for stack as it is suspected to be the pair that paw detected during walk
      bool isBulge = TestAfterPair_Bulge(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->currentNuc_d);
      
      if (isBulge==TRUE)
      {
        SetStructureCondidence(TRUE, TRUE, FALSE,
                              pknotData_mainStruct->isBulge_suspected, 
                              pknotData_mainStruct->isBulge_confident, 
                              pknotData_mainStruct->isBulge_confirmed);
      }
      else
      {
        SetStructureCondidence(FALSE, FALSE, FALSE,
                              pknotData_mainStruct->isBulge_suspected, 
                              pknotData_mainStruct->isBulge_confident, 
                              pknotData_mainStruct->isBulge_confirmed);
      }
    }

    if (pknotData_mainStruct->isLoop_suspected==TRUE)
    {
      //test for stack as it is suspected to be the pair that paw detected during walk
      bool isLoop = TestAfterPair_Loop(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->currentNuc_d);
      
      if (isLoop==TRUE)
      {
        SetStructureCondidence(TRUE, TRUE, FALSE,
                              pknotData_mainStruct->isLoop_suspected, 
                              pknotData_mainStruct->isLoop_confident, 
                              pknotData_mainStruct->isLoop_confirmed);
      }
      else
      {
        SetStructureCondidence(FALSE, FALSE, FALSE,
                              pknotData_mainStruct->isLoop_suspected, 
                              pknotData_mainStruct->isLoop_confident, 
                              pknotData_mainStruct->isLoop_confirmed);
      }
    }
    

    if (pknotData_mainStruct->isPknot_suspected==TRUE)
    {
      //test for pknot
    }
    



    
    
  }  
};

bool TestAfterPair_LoopOrGap(int start_y, int start_d, fold *thefold)
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
};

bool TestAfterPair_Bulge(int start_y, int start_d, int *unpairedNucStart, fold *thefold)
{
  //is BULGE if i_second IS paired but j_second IS NOT or i_second IS NOT paired but j_second IS  
  bool isBulge = FALSE;

  //we know that start_y i spaired so assume that
  int i_second = start_y+1;
  int j_second = start_d-1;
  if (i_second != -1 && j_second == -1)
  {
    unpairedNucStart=j_second;
    isBulge = TRUE
  }
  else if (i_second == -1 && j_second != -1)
  {
    unpairedNucStart=i_second;
    isBulge = TRUE
  }
  return isBulge;
}; 

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
};

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
  bool isValidFront_currentNuc_y;
  int currentNuc_d;
  bool isValidFront_currentNuc_d;
  int endNuc_y;
  bool isPaired_current_y;
  bool is_yGRTd_current;
  bool is_dGRTy_current;
  
  int nextNuc_y;
  bool isValidFront_nextNuc_y;
  int nextNuc_d;
  bool isValidFront_nextNuc_d;
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
  *tempPknot.small_behind_y_trackerList_Count = 0;
  *tempPknot.small_front_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.small_front_y_trackerList_Count = nucsLenght;
  
  //large to small trackers
  *tempPknot.large_behind_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.large_behind_y_trackerList_Count = 0;
  *tempPknot.large_front_y_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.large_front_y_trackerList_Count = nucsLenght;
  
  
  //gap and trackers 
  *tempPknot.gapNucs_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.gapNucs_trackerList_Count = 0;
  *tempPknot.pairsNucs_trackerList = (int*) calloc( thefold->seqlength, sizeof(int));
  *tempPknot.pairsNucs_trackerList_Count = 0;


 int indexToAdd = -1;
  //initialize nucs front list
  for (int index = 0; index < nucsLenght, index++)
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
  *tempPknot.isValidFront_nextNuc_y = FALSE;
  *tempPknot->isValidFront_nextNuc_d=FALSE;
  *tempPknot->isValidFront_currentNuc_y=FALSE;
  *tempPknot->isValidFront_currentNuc_d=FALSE;
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

void ResetTracker_Back(int *trackerList, int *trackerList_Count,  int sequenceLength)
{
  for (int index =0; index < sequenceLength; index++)
  {
    *trackerList[index]=-1;
  }
  *trackerList_Count = 0; 
};

void ResetTracker_Front(int *trackerList, int *trackerList_Count, int sequenceLength)
{
  for (int index =0;index < sequenceLength;index++)
  {
    *trackerList[index]=index+1;
  }
  *trackerList_Count = sequenceLength;
};

void SetStructureCondidence(bool suspected, bool confident, bool confirmned,
                           bool *suspectedTracker, bool *confidentTracker, bool *confirmnedTracker)
{
  *suspectedTracker =  suspected;
  *confidentTracker = confident;
  *confirmnedTracker = confirmned;
}; 


void LogFrontBackTrackers(int current_y, struct PknotDetectionData *testStructDatad)
{
  //first go through and log front tracekrs
  //start nuc is the begining of the struct or the last nuc tested
  //current  nuc is the nuc that you are currently at now
 
  RemoveNuc_TrackerList(current_y, testStructData->large_front_y_trackerList, testStructData->large_front_y_trackerList_Count); 

  AddNuc_TrackerList(current_y, testStructData->large_front_y_trackerList, testStructData->large_front_y_trackerList_Count);  
}

bool TestAfterPairFound_IsPknot(int start_y, int start_d, struct PknotDetectionData *mainStructData, struct PknotDetectionData *testStructData, fold *thefold)
{
  //the concept here is that a pknot basically can be thought of as a wormhole for RNA. Each loop is region in space and a pknot linkes those loops together
  //if the pair is in a loop then it just forms stacks but the pairs transect loops for lack of a better word
  //thus if you have traversed into a pknot you should not have any information about where you are and you fron and back estimates will be off when you poke around
  //and compare a few different itterations. You go past each hit for a stack each iteration as well as going in the opposite direction at the same time. You shoud eventually
  //be able to say it is a pknot quickly or you end up building the whole strucutre and you can then make a very good educated decision.

  //pair has been found so now first test the complentary nuc for the nuc under test and this is start_d
  //test if in the forward tracker

  //when find a pair travers it and record. now test walker continues on and ignores any pairs that make sense based on the pair found. if I finds a pair that does not match then it knows
  //there is pknot between one of them.

  bool isPknot = FALSE;

  bool isValid_start_d = NucInTrackerList(start_d, mainStructData->small_front_y_trackerList, mainStructData->large_front_y_trackerList_Count);

  if (isValid_start_d==TRUE)
  {
    //no make sure that it is not a pknot that bind s just before a stack so it does look like it is valid in front but there is a stack later that makes it not valid
    //idea
    //do a pknot test jump. the front should now be 21 and above now for exampleif y was 10. so you say ok what if I go to the next 
    //nuc in the number if its a gap and the num is not in the the test front (instead nowin the back) but in the main seq front 
    //then you know you just hit a pknot. now also if you hit a good stack then you may just be in a pknot so you go until you hit a gap. 
    //if you were in a pknot as before you would now be in a unique front state. f you were not in a pknot then the first
    //gap you hit teh nuc should match the test back tracker and not be in teh front or in main sequence front but 
    //should be I think in the main seq back now as well
    *testStructData = mainStructData;

    
    //now walk and test
    //first initialize walker logic pairs
    //we know start_y and d are pairs
    //walker_y is the nuc after the jump
    //now work with the test structure data
    *testStructData->currentNuc_y = start_d;
    *testStructData->currentNuc_d = thefold->pairs[start_d];

    //loop until find a gap
    int walker_y = start_d;
    int walker_d = thefold->pairs[walker_y];

    bool inGapNow = FALSE;

    //setup trackers
    ResetTracker_Front(testStructData->small_front_y_trackerList,
                       testStructData->large_front_y_trackerList_Count, testStructData->fullSequenceLenght);

    ResetTracker_Back(testStructData->small_behind_y_trackerList,
                       testStructData->small_behind_y_trackerList_Count, testStructData->fullSequenceLenght);
    
    //npw set front and back trackers to current Y poisiton after the jump
    int index_y = 0;
    int index_d = 0;
    for (int index_y = 0; index_y <= start_y; index_y++)
    {
      LogFrontBackTrackers(index_y, testStructData);
    }

    //now that all the nucs for front and back are logged upto teh start_y before teh jump we now log the jumped nuck
    LogFrontBackTrackers(start_d, testStructData);

    bool inFrontTracker=FALSE;
    //now record front and back
    while (inGapNow==FALSE)
    {
      walker_y++;
      walker_d = thefold->pairs[walker_y];
      //test if walker nucs are paired and figure out if walker_d is in the front
      if(walker_d != -1)
      {
        inFrontTracker = NucInTrackerList(walker_d, testStructData->large_front_y_trackerList,
                                          testStructData->large_front_y_trackerList_Count);

        if (inFrontTracker==TRUE)
        {
          //this was a good gap and keep going
        }
        else
        {
          //this is not in the front tracker so it might 
        }
      }
      else
      {
        inGapNow=TRUE;
      }

    }
  }
  else
  {
    //if its not valid then it be default means teh fiold is nopt in the normal folding prgression of a zipper and the fold in question is bending backward which is kinda teh definition of pknot
    isPknot = TRUE;
    return isPknot;
  }
};


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
  for (int index = 1; index <= nucsLenght; nucsLenght++)
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


