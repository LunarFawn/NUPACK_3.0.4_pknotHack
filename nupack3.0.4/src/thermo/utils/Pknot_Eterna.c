/* 
    Pknot_Eterna.c is part of modifactions to NUPACK 3.0.4 to allow for extended abilty to identify and process 
    pseudoknots outside of the limited subset designated by Nupack in their paper.
    It attempts to identify and label structures as it walkes a secondary structure and allow for identification of
    any pknot type it might encounter. This is being developed and maintained by volunteers in the Eterna community 
    and the intial algorithim for this new process was developed by Jennifer Pearl .
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

bool WalkAndTest_Structure(int startNuc_y, int endNuc_y, bool doSmallToLargeNuc, struct PknotDetectionData *pknotData_mainStruct, fold *thefold)
{
  bool isPknot_finalAnswer = FALSE;
  bool PknotFoundInSequence= FALSE;
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

  *pknotData_mainStruct->endNuc_y = endNuc_y;
  *pknotData_mainStruct->currentNuc_y = startNuc_y;
  *pknotData_mainStruct->currentNuc_d = thefold->pairs[startNuc_y];
  
  *pknotData_mainStruct->nextNuc_y = startNuc_y+1;
  *pknotData_mainStruct->nextNuc_d = thefold->pairs[startNuc_y];
  
  
  *pknotData_mainStruct->nucsLenght = thefold->seqlength;
    fullSeqLength = pknotData_mainStruct->nucsLenght;

    while (pknotData_mainStruct->currentNuc_y<=pknotData_mainStruct->endNuc_y)
    {
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

        bool foundStructureType=FALSE;
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

                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isBulge_suspected, 
                                        pknotData_mainStruct->isBulge_confident, 
                                        pknotData_mainStruct->isBulge_confirmed);
                    
                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isLoop_suspected, 
                                        pknotData_mainStruct->isLoop_confident, 
                                        pknotData_mainStruct->isLoop_confirmed);
                    foundStructureType = TRUE;
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

                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isStack_suspected, 
                                        pknotData_mainStruct->isStack_confident, 
                                        pknotData_mainStruct->isStack_confirmed);

                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isLoop_suspected, 
                                        pknotData_mainStruct->isLoop_confident, 
                                        pknotData_mainStruct->isLoop_confirmed);
                    foundStructureType = TRUE;
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
                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isBulge_suspected, 
                                        pknotData_mainStruct->isBulge_confident, 
                                        pknotData_mainStruct->isBulge_confirmed);
                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isStack_suspected, 
                                        pknotData_mainStruct->isStack_confident, 
                                        pknotData_mainStruct->isStack_confirmed);
                    foundStructureType = TRUE;
                }
                else
                {
                    SetStructureCondidence(FALSE, FALSE, FALSE,
                                        pknotData_mainStruct->isLoop_suspected, 
                                        pknotData_mainStruct->isLoop_confident, 
                                        pknotData_mainStruct->isLoop_confirmed);
                }
            }
            
            
            if (pknotData_mainStruct->isPknot_suspected==TRUE || foundStructureType == FALSE;)
            {
                //test for pknot      
                isPknot_finalAnswer = TestAfterPairFound_IsPknot(pknotData_mainStruct->currentNuc_y, pknotData_mainStruct->currentNuc_d, pknotData_mainStruct, thefold );

                if (isPknot_finalAnswer==TRUE)
                {
                    PknotFoundInSequence=TRUE;
                    *thefold->actualPknots[pknotData_mainStruct->currentNuc_y]=pknotData_mainStruct->currentNuc_d; 
                }
            }                
        } 
        else
        {
            //current_y was not paired so go to the next nuc            
        }

        //now go to next y in sequnce for next nuc inspection
        *pknotData_mainStruct->currentNuc_y++;
    }  
  //need to change this to now walk the whole structre and log pknots then pass whether one was found or not
  return PknotFoundInSequence; 
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

bool NucInTrackerList(int nucIndexNum, int *trakerList, int trackerList_Count)
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

bool TestAfterPairFound_IsPknot(int start_y, int start_d, struct PknotDetectionData *mainStructData, fold *thefold)
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

    //new thought 6-25-22
    //walk until gap is hi and record. now begin walkig the main struct to look for the infamous stack that is the source of all our troubles
    //if a stack is there and it influences the trackers wrong then it is a pknot. basically need to walk the main struct now disregarding all pairs but
    //doing the jumps to understand the loop and whatnot. stop either when the frist stack is hit or the logic has looped back on teh original start_y that was passed at time of call
    
    //now walk and test
    //first initialize walker logic pairs
    //we know start_y and d are pairs
    //walker_y is the nuc after the jump
    //now work with the test structure data

    struct PknotDetectionData *mainPknotStructData
    struct PknotDetectionData *testPknotStructData

    *mainPknotStructData = mainStructData;
    *testPknotStructData = mainStructData;

    *testStructData->currentNuc_y = start_d;
    *testStructData->currentNuc_d = thefold->pairs[start_d];
    

    //new thought 6-25-22
    //make single walk of walker_y and record now step start_y one step and stop untul wlaker_d is hit or a pair is found.
    // the pair is jumped and the section betwwen but not including y to d of this search is removed from a test tracker for the main seq front tracker.
    //walker_y is now investigated.
    //if hit walker_y then the pair should  just be a normal pair 
    //if walker_y is between start y and start d then it is a pknot
    //if start_y pair as it is incremented adn tested is paired and the d for it is less than walker_d then it is a pkont
    //loop until find a gap
    
    int walker_y_entrance = start_d;
    int walker_d_entrance = thefold->pairs[walker_y_entrance];

    int walker_y_oneStep = walker_y_entrance+1;
    int walker_d_oneStep = thefold->pairs[walker_y_oneStep];

    //now walk nucs from start_y
    int mainStruct_y=start_y;
    int mainStruct_d=start_d;
    bool stopWalkingMainStruct = FALSE;

    bool isPair = FALSE;
    bool isWalker_y_oneStep=FALSE;

    
    while (stopWalkingMainStruct==FALSE)
    {
      //advance the main struct nuc 1 step at a time until a pair is hit or walker_y_onstep is
      mainStruct_y++;
      mainStruct_d = thefold->pairs[mainStruct_y];
    
      //
    
      //is this a pair or onestep walker y
    

      //if have his the onesetp walker than you have gone all the way around and the main struct is in the same loop as the test nuc
      //this means that it is not a pknot
      if(mainStruct_y == walker_y_oneStep)
      {
        isWalker_y_oneStep=TRUE;   
        stopWalkingMainStruct=TRUE;     
      }
      else
      {
        if(mainStruct_y != -1)
        {
          isPair=TRUE;
          stopWalkingMainStruct=TRUE;  
        }
        else
        {
          //if its not a pair then keep going
        }
      }
    }

    //now we know what all the nuc numbers are of those in play and we need to do logic on them
    bool inTracker = FALSE;
    
    if (isWalker_y_oneStep==TRUE)
    {
      //it went all the way around and came to the walker form teh main struct so it is expected and is not pknot
      isPknot=FALSE;
    }
    else
    {
      if (isPair==TRUE)
      {
        int index_y=-1;
        for (index_y=mainStruct_y+1; index_y < mainStruct_d; index_y++)
        {
          //remove all the nucs betwwen start_y and start_d from front tracker for main struct 
          RemoveNuc_TrackerList(index_y, mainPknotStructData->small_front_y_trackerList, mainPknotStructData->small_front_y_trackerList_Count);
        }

        inTracker = NucInTrackerList(walker_y_oneStep, mainPknotStructData->small_front_y_trackerList,
                                      mainPknotStructData->small_front_y_trackerList_Count);
        if (inTracker==FALSE)
        {
          //its not in the tracker so it is not expected and thus it is a pknot
          isPknot=TRUE;
        }
      }
    }
  }
  else
  {
    //if its not valid then it be default means teh fold is not in the normal folding prgression of a zipper and the fold in question is bending backward which is kinda teh definition of pknot
    isPknot = TRUE;
  }

  return isPknot;
};
   
