
#ifndef PKNOT_ETERNA_H
#define PKNOT_ETERNA_H

#ifdef __cplusplus
extern "C" {
#endif

#include "pfuncUtilsConstants.h"
#include "runtime_constants.h"
#include "physical_constants.h"
#include "../../shared/utilsHeader.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

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

bool WalkAndTest_Structure(int startNuc_y, int endNuc_y, bool doSmallToLargeNuc, struct PknotDetectionData *pknotData_mainStruct, fold *thefold);
bool TestAfterPair_LoopOrGap(int start_y, int start_d, fold *thefold);
bool TestAfterPair_Bulge(int start_y, int start_d, int *unpairedNucStart, fold *thefold);
bool TestAfterPair_Stack(int start_y, int start_d, fold *thefold);
void InitalizePknotStruct(fold *thefold, struct PknotDetectionData *tempPknot);
void AddNuc_TrackerList(int y, int *trackerList, int *trackerList_Count);
void RemoveNuc_TrackerList(int y, int *trackerList, int *trackerList_Count);
bool NucInTrackerList(int nucIndexNum, int *trakerList, int trackerList_Count);
void ResetTracker_Back(int *trackerList, int *trackerList_Count,  int sequenceLength);
void ResetTracker_Front(int *trackerList, int *trackerList_Count, int sequenceLength);
void SetStructureCondidence(bool suspected, bool confident, bool confirmned,
                           bool *suspectedTracker, bool *confidentTracker, bool *confirmnedTracker);
void LogFrontBackTrackers(int current_y, struct PknotDetectionData *testStructDatad);
bool TestAfterPairFound_IsPknot(int start_y, int start_d, struct PknotDetectionData *mainStructData, fold *thefold);

#ifdef __cplusplus
}
#endif

#endif
