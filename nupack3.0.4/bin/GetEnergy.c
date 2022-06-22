/* 
   GetEnergy.c  by Robert Dirks.  

   This program determines the energies of
   substructures, mirroring the loops of Fold.out
   and includes the function for determining energies.  03/15/2001  
   
*/
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "pfuncUtilsHeader.h"
#include "DNAExternals.h"


/* ******************************** */
DBL_TYPE GetEnergy( fold *thefold) {
  DBL_TYPE energy;
  
  LoadEnergies();
  energy = EnergyF( 0, thefold->seqlength - 1, thefold);
  
  return energy;
}

/* ********************************* */
DBL_TYPE EnergyF( int start, int stop, fold *thefold) {
  
  DBL_TYPE energy = 0.0;
  int d; //Left end of rightmost pair or pk
  
  int j;
  DBL_TYPE bp_penalty;
  
  j = stop; 
  while( j >= start) {
    
    if( thefold->isNicked[j]) {
#ifdef STRUCTURE_WARNINGS
      printf("Warning, disconnected complex for:%s\n", thefold->seq);
#endif
        return NAD_INFINITY;
    }
    
    if( thefold->pknots[ j] != -1) {
      d = thefold->pknots[ j];      
      
      energy +=
        EnergyPk( d, j, thefold) + BETA_1 +
        DangleEnergyWithPairs( j+1, stop, thefold);
      
      j = d-1;
      stop = j;
      
    }
    else if( thefold->pairs[ j] != -1) {
      d = thefold->pairs[ j];
      if( d > j) {
        printf("Error: Unclassified pknot!\n");
        exit(1);
      }
      bp_penalty = 0;
      if( thefold->seq[d] != BASE_C && thefold->seq[j] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }
      //EnergyFb( d, j, thefold);
      energy +=
        EnergyFb( d, j, thefold) + 
        DangleEnergyWithPairs( j+1, stop, thefold) +
        bp_penalty;
      
      j = d-1;
      stop = j;
    }
    else {
      j--;
    }
  }
  
  
  energy += DangleEnergyWithPairs( start, stop, 
                                  thefold);
  /*  
  energy += DangleEnergy( start, stop, 
  thefold->seq, thefold->seqlength);
  */  
  
  return energy;
}

/* *************************************** */

//this is a hack of EnergyFb that will ignore a stack between a large interior loop and will just compute the energies as if
//the two loops were a big main loop like figure 15
//Fb walkes the stacks and sums energy
DBL_TYPE EnergyFb( int start, int stop, fold *thefold) {
  char mystr[10];
  int *specialCasePknotStackStart_nuc1 = calloc((thefold->seqlength), sizeof( int));
  int *specialCasePknotStackStart_nuc2 = calloc((thefold->seqlength), sizeof( int));
  int *specialCasePknotStackEnd_nuc1 = calloc((thefold->seqlength), sizeof( int));
  int *specialCasePknotStackEnd_nuc2 = calloc((thefold->seqlength), sizeof( int));
  int numStacksFound = -1;

  bool use_SpecialCase_DoubleInteriorLoop = TRUE;

  DBL_TYPE energy = 0.0;
  int d=-1; //Left end of rightmost pair or pk
  
  int i, j;
  
  //int hairpin = TRUE;
  //int interiorLoop = FALSE;
  DBL_TYPE bp_penalty;
  int firstStop = stop;
  
  int p1, p2;
  
  int *pairs;
  int nPairs;
  int nNicks = 0;
  if(!CanPair(thefold->seq[start],thefold->seq[stop])) {
    return NAD_INFINITY;
  }
  
  pairs = (int*) malloc( (thefold->seqlength)*sizeof( int) );
  pairs[0] = start;
  pairs[1] = stop;
  nPairs = 1;
  for( i = 2; i <= thefold->seqlength - 1; i++) {
    pairs[i] = -1;
  }
  
  j = stop - 1; 
  stop = j;
  while( j >= start + 1) { //determine loop type, save pairs
    if( thefold->isNicked[j]) nNicks++;
    
    
    if( thefold->pknots[ j] != -1) {
      
      //hairpin = interiorLoop = FALSE;
      energy += 2*ALPHA_2;
      
      d = thefold->pknots[ j];
      
      
      if(!CanPair(thefold->seq[d],thefold->seq[thefold->pairs[d]])
          || !CanPair(thefold->seq[j], thefold->seq[thefold->pairs[j]])) {
        free(pairs);
        pairs = NULL;
        return NAD_INFINITY;
      }
      pairs[2*nPairs] =  d;
      pairs[2*nPairs + 1] = j;
      nPairs++;
      
      //pairs[2*nPairs] = j;
      //pairs[2*nPairs + 1] = thefold->pairs[j];
      
      j = d-1;
      stop = j;
    }
    else if( thefold->pairs[ j] != -1) {    
      d = thefold->pairs[ j];
  
      pairs[2*nPairs] =  d;
      pairs[2*nPairs + 1] = j;
      if(!CanPair(thefold->seq[d],thefold->seq[j])) {
        free(pairs);
        pairs = NULL;
        return NAD_INFINITY;
      }
      nPairs++;

      if( d > j) 
      {   
        //this is the special case code to do the doubel interiro loop with a pknot
        if (use_SpecialCase_DoubleInteriorLoop == TRUE)
        {
          
          //test if pair feeds into another loop and if so then go on next nuc in order after teh stack ends.
          //treat these stacks as unpaired for now

          int pairs_nuc1 = j;
          int nextPair_nuc1 = pairs_nuc1;
          int pairs_nuc2 = d;
          int nextPair_nuc2 = pairs_nuc2;

          bool stopChecking=FALSE;
          //we know that this is not the start nuc as it would have already stopped at the while loop at the start
          //so check if there is another nuc pair and if so then walk until you find the end and resume the rest of the logic
          //even if you can consider a pknot bewteen two independant internal loops a multiploop. the multiloop energy is only
          //determined by the base pairs at each branch and does not care about the stack energies.

          //make sure not to go past start so if the check shows teh next nuc is start nuc then bail and exit the pairs function and 
          //let the main while loop finish
          bool pastFirstPair=FALSE;
          
          while (stopChecking == FALSE)
          {
            nextPair_nuc1 = nextPair_nuc1-1;
            nextPair_nuc2 = thefold->pairs[nextPair_nuc1];

            if (nextPair_nuc1 <= start)
            {
              //you are at the end of this main loop so exit the parent pairs function
              //printf("I am at the end of the main loop!\n");
              //exit(1);
              j=nextPair_nuc1;
              stop=j;
              break;
            }
            else
            {
              //folow the nucs until you either come to start or you hit a gap
              //this step is recursive so need to call a recursive algorithm

              bool isPaired = FALSE;
              if ((thefold->pairs[thefold->pairs[nextPair_nuc1]]) == nextPair_nuc1)
              {
                //printf("I found a pair!\n");
                //exit(1);
                isPaired = TRUE;
                
              }

              if (isPaired == FALSE)
              {
                //this means the the next nuc is a gap and not the start so that is all we need to know so break out of this
                //and contniue walking the sequence
                //record the last nuc pair in the stack which is current nuc 1 minus 1
                specialCasePknotStackEnd_nuc1[numStacksFound]=nextPair_nuc1;
                specialCasePknotStackEnd_nuc2[numStacksFound]=nextPair_nuc2;
                stopChecking=TRUE;
                //printf("I stopped checking!\n");
                //exit(1);              
                j=nextPair_nuc1;
                stop=j;
                //char mystr[10];
                //printf("I should have exited teh pknot!\n");
                
                //exit(1);
              }
              else
              {
                //the next pair is paired so continue down the line checking
                if (pastFirstPair==FALSE)
                {                    
                  numStacksFound++;
                  specialCasePknotStackStart_nuc1[numStacksFound]=nextPair_nuc1;
                  specialCasePknotStackStart_nuc2[numStacksFound]=nextPair_nuc2;
                  pastFirstPair=TRUE;
                  nPairs--;
                }
                else
                {
                  //it is in the stack so just keep going through the loop until reach the end
                  //aka do nothing
                }                
              }
            }
          }          
        }
        else
        {            
          printf("Error: Unclassified pknot!\n");
          exit(1); 
          //j = d-1;
          //stop = j;
        }
      }
      else
      {
       
      }  
    }
    else {
      j--;
    }
  }
  if( thefold->isNicked[start]) nNicks++;
 
  
  if( nNicks >= 2) {
#ifdef STRUCTURE_WARNINGS	
    printf("Warning!, disconnected structure.\n");
#endif
      free(pairs);
      pairs = NULL;
      return NAD_INFINITY;
  }
  
  if( nNicks == 0) {
    if( nPairs == 1) { //hairpin
    sprintf(mystr, "%d", nPairs);
    printf(mystr);
    printf("\n");
    printf("I did a hairpin\n");
    //exit(1);
      energy = HairpinEnergy( start, firstStop, thefold->seq);			
      //return energy;
    }
    else if( nPairs == 2) { //interior loop
      sprintf(mystr, "%d", start);
      printf(mystr);
      printf("\n");
      printf("I did a interior\n");
      //exit(1);
      energy = InteriorEnergy( start, firstStop, pairs[2], pairs[3], thefold->seq)
        + EnergyFb( pairs[2], pairs[3], thefold);
      //return energy;
    }
    else if( nPairs >= 3) { //multiloop
      sprintf(mystr, "%d", nPairs);
      printf(mystr);
      printf("\n");
      printf("I did a multi\n");
      //exit(1);
      energy = ALPHA_1 + ALPHA_2 + ALPHA_3 * (pairs[1]-pairs[3]-1) + 
        DangleEnergyWithPairs(pairs[3]+1, pairs[1]-1, thefold);
      
      if( thefold->seq[ pairs[0]] != BASE_C && thefold->seq[ pairs[1]] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }
      else {
        bp_penalty = 0;
      }
      
      for( i = 1; i < nPairs; i++) {
        if( thefold->pknots[ pairs[2*i]] == -1) { //not a pseudoknot
          if( thefold->seq[ pairs[2*i]] != BASE_C && thefold->seq[ pairs[2*i+1]] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          
          energy += EnergyFb( pairs[2*i], pairs[ 2*i+1], thefold) + ALPHA_2;
          sprintf(mystr, "%d", nPairs);
          printf(mystr);
          printf("\n");
          printf("I did a not pseudo\n");
          //exit(1);
        }
        else { //a pseudoknot

           printf("I did a special\n");
           //exit(1);
          p1 = pairs[2*i];
          p2 = thefold->pairs[ p1];
          if( thefold->seq[ p1 ] != BASE_C && thefold->seq[ p2] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          
          p2 = pairs[2*i+1];
          p1 = thefold->pairs[ p2 ];
          if( thefold->seq[ p1 ] != BASE_C && thefold->seq[ p2] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          
          energy += EnergyPk( pairs[2*i], pairs[2*i + 1], thefold) + BETA_1M + 2*ALPHA_2;

          //now compensate for removing the stack by adding it and having the energy added
          //here start means as it is found durung pknot adn that is actually last in this case
        }	
        
        //reevaluate this
        p1 = 2*i;
        if( i != nPairs - 1) p2 = 2*i+3;
        else p2 = 0;
        
        energy += ALPHA_3*(pairs[p1] - pairs[ p2]-1) + 
          DangleEnergyWithPairs( pairs[ p2]+1, pairs[p1]-1, thefold);
      }

        //if (use_SpecialCase_DoubleInteriorLoop == TRUE)
        //{        
        //  printf("I did a special\n");
        //  //exit(1);
        //
        //  int special_end=specialCasePknotStackStart_nuc2[0];
        //  int special_start = specialCasePknotStackEnd_nuc2[0];
        //  //energy += EnergyFb( special_end, special_start , thefold);
        //  
        //   sprintf(mystr, "%d", energy);
        //  printf(mystr);
        //  printf("\n");
        //  printf("I did a energy\n");
        //  //exit(1);
        //}

      energy += bp_penalty;
    }
    else {
      printf("Error in calculation of EnergyFb! %d\n", nPairs);
      exit(1);
    }
  }
  else if( nNicks == 1) { //nNicks
    if( thefold->seq[ pairs[0]] != BASE_C && thefold->seq[ pairs[1]] != BASE_C) {
      bp_penalty = AT_PENALTY;
    }
    else {
      bp_penalty = 0;
    }
    
    if( nPairs == 1) //nicked hairpin
      energy = DangleEnergyWithPairs(pairs[0]+1, pairs[1]-1, thefold) + bp_penalty;
    else {
      energy = DangleEnergyWithPairs(pairs[3]+1, pairs[1]-1, thefold);
      
      for( i = 1; i < nPairs; i++) {
        if( thefold->seq[ pairs[2*i]] != BASE_C && thefold->seq[ pairs[2*i+1]] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }
        
        energy += EnergyFb( pairs[2*i], pairs[ 2*i+1], thefold);
        
        p1 = 2*i;
        if( i != nPairs - 1) p2 = 2*i+3;
        else p2 = 0;
        
        energy += DangleEnergyWithPairs( pairs[ p2]+1, pairs[p1]-1, thefold);
      }
      energy += bp_penalty;
      //return energy;
    }
    
  }
  
  //#ifdef SHOWFB
  //printf("start = %d stop = %d Fb = %f\n", start, firstStop, (double) energy);
  //#endif
  
  free( pairs);
  pairs = NULL;
  
  return energy;
}

/* *************************************************** */
DBL_TYPE EnergyPk( int i, int j, fold *thefold) {
  
  int a=0,b=0,c=0,d=0,e=0,f=0; 
  //along with i,j these are the 
  //eight key points of a pseudoknots (see figure 20 of paper)
  
  DBL_TYPE energy;
  int findCF, findAD;
  DBL_TYPE bp_penalty;
  int *seq = thefold->seq;
  int rightSingle;
  int leftSingle;
  
  //  printf("Pk! %d %d\n", i, j);
  
  e = thefold->pairs[i]; //guaranteed
  b = thefold->pairs[j]; //guaranteed
  
  rightSingle = FALSE;
  findCF = FALSE;
  f = e+1; //initial guess
  while( findCF == FALSE && rightSingle == FALSE) {
    if( thefold->pairs[ f] != -1 && thefold->pairs[f] < e) {
      c = thefold->pairs[ f];
      findCF = TRUE;
    }
    else {
      f++;
      if( f >= j) {
        rightSingle = TRUE;
        f = j;
        c= b;
      }
    }
  }
  
  leftSingle = FALSE;
  findAD = FALSE;
  a = b - 1; //initial guess
  while( findAD == FALSE) {
    
    if( thefold->pairs[ a] != -1 && thefold->pairs[ a] >  c) {
      d = thefold->pairs[ a];
      findAD = TRUE;
    }
    else {
      a--;
      if( a <= i) {
        leftSingle = TRUE;
        a = i;
        d = e;
      }
    }
  }
  
  bp_penalty = 0;
  if( seq[ thefold->pairs[ i]] != BASE_C && seq[ thefold->pairs[ e]] != BASE_C) {
    bp_penalty += AT_PENALTY;
  }
  if( seq[ thefold->pairs[ a]] != BASE_C && seq[ thefold->pairs[ d]] != BASE_C) {
    bp_penalty += AT_PENALTY;
  }
  if( seq[ thefold->pairs[ b]] != BASE_C && seq[ thefold->pairs[ j]] != BASE_C) {
    bp_penalty += AT_PENALTY;
  }  
  if( seq[ thefold->pairs[ c]] != BASE_C && seq[ thefold->pairs[ f]] != BASE_C) {
    bp_penalty += AT_PENALTY;
  }   
  
  energy = EnergyFg( i, a, d, e, thefold) + EnergyFg( b, c, f, j, thefold) +
    EnergyFz( a+1, b-1, thefold) + EnergyFz( c+1, d-1, thefold) + 
    EnergyFz( e+1, f-1, thefold) + 2*BETA_2 + bp_penalty;
  
  //printf("PkEnergy = %f\n", energy);
  return energy;
}

/* ****************************  */
DBL_TYPE EnergyFg( int i, int d, int e, int j, fold *thefold) {
  
  DBL_TYPE energy = 0.0;
  int c=-1, f=-1; //end of rightmost pair or pk
  
  DBL_TYPE multi_bp_penalty;  //bp_penalty, for the multiloop case
  int multiloop = FALSE; //extended gap
  int interiorLoop = TRUE; //regular gap
  int noPairs = TRUE; //empty gap matrix
  DBL_TYPE energyRight=0; //energy of right gap, in case it is a pk
  
  int span1=-1, span2=-1;  //pair that spans the gap
  int side = 1; //1 = between e and j, -1 = between i and d
  int stop;
  
  //  printf("Fg! %d %d %d %d\n", i,d,e,j);
  
  //check if a single pair
  if( i == d) {
    if( e == j) {
      return 0;
    }
    else {
      printf("Error in EnergyFg %d %d %d %d\n", i,d,e,j);
      exit(1);
    }
  }
  else {
    if( e == j) {
      printf("Error in EnergyFg %d %d %d %d\n", i,d,e,j);
      exit(1);
    }
  }
  
  
  
  multi_bp_penalty = 0;
  f = j - 1;  
  stop = f;
  
  while( f >= i + 1) {
    //printf("f = %d\n", f);
    if( thefold->pknots[ f] != -1) {
      noPairs = FALSE;
      interiorLoop = FALSE;
      multiloop = TRUE;
      energy += 2*ALPHA_2;
      c = thefold->pknots[ f];
      if( c <= d && side == 1) {
        printf("Error: Pseudoknot spans a gap matrix!\n");
        exit(1);
      }
      
      //printf("- %d %d\n", f+1, stop);
      
      energy +=
        EnergyPk( c, f, thefold) + BETA_1M +
        DangleEnergyWithPairs( f+1, stop, thefold) + 
        (stop - f)*ALPHA_3;
      f = c-1;
      stop = f;
    }
    else if( thefold->pairs[ f] != -1) {
      c = thefold->pairs[ f];
      if( c > f) {
        printf("Error: Unclassified pknot!\n");
        exit(1);
      }
      else if( side == 1 && c < d) {
        span1 = c; 
        span2 = f;
        noPairs = FALSE;
        side = -1;
        energy += EnergyFg( c,d,e,f, thefold);
        if( thefold->seq[c] != BASE_C && thefold->seq[f] != BASE_C) {
          multi_bp_penalty += AT_PENALTY;
        }
        //	printf("-- %d %d\n", f+1, stop);
        energyRight = DangleEnergyWithPairs( f+1, stop, thefold)
          + (stop - f)*ALPHA_3;
        f = c-1;
        stop = f;
      }
      else if( side == 1 && c > e) {
        noPairs = FALSE;
        interiorLoop = FALSE;
        multiloop = TRUE;
        //printf("c = %d, f = %d, stop = %d\n", c,f,stop);
        printf("-- %d %d\n", f+1, stop);
        
        energy += EnergyFb( c, f, thefold) + 
          DangleEnergyWithPairs( f+1, stop, thefold)
          + (stop - f)*ALPHA_3 + ALPHA_2;
        if( thefold->seq[c] != BASE_C && thefold->seq[f] != BASE_C) {
          multi_bp_penalty += AT_PENALTY;
        }
        f = c-1;
        stop = f;
      }
      else if( side == -1) {
        noPairs = FALSE;
        interiorLoop = FALSE;
        multiloop = TRUE;
        
        //printf("--- %d %d\n", f+1, stop);
        
        energy += EnergyFb( c, f, thefold) + 
          DangleEnergyWithPairs( f+1, stop, thefold)
          + (stop - f)*ALPHA_3 + ALPHA_2;
        if( thefold->seq[c] != BASE_C && thefold->seq[f] != BASE_C) {
          multi_bp_penalty += AT_PENALTY;
        }
        f = c-1;
        stop = f;
      }
      else if( c == d) {
        //printf("---- %d %d\n", f+1, stop);
        
        energyRight = 
          DangleEnergyWithPairs( f+1, stop, thefold)
          + (stop - f)*ALPHA_3;
        
        side = -1;
        f = c - 1;
        stop = f;
      }
      else {
        printf("Impossible construction in pknot!\n");
        exit(1);
      }
    }
    else {
      f--;
    }
    
  }
  
  if( noPairs == TRUE) {
    energy +=  InteriorEnergy( i, j, d, e, thefold->seq);
  }
  else if( interiorLoop == TRUE) {
    energy += InteriorEnergy( i, j, span1, span2, thefold->seq);
  }
  else if( multiloop == TRUE) {
    if( thefold->seq[ i] != BASE_C && thefold->seq[j] != BASE_C) {
      multi_bp_penalty += AT_PENALTY;
    }
    //printf("---!\n");
    
    energy += ALPHA_1 + multi_bp_penalty + 2*ALPHA_2 + 
      (stop-i)*ALPHA_3 + 
      energyRight + DangleEnergyWithPairs( i+1, c-1, thefold);
  }
  else {
    printf("Error in Fg!\n");
    exit(1);
  }
  
  return energy;
}

/* ******************************************** */
DBL_TYPE EnergyFz( int start, int stop, fold *thefold) {
  
  DBL_TYPE energy = 0.0;
  int d; //Left end of rightmost pair or pk
  
  DBL_TYPE bp_penalty;
  int j; //right end of pair or pk
  
  //printf("Fz %d %d\n", start, stop);
  
  j = stop; 
  while( j >= start) {
    if( thefold->pknots[ j] != -1) {
      d = thefold->pknots[ j];
      
      energy +=
        EnergyPk( d, j, thefold) + BETA_1P + 
        DangleEnergyWithPairs( j+1, stop, thefold) +
        BETA_3*( stop - j) + 2*BETA_2;
      
      j = d-1;
      stop = j;
    }
    else if( thefold->pairs[ j] != -1) {
      d = thefold->pairs[ j];
      if( d > j) {
        printf("Error: Unclassified pknot!\n");
        exit(1);
      }
      bp_penalty = 0;
      if( thefold->seq[d] != BASE_C && thefold->seq[j] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }
      
      energy +=
        EnergyFb( d, j, thefold) + 
        DangleEnergyWithPairs( j+1, stop, thefold) +
        bp_penalty + BETA_3*( stop - j) + BETA_2;
      
      j = d-1;
      stop = j;
    }
    else {
      j--;
    }
  }
  
  energy += DangleEnergyWithPairs( start, stop, thefold) +
    BETA_3*(stop - start + 1);
  
  return energy;
}

/* *************************************** */
#ifdef COAXIAL
DBL_TYPE minCoaxialStacking( fold *thefold, int nPairs, 
                            int *pairs, DBL_TYPE nonStackDangle) {
  /*This function will determine the minimum coaxially stacking
  configuration for a multiloop
  pairs should be of the format
  for all even a:
  pairs[a] is paired to pairs[a+1] with pairs[a] < pairs[a+1]
  pairs[0]-pairs[1] is closing pair of multiloop, and remaining pairs
  pairs[2]-pairs[3] is the rightmost (3') pair in the multiloop, and
  successive pairs move from right to left
  
  nPairs is the number of pairs (including the closing one
  in the multiloop.
  */
  
  char *seq = thefold->seq;
  DBL_TYPE minEnergy = nonStackDangle;
  DBL_TYPE energy;
  
  if( nPairs <= 2) {
    printf("Error!  Non multiloop sent to coaxially stacking subroutine!");
    exit(1);
  }
  
  //first check if closing pair can stack with 5' pair
  if( pairs[0] == pairs[nPairs*2-2] - 1) {
    
    energy = CoaxEnergy( seq[pairs[0]], seq[pairs[1]],
                        seq[pairs[nPairs*2-2]], seq[pairs[nPairs*2-1]]) +
      minCoax( 2, TRUE, 5, nPairs, pairs, thefold);
    
    minEnergy = MIN( minEnergy, energy);
    //printf("+1 %f %f\n", energy, minEnergy);
  }
  
  //next check if closing pair can stack with 3' pair
  if( pairs[1] == pairs[ 3] + 1) {
    
    energy = CoaxEnergy( seq[pairs[0]], seq[pairs[1]], 
                        seq[pairs[2]], seq[pairs[3]] ) + 
      minCoax( 3, TRUE, 3, nPairs, pairs, thefold);
    
    minEnergy = MIN( minEnergy, energy);
    //printf("+2 %f %f\n", energy, minEnergy);
  }
  
  //next check the case where closing pair is not stacked
  energy = minCoax(2, FALSE, 0, nPairs, pairs, thefold);
  minEnergy = MIN( minEnergy, energy);
  //printf("+3 %f %f\n", energy, minEnergy);
  
  return minEnergy;
}

/* ******************************** */
DBL_TYPE minCoax( int startPair, int isPreviousPairStacked, 
                 int closingPairState,
                 int nPairs, int *pairs, fold *thefold) {
                   /* This function will determine whether or not the startPair should be
                   stacked with the pair on its 5' side (find min energy).
                   
                   startPair is the current pair number being considered for stacking
                   with a pair on its left.
                   
                   isPreviousPairStacked indicates whether the previous pair is stacked 
                   or not.
                   
                   closingPairState = 0 if not stacked
                   3 if stacked on the 3' side
                   5 if stacked on the 5' side
                   */
                   
                   
                   DBL_TYPE minEnergy = NAD_INFINITY;
                   DBL_TYPE energy;
                   int i, j, h, m;
                   int whichDangle; //5 = 5' end only; 3 = 3' only, 53 = both  
                   char *seq = thefold->seq;
                   
                   //printf("%d %d %d %d\n", startPair, isPreviousPairStacked, 
                   //	 closingPairState, nPairs);
                   
                   if( startPair == nPairs) {
                     h = pairs[startPair*2 - 1];
                     if( closingPairState == 5) {
                       energy = 0;
                       if( isPreviousPairStacked == FALSE) {
                         energy = CoaxDangle( 3, h+1, pairs[ startPair*2-4] - 1,
                                             thefold->pairs, thefold->seq, 
                                             thefold->seqlength);
                       }
                       //printf("- %d %d %f\n", h, pairs[ startPair*2-4], energy);
                     }
                     else {
                       whichDangle = 53;
                       if( isPreviousPairStacked == TRUE) {
                         whichDangle = 5;
                       }
                       
                       energy = CoaxDangle( whichDangle, h+1, pairs[ startPair*2-4] - 1,
                                           thefold->pairs, thefold->seq, 
                                           thefold->seqlength);
                       //printf("energy= %d %d %d %f\n", whichDangle, h+1, 
                       //pairs[ startPair*2-3] - 1, energy);
                       
                       whichDangle = 53;
                       if( closingPairState == 3) {
                         whichDangle = 3;
                       }
                       energy += CoaxDangle( whichDangle, pairs[0]+1, pairs[ startPair*2-2] - 1,
                                            thefold->pairs, thefold->seq, 
                                            thefold->seqlength);
                       
                     }
                     minEnergy = MIN( energy, minEnergy);
                     
                   }
                   else {
                     if( pairs[startPair*2-2] == pairs[ startPair*2+1] + 1 &&
                        (startPair != nPairs - 1 || closingPairState != 5) ) {
                          //consider stacked case
                          i = pairs[ startPair*2];
                          j = pairs[ startPair*2+1];
                          h = pairs[ startPair*2-2]; 
                          m = pairs[ startPair*2-1];
                          
                          energy = CoaxEnergy( seq[m], seq[h], seq[i], seq[j]);
                          
                          if( isPreviousPairStacked == FALSE && startPair != 2) {
                            energy += CoaxDangle( 3, m+1, pairs[ startPair*2-4] - 1,
                                                 thefold->pairs, thefold->seq, 
                                                 thefold->seqlength);
                          }
                          else if( isPreviousPairStacked == FALSE) {
                            energy += CoaxDangle( 3, m+1, pairs[1] - 1,
                                                 thefold->pairs, thefold->seq, 
                                                 thefold->seqlength);
                          }
                          //if previous pair is stacked, there is no dangle
                          
                          if( startPair != nPairs - 1) {
                            energy += minCoax( startPair+2, TRUE, closingPairState, nPairs,
                                              pairs, thefold);
                          }
                          else if( closingPairState == 0) {
                            energy += CoaxDangle( 5, pairs[0] + 1, pairs[ nPairs*2-2] - 1,
                                                 thefold->pairs, thefold->seq, 
                                                 thefold->seqlength);
                          }
                          //closingPairState == 3 then no dangle
                          
                          minEnergy = MIN( minEnergy, energy);
                          //printf("test\n");
                        }
                     
                     energy = 0;
                     whichDangle = 53;
                     if( isPreviousPairStacked == TRUE) {
                       whichDangle = 5;
                     }
                     
                     h = pairs[ startPair*2 -1]; //+1
                     //consider unstacked case
                     if( startPair != 2) {
                       energy = CoaxDangle( whichDangle, h+1, pairs[ startPair*2-4] - 1,//-2
                                           thefold->pairs, thefold->seq, 
                                           thefold->seqlength);
                     }
                     else {
                       energy = CoaxDangle( whichDangle, h+1, pairs[1] - 1,
                                           thefold->pairs, thefold->seq, 
                                           thefold->seqlength);
                     }
                     
                     energy += minCoax( startPair+1, FALSE, closingPairState, nPairs,
                                       pairs, thefold);
                     
                     minEnergy = MIN( minEnergy, energy);
                   }
                   
                   //printf("return %f\n", energy);
                   return minEnergy;
                 }

/* ***************************** */

DBL_TYPE CoaxEnergy( char i, char j, char h, char m) {
  DBL_TYPE energy;
  DBL_TYPE bp_bonus = 0;
	
  //remove any AT_PENALTY if bases are stacked
  if( i != BASE_C && j != BASE_C) {
    bp_bonus -= AT_PENALTY;
  }
  if( h != BASE_C && m != BASE_C) {
    bp_bonus -= AT_PENALTY;
  }
	
  energy = HelixEnergy(i,j,h,m) + bp_bonus;  
  
  return energy;
}

/* ******************************** */
#endif
