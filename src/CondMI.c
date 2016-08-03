/*******************************************************************************
** CondMI.c, implements the CMI criterion using a greedy forward search
**
** Initial Version - 19/08/2010
** Updated - 23/06/2011
** Updated - 22/02/2014 - Patched calloc.
** Updated - 12/03/2016 - Changed initial value of maxMI to -1.0 to prevent segfaults when I(X;Y) = 0.0 for all X.
**
** Author - Adam Pocock
** 
** Part of the FEAture Selection Toolbox (FEAST), please reference
** "Conditional Likelihood Maximisation: A Unifying Framework for Information
** Theoretic Feature Selection"
** G. Brown, A. Pocock, M.-J. Zhao, M. Lujan
** Journal of Machine Learning Research (JMLR), 2012
**
** Please check www.cs.manchester.ac.uk/~gbrown/fstoolbox for updates.
** 
** Copyright (c) 2010-2013, A. Pocock, G. Brown, The University of Manchester
** All rights reserved.
** 
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
** 
**   - Redistributions of source code must retain the above copyright notice, this 
**     list of conditions and the following disclaimer.
**   - Redistributions in binary form must reproduce the above copyright notice, 
**     this list of conditions and the following disclaimer in the documentation 
**     and/or other materials provided with the distribution.
**   - Neither the name of The University of Manchester nor the names of its 
**     contributors may be used to endorse or promote products derived from this 
**     software without specific prior written permission.
** 
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
** ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
** WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
** DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
** ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
** (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
** LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
** ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
*******************************************************************************/

#include "FSAlgorithms.h"
#include "FSToolboxR.h"

/* for memcpy */
#include <string.h>

/* MIToolbox includes */
#include "MutualInformation.h"
#include "ArrayOperations.h"

void CondMI(int k, int noOfSamples, int noOfFeatures, double *featureMatrix, double *classColumn,int noInitialFeatures,double *initialFeatures, int *outputFeatures,int *nselectedFeatures)
{
  /*holds the class MI values*/
  double *classMI = (double *)checkedCalloc(noOfFeatures,sizeof(double));
  
  char *selectedFeatures = (char *)checkedCalloc(noOfFeatures,sizeof(char));
  
  /*holds the intra feature MI values*/
  int sizeOfMatrix = k*noOfFeatures;
  double *featureMIMatrix = (double *)checkedCalloc(sizeOfMatrix,sizeof(double));
  
  /*Changed to ensure it always picks a feature*/
  double maxMI = -1.0;
  int maxMICounter = -1;
    
  double score, currentScore, totalFeatureMI;
  int currentHighestFeature;
  
  double *conditionVector = (double *) checkedCalloc(noOfSamples,sizeof(double));
  
  int arrayPosition;
  double mi, tripEntropy;
  
  int i,j,x;
    
  /*holds the first element of each sample*/
  double **feature2D = (double **) checkedCalloc(noOfFeatures,sizeof(double *));
  firstElementOfEachSample(feature2D,featureMatrix,noOfFeatures,noOfSamples);
  
  for (i = 0; i < sizeOfMatrix;i++)
  {
    featureMIMatrix[i] = -1;
  }/*for featureMIMatrix - blank to -1*/

  if(noInitialFeatures==0){
    for (i = 0; i < noOfFeatures;i++)
    {    
      /*calculate mutual info
      **double calculateMutualInformation(double *firstVector, double *secondVector, int vectorLength);
      */
      
      classMI[i] = calculateMutualInformation(feature2D[i], classColumn, noOfSamples);
      
      if (classMI[i] > maxMI)
      {
        maxMI = classMI[i];
        maxMICounter = i;
      }/*if bigger than current maximum*/
    }/*for noOfFeatures - filling classMI*/
    printf("Maximum MI before start:%f\n",maxMI);
    selectedFeatures[maxMICounter] = 1;
    outputFeatures[0] = maxMICounter;
    (*nselectedFeatures)=1;
    memcpy(conditionVector,feature2D[maxMICounter],sizeof(double)*noOfSamples);
  }else{
    for (i = 0; i < noInitialFeatures;i++)
    { 
      /* Minus one because matlab matrix position starts on 1 , not on 0 like C. */
      j=(int) initialFeatures[i]-1;
      selectedFeatures[j] = 1;
      outputFeatures[i]=j;
      (*nselectedFeatures)++;
      mergeArrays(feature2D[j],conditionVector,conditionVector,noOfSamples);
    }
    /*
    printf("arrayOutputFeatures:");
    for (i = 0; i < k;i++)
    { 
      printf("Pos matlab(%d)real(%d)",(outputFeatures[i]+1),outputFeatures[i]);
    }
    printf("\n");
    printf("Selected Features:");
    for (i = 0; i < noOfFeatures;i++)
    { 
      if(selectedFeatures[i]>0)
        printf("Pos matlab(%d)real(%d)",(i+1),i);
    }
     printf("\n");
     */
  }

  /*****************************************************************************
  ** We have populated the classMI array, and selected the highest
  ** MI feature as the first output feature
  ** Now we move into the CondMI algorithm
  *****************************************************************************/
  printf("Scores current|highest: ");
  for (i = noInitialFeatures; i < k; i++)
  {
    score = 0.0;
    currentHighestFeature = -1;
    currentScore = 0.0;
    totalFeatureMI = 0.0;
    
    for (j = 0; j < noOfFeatures; j++)
    {
      /*if we haven't selected j*/
      if (selectedFeatures[j] == 0)
      {
        currentScore = 0.0;
        totalFeatureMI = 0.0;
        
        /*double calculateConditionalMutualInformation(double *firstVector, double *targetVector, double *conditionVector, int vectorLength);*/
        currentScore = calculateConditionalMutualInformation(feature2D[j],classColumn,conditionVector,noOfSamples);
        printf("%f|%f ",currentScore,score);
			  if (currentScore > score)
			  {
				  score = currentScore;
				  currentHighestFeature = j;
			  }
			}/*if j is unselected*/
	  }/*for number of features*/
  
    outputFeatures[i] = currentHighestFeature;
    (*nselectedFeatures)++;
    if (currentHighestFeature != -1)
    {
      selectedFeatures[currentHighestFeature] = 1;
      mergeArrays(feature2D[currentHighestFeature],conditionVector,conditionVector,noOfSamples);
    }
  
  }/*for the number of features to select*/
  printf("\n");
  FREE_FUNC(classMI);
  FREE_FUNC(conditionVector);
  FREE_FUNC(feature2D);
  FREE_FUNC(featureMIMatrix);
  FREE_FUNC(selectedFeatures);
  
  classMI = NULL;
  conditionVector = NULL;
  feature2D = NULL;
  featureMIMatrix = NULL;
  selectedFeatures = NULL;

}/*CondMI(int,int,int,double[][],double[],double[])*/

