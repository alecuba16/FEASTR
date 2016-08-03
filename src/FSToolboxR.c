/*
 *gcc -std=gnu99 -I/usr/include/R  -I/usr/local/include -fpic -O2 -g -pipe -Wall -fexceptions -fstack-protector --param=ssp-buffer-size=4 -mtune=atom -fasynchronous-unwind-tables -o FSToolboxR.o -I../MIToolbox FSToolboxR.c BetaGamma.c CMIM.c CondMI.c DISR.c ICAP.c JMI.c mRMR_D.c ../MIToolbox/MutualInformation.c ../MIToolbox/Entropy.c ../MIToolbox/CalculateProbability.c ../MIToolbox/ArrayOperations.c -DCOMPILE_C -lm -shared -fpic
 *gcc -std=gnu99 -shared -L/usr/local/lib -o FSToolboxR.so FSToolboxR.o -L/usr/lib/R/lib -I../MIToolbox FSToolboxR.c BetaGamma.c CMIM.c CondMI.c DISR.c ICAP.c JMI.c mRMR_D.c ../MIToolbox/MutualInformation.c ../MIToolbox/Entropy.c ../MIToolbox/CalculateProbability.c ../MIToolbox/ArrayOperations.c -DCOMPILE_C -lm -shared -fpic
 */

/*******************************************************************************
** FSToolboxR.c is the entry point for the feature selection toolbox.
** It provides a R interface to the various selection algorithms.
**
** Initial Version - 27/06/2011
** Updated         - 22/02/2014 - Moved increment of feature numbers here.
** Port to R       - 30/06/2016 - Modified from matlab interface to R interface.
**
** Author - Adam Pocock
** Code Port to R - Alejandro Blanco Martinez
**
** Part of the FEAture Selection Toolbox (FEAST), please reference
** "Conditional Likelihood Maximisation: A Unifying Framework for Information
** Theoretic Feature Selection"
** G. Brown, A. Pocock, M.-J. Zhao, M. Lujan
** Journal of Machine Learning Research (JMLR), 2012
**
** Please check www.cs.manchester.ac.uk/~gbrown/fstoolbox for updates.
**
** Copyright (c) 2010-2014, A. Pocock, G. Brown, The University of Manchester
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

#include "FSToolboxR.h"
#include "FSAlgorithms.h"
#include "ArrayOperations.h"
#include "Entropy.h"
#include <string.h>
#include <unistd.h> 
#include <fcntl.h>

#define DEBUG 0
#define STDOUT_TO_FILE 0

//Alecuba16 look up table
#define mifs 1
#define mrmr 2
#define cmim 3
#define jmi 4
#define disr 5
#define cife 6
#define icap 7
#define condred 8
#define betagamma 9
#define cmi 10
#define mim 11
#define badkey -1

typedef struct { char *key; int val; } t_symstruct;

static t_symstruct lookuptable[] = {
  { "MIFS", mifs },
  { "mRMR", mrmr },
  { "CMIM", cmim },
  { "MIM", mim },
  { "JMI", jmi },
  { "DISR", disr },
  { "JMI", jmi },
  { "CIFE", cife },
  { "ICAP", icap },
  { "CondRed", condred },
  { "BetaGamma", betagamma },
  { "CMI", cmi },
};

#define NKEYS (sizeof(lookuptable)/sizeof(t_symstruct))
int keyfromstring(char *key)
{
  int i;
  for (i=0; i < NKEYS; i++) {
    if (strcmp(lookuptable[i].key, key) == 0)
      return lookuptable[i].val;
  }
  return badkey;
}

void trasposteMatrix(double *original,double *vectorizedMatrix,int rows,int cols){
  int i=0;
  int j=0;  
  for(i=0; i<rows; i++){    /* Iterate of each row */
    for(j=0; j<cols; j++) /* In each row, go over each col element  */
      vectorizedMatrix[(j*rows)+i]=original[(i*cols)+j];
  }
}

void print_matrix(double *mat,int rows,int cols,char *title){
  int i=0;
  int j=0;
  printf("\n -------------- Matrix:%s ------------- \n",title);
  printf("rows:%d cols:%d\n",rows,cols);
  for(i=0; i<rows; i++){    /* Iterate of each row */
     for(j=0; j<cols; j++){  /* In each row, go over each col element  */
      printf("%f ",*(mat + ((int)i*cols)+j)); /* Print each row element */
     }
     printf("\n");
  }
  printf("-----------------------------------\n");
}

void print_vector(double *vector,int length,char *title){
  int i=0;
  printf("\n -------------- Vector:%s ------------- \n",title);
  printf("length:%d\n",length);
  for(i=0; i<length; i++) printf("%f ",vector[i]); /* Print each element */
  printf("\n-----------------------------------\n");
}

/******************************************************************************
** entry point for the R call
******************************************************************************/
void FSToolboxR(char **algorithm,int *k, double *vectorizedfeatureMatrix,int *featureAsColumn,int *numberOfFeatures,int *numberOfSamples,double *initialFeatures,int *noInitialFeatures, double *vectorizedClassColumn,int *numberOfTargets, double *optionalParam1,double *optionalParam2,int *outputFeatures,int *noutputFeatures,int *error_code,char **error_msg)
{
  /*************************************************************
  ** this function takes 4-6 arguments:
  ** algorithm = which algorithm to use,
  ** k = number of features to select,
  ** vectorizedfeatureMatrix[]= matrix of features, in the following format (featureAsColumn=1):
  **  f e a t u r e s
  ** s
  ** a
  ** m
  ** p
  ** l
  ** e
  ** s
  ** Or (featureAsColumn=0):
  **  s a m p l e s
  ** f
  ** e
  ** a
  ** t
  ** u
  ** r
  ** e
  ** s
  ** vectorizedClassColumn[] = targets,
  ** optionalParam1 = (path angle or beta value),
  ** optionalParam2 = (gamma value),
  ** the arguments should all be discrete integers.
  ** and has one output:
  ** outputFeatures[] of size k
  ** error_code returns -1 if error, 0 otherwise
  ** error_msg returns Error message if error, empty string otherwise
  *************************************************************/
  
  double *outputFeatures2,entropyTest,*featureMatrix;
  //featureMatrix=(double *) malloc(sizeof(double)*(*numberOfFeatures)*(*numberOfSamples));
  //classVector=(double *) malloc(sizeof(double)*(*numberOfSamples));

  int savedSTDOUT,f;
  if(DEBUG&&STDOUT_TO_FILE){
  	f = open("log.txt", O_CREAT | O_RDWR | O_TRUNC,0777);
  	savedSTDOUT = dup(STDOUT_FILENO);
	dup2(f,STDOUT_FILENO);
	close(f);
  }

  //featureMatrix=unVectorizeMatrix(vectorizedfeatureMatrix,*numberOfSamples,*numberOfFeatures);
  if(*featureAsColumn>0){
  	 if(DEBUG) print_matrix(vectorizedfeatureMatrix,*numberOfSamples,*numberOfFeatures,"vectorizedfeatureMatrix Before");
	  //memcpy(vectorizedfeatureMatrix, featureMatrix, sizeof(double)*(*numberOfSamples)*(*numberOfFeatures));
  	  featureMatrix=(double *) malloc(sizeof(double)*(*numberOfFeatures)*(*numberOfSamples));
  	  memcpy(featureMatrix,vectorizedfeatureMatrix, sizeof(double)*(*numberOfSamples)*(*numberOfFeatures));
	  trasposteMatrix(featureMatrix,vectorizedfeatureMatrix,*numberOfSamples,*numberOfFeatures);
	  free(featureMatrix);
	  //free(vectorizedfeatureMatrix);//Problems free space with R?
  }
  // }else{
  // 	memcpy(featureMatrix,vectorizedfeatureMatrix, sizeof(double)*(*numberOfSamples)*(*numberOfFeatures));
  // }

  if(DEBUG) print_matrix(vectorizedfeatureMatrix,*numberOfSamples,*numberOfFeatures,"vectorizedfeatureMatrix");
  
  //memcpy(classVector,vectorizedClassColumn, sizeof(double)*(*numberOfSamples));
  if(DEBUG) print_vector(vectorizedClassColumn,*numberOfTargets,"vectorizedClassColumn");

  /************************************************************
  ** number to function map
  ** 1 = MIFS
  ** 2 = mRMR
  ** 3 = CMIM
  ** 4 = JMI
  ** 5 = DISR
  ** 6 = CIFE
  ** 7 = ICAP
  ** 8 = CondRed
  ** 9 = BetaGamma
  ** 10 = CMI
  *************************************************************/
  if(algorithm == NULL ||*algorithm == NULL){
    *error_code=-1;
    *error_msg="No algortihm specified or null,accepted algorithms: MIFS , mRMR,CMIM, MIM, JMI, DISR, CIFE, ICAP, CondRed,BetaGamma,CMI";
    return;
  }else{
    if(!(strcmp(*algorithm,"MIFS")==0||strcmp(*algorithm,"mRMR")==0||strcmp(*algorithm,"MIM")==0||strcmp(*algorithm,"CMIM")==0||strcmp(*algorithm,"JMI")==0
           ||strcmp(*algorithm,"DISR")==0||strcmp(*algorithm,"CIFE")==0||strcmp(*algorithm,"ICAP")==0
           ||strcmp(*algorithm,"CondRed")==0||strcmp(*algorithm,"BetaGamma")==0||strcmp(*algorithm,"CMI")==0)){
           *error_code=-1;
           *error_msg="Invalid algortihm specified,accepted algorithms: MIFS , mRMR,CMIM, MIM, JMI, DISR, CIFE, ICAP, CondRed,BetaGamma,CMI";
           return;
           
    }
    
    if (*numberOfTargets != *numberOfSamples)
    {
      sprintf(*error_msg,"Number of targets must match number of samples\n");
      sprintf(*error_msg,"Number of targets = %d, Number of Samples = %d, Number of Features = %d\n",*numberOfTargets,*numberOfSamples,*numberOfFeatures);
      *error_code=-1;
      return;
    }/*if size mismatch*/
  else if ((*k < 1) || (*k > *numberOfFeatures))
  {
    if(*k < 1)
      sprintf(*error_msg,"You have requested k = %d features, which is not possible\n",*k);
    else
      sprintf(*error_msg,"You have requested k = %d features, which is not possible because the available features are %d\n",*k,*numberOfFeatures);
    *error_code=-1;
    return;
  }
  else
  {
    /*double calculateEntropy(double *dataVector, int vectorLength)*/
    entropyTest = calculateEntropy(vectorizedClassColumn,*numberOfSamples);
    if (entropyTest < 0.0000001)
    {
      sprintf(*error_msg,"The class label Y has entropy of 0, therefore all mutual informations containing Y will be 0. No feature selection is performed\n");
      *error_code=-1;
      return;
    }
    else
    {
      switch (keyfromstring(*algorithm))
      {
      case mifs: /* MIFS*/
      {
        if (optionalParam1 == NULL) *optionalParam1=1.0;
        if (optionalParam2 == NULL) *optionalParam2=0.0;
        
        /*void BetaGamma(int k, long noOfSamples, long noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures, double beta, double gamma)*/
        BetaGamma(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures,*optionalParam1,*optionalParam2);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case mrmr: /* mRMR*/
      {
        /*void mRMR_D(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures)*/
        mRMR_D(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case cmim: /* CMIM*/
      {
        /*void CMIM(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn,double *initialFeatures,int noInitialFeatures, int *outputFeatures,int *noutputFeatures)*/
        CMIM(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,initialFeatures,*noInitialFeatures,outputFeatures,noutputFeatures);

        incrementVector(outputFeatures,*noutputFeatures);
        if(DEBUG){
        printf("Selected features(%d):",*noutputFeatures);
		int j=0;
		for (j=0;j<*noutputFeatures;j++){
		    printf(" outputFeatures[%d]:%d",j,outputFeatures[j]);
        }
        printf("\n");
        }
        break;
      }
      case mim: /* MIM*/
      {
        /*void MIM(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures)*/
        MIM(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case jmi: /* JMI */
      {
        /*void JMI(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures)*/
        JMI(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case disr: /* DISR */
      {
        /*void DISR(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures,int *noutputFeatures)*/
        DISR(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case cife: /* CIFE */
      {
        /* CIFE is Beta = 1, Gamma = 1*/
        *optionalParam1 = 1.0;
        *optionalParam2 = 1.0;
        
        /*void BetaGamma(int k, long noOfSamples, long noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures, double beta, double gamma)*/
        BetaGamma(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures,*optionalParam1,*optionalParam2);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case icap: /* ICAP */
      {
        /*void ICAP(k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,targets,output);*/
        ICAP(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case condred: /* CondRed */
      {
        /* CondRed is Beta = 0, Gamma = 1*/
        *optionalParam1 = 0.0;
        *optionalParam2 = 1.0;
        
        /*void BetaGamma(int k, long noOfSamples, long noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures, double beta, double gamma)*/
        BetaGamma(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures,*optionalParam1,*optionalParam2);
        
        incrementVector(outputFeatures,*noutputFeatures);
        break;
      }
      case betagamma: /* BetaGamma */
      {
        if (optionalParam1 == NULL || optionalParam2 == NULL)
      {
        sprintf(*error_msg,"Insufficient arguments specified for Beta Gamma FS\n");
        *error_code=-1;
        return;
      }
        else
        {
          /*void BetaGamma(int k, long noOfSamples, long noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn, int *outputFeatures,int *noutputFeatures, double beta, double gamma)*/
          BetaGamma(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,outputFeatures,noutputFeatures,*optionalParam1,*optionalParam2);
          incrementVector(outputFeatures,*noutputFeatures);
        }
        break;
      }
      case cmi: /* CMI */
      {
        //    output = (double *)mxCalloc(*k,sizeof(double));
        
        //    /*void CondMI(int k, int noOfSamples, int noOfFeatures,double *vectorizedfeatureMatrix, double *vectorizedClassColumn,int noInitialFeatures,double *initialFeatures, int *outputFeatures,int *noutputFeatures)*/
        //    CondMI(*k,*numberOfSamples,*numberOfFeatures,vectorizedfeatureMatrix,vectorizedClassColumn,noInitialFeatures,initialFeatures,outputFeatures,noutputFeatures);
        
        // i = 0;
        
        // while((outputFeatures[i] != -1) && (i < *k))
        // {
        //     i++;
        // }
        
        // plhs[0] = mxCreateDoubleMatrix(i,1,mxREAL);
        // outputFeatures = (double *)mxGetPr(plhs[0]);
        
        // for (j = 0; j < i; j++)
        // {
        //     outputFeatures[j] = outputFeatures[j] + 1; /*C indexes from 0 not 1*/
        // }/*for number of selected features
        
        // mxFree(outputFeatures);
        // outputFeatures = NULL;
        break;
      }
      }/*switch on algorithm*/
    }
  }
  }
  *error_code=0;
	if(DEBUG&&STDOUT_TO_FILE){
	        dup2(savedSTDOUT,STDOUT_FILENO);
	}
	return;
}
