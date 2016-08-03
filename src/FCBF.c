
#include <stdio.h>
#include <stdlib.h>
int *array;

int cmp(const void *a, const void *b){
    int ia = *(int *)a;
    int ib = *(int *)b;
    return array[ia] < array[ib] ? -1 : array[ia] > array[ib];
}

double SU(double *firstVector,int firstVectorLength,double *secondVector,int secondVectorLength);

void FCBF(int k, int noOfSamples, int noOfFeatures,double *featureMatrix, double *classColumn,double *initialFeatures,int noInitialFeatures,double threshold, int *outputFeatures,int *noOfOutput)
{
//function [selectedFeatures] = FCBF(featureMatrix,classColumn,threshold)

double *classScore = (double *)calloc(noOfFeatures,sizeof(double));
double *entireColumn=(double *)calloc(noOfSamples,sizeof(double));
double *entireColumn2=(double *)calloc(noOfSamples,sizeof(double));
int *indexScore = (double *)calloc(noOfFeatures,sizeof(double));
int i,j;

for(i=0;i<noOfFeatures;i++){
    for(j=0;j<noOfSamples;j++)
    	entireColumn[j]=featureMatrix[i][j];
    classScore[i] = SU(entireColumn,noOfSamples,classColumn,noOfSamples);
}
//[classScore indexScore] = sort(classScore,1,'descend');

int size = sizeof(double)*noOfFeatures;
int index[size];//use malloc to large size array

for(i=0;i<size;i++){
    index[i] = i;
}
array = featureMatrix;
qsort(index, size, sizeof(*index), cmp);
// indexScore = indexScore(classScore > threshold);
// classScore = classScore(classScore > threshold);
printf("\n\tclassScore\tindex\n");
lengthIndexScore=0;
double *classScoreTmp = (double *)calloc(noOfFeatures,sizeof(double));
for(i=0;i<size;i++){
    printf("%d\t%d\n", classScore[index[i]], index[i]);
    if(classScore[index[i]]>threshold){
    	indexScore[lengthIndexScore]=index[i];
    	classScoreTmp[lengthIndexScore]=classScore[index[i]];
        lengthIndexScore++;
    }
}
free(classScore);
indexScore=realloc(indexScore, lengthIndexScore*sizeof(int));
classScoreTmp=realloc(classScoreTmp,lengthIndexScore*sizeof(double));


if(lengthIndexScore>0)
    int curPosition = 1;
 else
   int  curPosition = 0;
// end
double curFeature,scoreij;
while (curPosition <= lengthIndexScore){
    j = curPosition + 1;
    curFeature = indexScore[curPosition];
    while (j <= lengthIndexScore){
    	for(i=0;i<noOfSamples;i++)
    		entireColumn[j]=featureMatrix[i][curFeature];
        for(i=0;i<noOfSamples;i++)
    		entireColumn2[j]=featureMatrix[i][indexScore[j]];
        scoreij = SU(entireColumn,noOfSamples,entireColumn2,noOfSamples);
        if (scoreij > classScore[j]){
            indexScore[j] = -1;
            classScore[j] = -1;
        }else{
            j = j + 1;
        }
    }
    curPosition = curPosition + 1;
}

for(i=0;i<size;i++){
	 outputFeatures[i]=indexScore[i];
	 (*noOfOutput)++;
}
}

double SU(double *firstVector,int firstVectorLength,double *secondVector,int secondVectorLength)
{
	double hX=calculateEntropy(firstVector,firstVectorLength);
	double hY=calculateEntropy(secondVector,secondVectorLength);
	double iXY=calculateMutualInformation(firstVector,secondVector,firstVectorLength)
	return (2 * iXY) / (hX + hY);
}
