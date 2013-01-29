// lda_mrf_traj.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "model.h"


int _tmain(int argc, _TCHAR* argv[])
{
	double beta = 1, sigma = 1, lambda = 0.02, angleThreshold = 0.5;//BETA=0.5,SIGMA=0.5. ETA=0.25, KAPPA=0.25
	//size_t topicNum = 25, imW = 720, imH = 480, sizeOfCell = 10, spatialRadius = 40, temporalRadius = 300, sizeOfBatch = 30, numIter = 500;
    size_t topicNum = 30, imW = 720, imH = 480, sizeOfCell = 10, spatialRadius = 40, temporalRadius = 1500, sizeOfBatch = 30, numIter = 800;
	
	double eta=1,kappa=1;
	size_t sourceNum=8,sinkNum=8,level_tree=8;

	FILE* fileptr;
	fileptr = fopen("model-parameter.txt","w");
	fprintf(fileptr,"topicNum=%d, spatialRaidius=%d, temporalRadius=%d, treeLevel=%d, itrNum=%d \n",topicNum,spatialRadius,temporalRadius,level_tree,numIter);
	fclose(fileptr);
	
	MRF_Traj *model;

	model = new MRF_Traj(beta, eta, kappa, sigma, lambda, topicNum,sourceNum, sinkNum, imW, imH, sizeOfCell, spatialRadius, temporalRadius, angleThreshold, sizeOfBatch,level_tree);
	model->loadData();
	model->initialization();
	//model->loadTree();
	//model->initializationTree();
	model->iterRandSample(numIter);
	model->saveData();
	delete model;
	return 0;
}

