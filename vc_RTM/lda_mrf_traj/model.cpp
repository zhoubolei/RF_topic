#include "stdafx.h"
#include <stdlib.h>
#include "math.h"
#include "model.h"
#include <vector>


using namespace std;

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)

int randmult(vector<double> pi, int length) {
  double sum = 0.0, mass;
  int i, cc = 0;

  for (i = 0; i < length; i++)
    sum += pi[i];
  mass = myrand() * sum;
  while (cc < length-1) {
    mass -= pi[cc];
    if ( mass <= 0.00000000001 ) break;
    cc ++;
  }
  return cc;
}

Trajectory::Trajectory(size_t topicNum, size_t sourceNum, size_t sinkNum, size_t len){

	Point p;

	p.x = 0; p.y = 0; p.t = 0;
	d_points.resize(len, p);
    source_gibbs.resize(len,0);
    sink_gibbs.resize(len,0);

	d_sources.resize(len,0);
	d_sinks.resize(len,0);
	d_topicNum = topicNum;
	d_sourceNum= sourceNum;
	d_sinkNum=sinkNum;
   
}

size_t Trajectory::length() { return d_points.size(); }

Point Trajectory::startPoint() { return d_points[0]; }

Point Trajectory::endPoint() { return d_points[d_points.size()-1]; }

Point Trajectory::getPoint(size_t i){return d_points[i];}

Velocity Trajectory::startVelocity(size_t sizeOfBatch){

	size_t len;
	Velocity v;

	len = d_points.size();
	if (sizeOfBatch > len){
		v.vx = double(d_points[len-1].x - d_points[0].x)/double(d_points[len-1].t - d_points[0].t);
		v.vy = double(d_points[len-1].y - d_points[0].y)/double(d_points[len-1].t - d_points[0].t);
	}else
	{
		v.vx = double(d_points[sizeOfBatch-1].x - d_points[0].x)/double(d_points[sizeOfBatch-1].t-d_points[0].t);
        v.vy = double(d_points[sizeOfBatch-1].y - d_points[0].y)/double(d_points[sizeOfBatch-1].t-d_points[0].t);
	}

	return v;
}

Velocity Trajectory::endVelocity(size_t sizeOfBatch){

	size_t len;
	Velocity v;

	len = d_points.size();
	if (sizeOfBatch > len){
		v.vx = double(d_points[len-1].x - d_points[0].x)/double(d_points[len-1].t - d_points[0].t);
		v.vy = double(d_points[len-1].y - d_points[0].y)/double(d_points[len-1].t - d_points[0].t);
	}else
	{
		v.vx = double(d_points[len-1].x - d_points[len-sizeOfBatch].x)/double(d_points[len-1].t - d_points[len-sizeOfBatch].t);
		v.vy = double(d_points[len-1].y - d_points[len-sizeOfBatch].y)/double(d_points[len-1].t - d_points[len-sizeOfBatch].t);
	}

	return v;
}

size_t Trajectory::word(size_t i){ return d_words[i]; }

size_t Trajectory::source(size_t i){ return d_sources[i]; }

size_t Trajectory::sink(size_t i) { return d_sinks[i]; }

size_t Trajectory::topic(size_t i){ return d_topics[i]; }

vector<size_t> Trajectory::neighborIds(){ return d_neighborIds; }

vector<double> Trajectory::neighborWeights() {return d_neighborWeights; }

void Trajectory::setPoint(size_t ii, Point p, size_t s1, size_t s2){	
	d_points[ii] = p; 
	d_sources[ii]=s1; 
	d_sinks[ii]=s2;
	
	if (d_sources[ii]>0){
		source_gibbs[ii]=0;
	}
	else
	{
		source_gibbs[ii]=1;
	}
	if (d_sinks[ii]>0){
		sink_gibbs[ii]=0;
	}
	else
	{
		sink_gibbs[ii]=1;
	}




}
void Trajectory::assignTopic(size_t ii, size_t newTopic){

	size_t oldTopic;

	oldTopic = d_topics[ii];
	d_topics[ii] = newTopic;
	d_topicsStat[oldTopic]--;
    d_topicsStat[newTopic]++;
}

void Trajectory::assignSource(size_t ii, size_t newSource){

	size_t oldSource;

	oldSource = d_sources[ii];
	d_sources[ii] = newSource;
	//d_sourceStat[oldTopic]--;
 //   d_topicsStat[newTopic]++;
}

void Trajectory::assignSink(size_t ii, size_t newSink){

	size_t oldSink;

	oldSink = d_sinks[ii];
	d_sinks[ii] = newSink;
	//d_sourceStat[oldTopic]--;
 //   d_topicsStat[newTopic]++;
}

void Trajectory::initializeTopicLabels(){

	size_t i, topic;  

	d_topics.resize(d_points.size(), 0);
	d_topicsStat.resize(d_topicNum, 0);



	for (i = 0; i < d_points.size();i++){
		topic = floor(myrand() * d_topicNum);
        d_topics[i] = topic;
        d_topicsStat[topic]++;


		if (source_gibbs[i]==1){
			d_sources[i]=floor(myrand() * d_sourceNum)+1;}
		if (sink_gibbs[i]==1){
			d_sinks[i]=floor(myrand() * d_sinkNum)+1;}



	}
}

size_t Trajectory::topicStat(size_t k){ return d_topicsStat[k]; }


void Trajectory::convert2Words(size_t imW, size_t imH, size_t sizeOfCell){

	size_t i, x, y, W, H,deltaTime=10,direction,len;
	Velocity v;

	W = imW/sizeOfCell;
	H = imH/sizeOfCell;

	d_words.resize(d_points.size(), 0);	
	len = d_points.size();

	
	for (i = 0; i < d_points.size(); i++){
		//printf(",%d",i);
		
		if (deltaTime>len){v.vx=0;v.vy=0;}
		else
		{
			if (i+deltaTime>=len-1){
				v.vx=double(d_points[len-1].x-d_points[len-1-deltaTime].x)/double(d_points[len-1].t-d_points[len-1-deltaTime].t);
				v.vy=double(d_points[len-1].y-d_points[len-1-deltaTime].y)/double(d_points[len-1].t-d_points[len-1-deltaTime].t);			
			}
			else{
				v.vx=double(d_points[i+deltaTime].x - d_points[i].x)/double( d_points[i+deltaTime].t-d_points[i].t);
				v.vy=double(d_points[i+deltaTime].y - d_points[i].y)/double( d_points[i+deltaTime].t-d_points[i].t);
			}
		}


		v.vx=0.7071*v.vx-0.7071*v.vy;
		v.vy=0.7071*v.vx+0.7071*v.vy;
		if (v.vx>=0 && v.vy>=0){direction=0;}
		if (v.vx<0 && v.vy<0){direction=1;}
		if (v.vx<0 && v.vy>0){direction=2;}
		if (v.vx>0 && v.vy<0){direction=3;}
		x = (d_points[i].x - 1)/sizeOfCell;
		y = (d_points[i].y - 1)/sizeOfCell;
		d_words[i] = x * H + y + direction * W * H;

	}
	

	


}

void Trajectory::addNeighbor(size_t id, double weight){

	d_neighborIds.push_back(id);
	d_neighborWeights.push_back(weight);
}

MRF_Traj::MRF_Traj(double beta, double sigma, double eta, double kappa, double lambda, size_t topicNum, size_t sourceNum, size_t sinkNum, size_t imW, size_t imH, size_t sizeOfCell, size_t spatialRadius, size_t temporalRadius, double angleThreshold, size_t sizeOfBatch,size_t level_tree){

    d_beta = beta;
	d_sigma = sigma;
    d_lambda = lambda;
	d_eta=eta;
	d_kappa=kappa;
	d_sourceNum=sourceNum;
	d_sinkNum=sinkNum;
	d_topicNum = topicNum;
	d_imW = imW;
	d_imH = imH;
	d_termNum = (imW / sizeOfCell) * (imH /sizeOfCell) * 4;
	d_sizeOfCell = sizeOfCell;
    d_spatialRadius = spatialRadius;
    d_temporalRadius = temporalRadius;
	d_angleThreshold = angleThreshold;
    d_sizeOfBatch = sizeOfBatch;
	d_level_tree=level_tree;


    


}

size_t MRF_Traj::classqq(size_t topic, size_t word){ return d_classqq[topic * d_termNum +word]; }

size_t MRF_Traj::classSource(size_t topic, size_t word){ return d_classSource[topic * d_sourceNum +word]; }

size_t MRF_Traj::classSink(size_t topic, size_t word){ return d_classSink[topic * d_sinkNum +word]; }

void MRF_Traj::increaseClassqq(size_t topic, size_t word){ d_classqq[topic * d_termNum +word]++; }

void MRF_Traj::decreaseClassqq(size_t topic, size_t word){ d_classqq[topic * d_termNum +word]--;}

void MRF_Traj::increaseClassSource(size_t topic, size_t word){ d_classSource[topic * d_sourceNum +word]++; }

void MRF_Traj::decreaseClassSource(size_t topic, size_t word){ d_classSource[topic * d_sourceNum +word]--;}

void MRF_Traj::increaseClassSink(size_t topic, size_t word){ d_classSink[topic * d_sinkNum +word]++; }

void MRF_Traj::decreaseClassSink(size_t topic, size_t word){ d_classSink[topic * d_sinkNum +word]--;}

void MRF_Traj::buildGraph(){

	size_t i, j, minj, maxj, range = 5000;
    Point p0, p1;
	Velocity v0, v1;
	double correlation, weight;
    minj=1;maxj=d_trajSet.size();
	printf("build MRF graph...\n");
	for (i = 0; i < d_trajSet.size(); i++){

		printf("current Trk No. %d...\n",i);
		

		/*if (i < range){ minj = 0; } else { minj = i - range; }
		if (i > d_trajSet.size() - range - 1){ maxj = d_trajSet.size()-1;} else {maxj = i + range; }*/
		if (d_trajSet[i].trkSource>0 && d_trajSet[i].trkSink==0){
			
			p1 = d_trajSet[i].endPoint();
			v1 = d_trajSet[i].endVelocity(d_sizeOfBatch);
	      
			for (j = minj; j < maxj; j++){
				if (j == i) {continue;}
				p0 = d_trajSet[j].startPoint();
				if ((p0.t <= p1.t) | (p0.t > (p1.t + d_temporalRadius))) {continue;}
				if ((abs(p0.x - p1.x) + abs(p0.y - p1.y)) > d_spatialRadius) {continue;}
				v0 = d_trajSet[j].startVelocity(d_sizeOfBatch);
				correlation = v0.vx * v1.vx + v0.vy * v1.vy;
				if (correlation != 0) {correlation = correlation/(sqrt(v0.vx*v0.vx+v0.vy*v0.vy)*sqrt(v1.vx*v1.vx+v1.vy*v1.vy));}
				if (correlation < d_angleThreshold) {continue;}
				//weight = 0;
				weight = exp((correlation-1)/d_sigma);
				d_trajSet[i].d_neighborIdsSink.push_back(j);
	            d_trajSet[i].d_neighborWeightsSink.push_back(weight);
				
			}
		}

		if (d_trajSet[i].trkSource==0 && d_trajSet[i].trkSink>0){
			
			p1 = d_trajSet[i].startPoint();
			v1 = d_trajSet[i].startVelocity(d_sizeOfBatch);
	      
			for (j = minj; j < maxj; j++){
				if (j == i) {continue;}
				p0 = d_trajSet[j].endPoint();
				if ((p0.t >= p1.t) | (p1.t > (p0.t + d_temporalRadius))) {continue;}
				if ((abs(p0.x - p1.x) + abs(p0.y - p1.y)) > d_spatialRadius) {continue;}
				v0 = d_trajSet[j].endVelocity(d_sizeOfBatch);
				correlation = v0.vx * v1.vx + v0.vy * v1.vy;
				if (correlation != 0) {correlation = correlation/(sqrt(v0.vx*v0.vx+v0.vy*v0.vy)*sqrt(v1.vx*v1.vx+v1.vy*v1.vy));}
				if (correlation < d_angleThreshold) {continue;}
				//weight = 0;
				weight = exp((correlation-1)/d_sigma);
				d_trajSet[i].d_neighborIdsSource.push_back(j);
	            d_trajSet[i].d_neighborWeightsSource.push_back(weight);				
			}
		}

		if (d_trajSet[i].trkSource==0 && d_trajSet[i].trkSink==0){
			
			p1 = d_trajSet[i].endPoint();
			v1 = d_trajSet[i].endVelocity(d_sizeOfBatch);
	      
			for (j = minj; j < maxj; j++){
				if (j == i) {continue;}
				p0 = d_trajSet[j].startPoint();
				if ((p0.t <= p1.t) | (p0.t > (p1.t + d_temporalRadius))) {continue;}
				if ((abs(p0.x - p1.x) + abs(p0.y - p1.y)) > d_spatialRadius) {continue;}
				v0 = d_trajSet[j].startVelocity(d_sizeOfBatch);
				correlation = v0.vx * v1.vx + v0.vy * v1.vy;
				if (correlation != 0) {correlation = correlation/(sqrt(v0.vx*v0.vx+v0.vy*v0.vy)*sqrt(v1.vx*v1.vx+v1.vy*v1.vy));}
				if (correlation < d_angleThreshold) {continue;}
				//weight = 0;
				weight = exp((correlation-1)/d_sigma);

				d_trajSet[i].d_neighborIdsSink.push_back(j);
	            d_trajSet[i].d_neighborWeightsSink.push_back(weight);
				
				
			}

			p1 = d_trajSet[i].startPoint();
			v1 = d_trajSet[i].startVelocity(d_sizeOfBatch);
	      
			for (j = minj; j < maxj; j++){
				if (j == i) {continue;}
				p0 = d_trajSet[j].endPoint();
				if ((p0.t >= p1.t) | (p1.t > (p0.t + d_temporalRadius))) {continue;}
				if ((abs(p0.x - p1.x) + abs(p0.y - p1.y)) > d_spatialRadius) {continue;}
				v0 = d_trajSet[j].endVelocity(d_sizeOfBatch);
				correlation = v0.vx * v1.vx + v0.vy * v1.vy;
				if (correlation != 0) {correlation = correlation/(sqrt(v0.vx*v0.vx+v0.vy*v0.vy)*sqrt(v1.vx*v1.vx+v1.vy*v1.vy));}
				if (correlation < d_angleThreshold) {continue;}
				//weight = 0;
				weight = exp((correlation-1)/d_sigma);
				d_trajSet[i].d_neighborIdsSource.push_back(j);
	            d_trajSet[i].d_neighborWeightsSource.push_back(weight);				
			}
		}






	}
}

void MRF_Traj::addTree(std::vector<size_t> curPath){
	treeIndexSet.push_back(curPath);
}

void MRF_Traj::buildForest(){

	
	int nGenSink=0,nGenSource=0,nSelf=0,nGenSinkSource=0;
	int addIndex;
	size_t i,k,w,sumTree=0,topic,word,sink,source;



	printf("build forest of spanning tree with source and sink...\n");
	for (i = 0; i < d_trajSet.size(); i++){
		printf("current trajectory is %d...\n",i);
		//d_trajSet[i].d_index=0;

		curPath.clear();
		curWeight.clear();
		
		addIndex=0;
		treeIndexSet.clear();
		treeWeightSet.clear();
		
		if (d_trajSet[i].trkSink>0 && d_trajSet[i].trkSource>0){
			nSelf=nSelf+1;
			d_trajSet[i].d_neighborIds.clear();		
            d_trajSet[i].d_neighborWeights.clear();
			

			
		}
		if (d_trajSet[i].trkSource>0 && d_trajSet[i].trkSink==0){
            nGenSink=nGenSink+genSink(i);
					
		}

		if (d_trajSet[i].trkSink>0 && d_trajSet[i].trkSource==0){
            nGenSource=nGenSource+genSource(i);		
		
		}
		
		if (d_trajSet[i].trkSink==0 && d_trajSet[i].trkSource==0){
            nGenSinkSource=nGenSinkSource+genSinkSource(i);		
		
		}


		if (d_trajSet[i].d_neighborIds.size()>0){
			d_trajSet[i].d_index=1;
			for (k = 0; k < d_trajSet[i].d_neighborIds.size(); k++){
				d_trajSet[d_trajSet[i].d_neighborIds[k]].d_index=1;
			}	
		}
	}
	size_t numInclude=d_trajSet.size();
	printf("exclude the outliers ... \n");

	for (i = 0; i < d_trajSet.size(); i++){
		

		if (d_trajSet[i].d_index==0){
			printf("current trak in outliers is %d... \n",i);

			for (k = 0; k < d_trajSet[i].length(); k++){
			source=d_trajSet[i].source(k);
            sink=d_trajSet[i].sink(k);
			word = d_trajSet[i].word(k);
			topic = d_trajSet[i].topic(k);			
			d_classqq[topic * d_termNum +word]--;
			d_classSource[topic * d_sourceNum +source-1]--;
			d_classSink[topic * d_sinkNum +sink-1]--;			
			d_topicTotal[topic]--;			
			d_sourceTotal[topic]--;
			d_sinkTotal[topic]--;
			

			}
			
			numInclude--;
		}
	}


    FILE* fileptr;
	fileptr = fopen("parameter.txt","w");
	fprintf(fileptr,"nGenSink=%d,nGenSource=%d,nSelf=%d,nGenSinkSource=%d,sumTree=%d,numTrks=%d,numInclude=%d...\n,",nGenSink,nGenSource,nSelf,nGenSinkSource,nGenSinkSource+nGenSink+nGenSource+nSelf,d_trajSet.size(),numInclude);
	fclose(fileptr);
	//system("pause");


	
    printf("save neigborhood ... \n");
	fileptr = fopen("neigborhood.txt","w");
	for (k = 0; k <d_trajSet.size(); k++){
		if (d_trajSet[k].d_index==1){
		fprintf(fileptr, "%d", k);
		for (w = 0; w < d_trajSet[k].d_neighborIds.size(); w++){
			
			fprintf(fileptr, "(%d,%f)", d_trajSet[k].d_neighborIds[w], d_trajSet[k].d_neighborWeights[w]);
		}
		fprintf(fileptr, "\n");
		}
	}
	fclose(fileptr);


}

int MRF_Traj::genSink(size_t i){
	treeIndexSet.clear();
	treeWeightSet.clear();

	
	seekSink(i,0);
	double treeWeight,maxwei=0,curwei;
	size_t k,j,treeIn=-1;

	//printf("begin choose the optimal tree,tree size=%d .",treeIndexSet.size());
	
	d_trajSet[i].d_neighborIds.clear();
	d_trajSet[i].d_neighborWeights.clear();

	if (treeIndexSet.size()>0){


		for (k=0;k<treeIndexSet.size();k++){
			curwei=0;
			for (j=1;j<treeIndexSet[k].size();j++){

				curwei=curwei+treeWeightSet[k][j];
			}
			curwei=curwei/(treeWeightSet[k].size()-1);
			if (curwei>maxwei){maxwei=curwei;treeIn=k;}
		}

		for (k=0;k<treeIndexSet[treeIn].size();k++){
			d_trajSet[i].addNeighbor(treeIndexSet[treeIn][k],treeWeightSet[treeIn][k]);
		}
		return 1;

	}
	return 0;
	//printf("choose %d\n",treeIn);system("pause");


}

int MRF_Traj::genSource(size_t i){
	treeIndexSet.clear();
	treeWeightSet.clear();
	
	seekSource(i,0);

	//printf("begin choose the optimal tree,tree size=%d...",treeIndexSet.size());
	
	double treeWeight,maxwei=-1,curwei;
	size_t k,j,treeIn=-1;
	d_trajSet[i].d_neighborIds.clear();
	d_trajSet[i].d_neighborWeights.clear();
	if (treeIndexSet.size()>0){

		for (k=0;k<treeIndexSet.size();k++){
			curwei=0;
			for (j=1;j<treeIndexSet[k].size();j++){

				curwei=curwei+treeWeightSet[k][j];
			}
			curwei=curwei/(treeWeightSet[k].size()-1);
			if (curwei>maxwei){maxwei=curwei;treeIn=k;}
		}

		
		for (k=0;k<treeIndexSet[treeIn].size();k++){
			d_trajSet[i].addNeighbor(treeIndexSet[treeIn][k],treeWeightSet[treeIn][k]);
		}
	
	return 1;
	}
	//printf("choose %d\n",treeIn);system("pause");
	return 0;




}

int MRF_Traj::genSinkSource(size_t i){
	vector<size_t> sinkNeigbor,sourceNeigbor;
	vector<double> sinkWeight,sourceWeight;
	int sinkI,sourceI,ii;

	sinkI=genSink(i);
	
	if (sinkI!=0){
		sinkNeigbor=d_trajSet[i].d_neighborIds;
		sinkWeight=d_trajSet[i].d_neighborWeights;
		
	}

	sourceI=genSource(i);

	if (sourceI!=0){
		sourceNeigbor=d_trajSet[i].d_neighborIds;
		sourceWeight=d_trajSet[i].d_neighborWeights;
		
	}
	d_trajSet[i].d_neighborIds.clear();
	d_trajSet[i].d_neighborWeights.clear();

	if (sourceI!=0 && sinkI!=0){//||

		for (ii=0;ii<sinkNeigbor.size();ii++){

				 d_trajSet[i].d_neighborIds.push_back(sinkNeigbor[ii]);
				 d_trajSet[i].d_neighborWeights.push_back(sinkWeight[ii]);
		}

		
		for (ii=0;ii<sourceNeigbor.size();ii++){

				 d_trajSet[i].d_neighborIds.push_back(sourceNeigbor[ii]);
				 d_trajSet[i].d_neighborWeights.push_back(sourceWeight[ii]);
		}
		
		return 1;
	}
	


	

	return 0;
}

void MRF_Traj::seekSink(size_t i,double wei){

	vector<size_t> neighborIds;
	vector<double> neighborWeights;
	size_t k;
	curPath.push_back(i);
	curWeight.push_back(wei);
	neighborIds = d_trajSet[i].d_neighborIdsSink;
	neighborWeights = d_trajSet[i].d_neighborWeightsSink;


	if (curPath.size()>d_level_tree){
		return;
	}
	if (neighborIds.size()>0){
		for (k=0;k<neighborIds.size();k++){
			
			seekSink(neighborIds[k],neighborWeights[k]);
			
		}
	}
	if (d_trajSet[i].trkSink>0 && d_trajSet[i].trkSink!=d_trajSet[curPath[0]].trkSource){
		treeIndexSet.push_back(curPath);
		treeWeightSet.push_back(curWeight);
		curPath.pop_back();
		curWeight.pop_back();
		return;
	}
	if (neighborIds.size()==0){
		curPath.pop_back();
		curWeight.pop_back();
		return;
	}
	curPath.pop_back();
	curWeight.pop_back();
	return;
}

void MRF_Traj::seekSource(size_t i,double wei){
	
	vector<size_t> neighborIds;
    vector<double> neighborWeights;
	size_t k;
	curPath.push_back(i);
	curWeight.push_back(wei);
	neighborIds = d_trajSet[i].d_neighborIdsSource;
	neighborWeights = d_trajSet[i].d_neighborWeightsSource;


	if (curPath.size()>d_level_tree){
		return;
	}
	if (neighborIds.size()>0){
		for (k=0;k<neighborIds.size();k++){
			
			seekSource(neighborIds[k],neighborWeights[k]);
			
		}
	}
	if (d_trajSet[i].trkSource>0 && d_trajSet[i].trkSource!=d_trajSet[curPath[0]].trkSink){
		treeIndexSet.push_back(curPath);
		treeWeightSet.push_back(curWeight);
		curPath.pop_back();
		curWeight.pop_back();
		return;
	}
	if (neighborIds.size()==0){
		curPath.pop_back();
		curWeight.pop_back();
		return;
	}
	curPath.pop_back();
	curWeight.pop_back();
	return;

}




void MRF_Traj::sampleOneTraj(size_t i){

	vector<size_t> neighborIds;
	vector<double> neighborWeights;
	vector<double> topicStat, weights;
    size_t len, j, k, m, oldTopic, newTopic, word,source, sink;
    double maxn, sumWeight;

	neighborIds = d_trajSet[i].neighborIds();
	neighborWeights = d_trajSet[i].neighborWeights();
	topicStat.resize(d_topicNum, 0);
	weights.resize(d_topicNum, 0);



	//for (k = 0; k < d_topicNum; k++){
	//	topicStat[k] = d_trajSet[i].topicStat(k);
	//}

	for (m = 0; m < neighborIds.size(); m++){
		for (k = 0; k < d_topicNum; k++){
			topicStat[k] += neighborWeights[m] * d_trajSet[neighborIds[m]].topicStat(k);
		}
	}

	len = d_trajSet[i].length();
	for (j = 0; j < len; j++){
		word = d_trajSet[i].word(j);
		source=d_trajSet[i].source(j);
		sink=d_trajSet[i].sink(j);
		oldTopic = d_trajSet[i].topic(j);
		topicStat[oldTopic]--;
		maxn = 0;
		for (k = 0; k < d_topicNum; k++){
			if (topicStat[k] > maxn) {maxn = topicStat[k];}
		}
        sumWeight = 0;
		for (k = 0; k < d_topicNum; k++){
			if (k == oldTopic){
				weights[k] = exp(d_lambda * (topicStat[k] - maxn)) * double(classqq(k, word)-1+d_beta)/double(d_topicTotal[k]-1+d_beta*d_termNum);
                weights[k]=weights[k] * double(classSource(k, source-1)-1+d_eta)/double(d_sourceTotal[k]-1+d_eta*d_sourceNum);
				weights[k]=weights[k] * double(classSink(k, sink-1)-1+d_kappa)/double(d_sinkTotal[k]-1+d_kappa*d_sinkNum);
				weights[k]=weights[k] * (d_trajSet[i].topicStat(k)+0.5*d_topicNum);
			}else
			{
				weights[k] = exp(d_lambda * (topicStat[k] - maxn)) * double(classqq(k, word)+d_beta)/double(d_topicTotal[k]+d_beta*d_termNum);
			    weights[k]=weights[k] * double(classSource(k, source-1)-1+d_eta)/double(d_sourceTotal[k]+d_eta*d_sourceNum);
				weights[k]=weights[k] * double(classSink(k, sink-1)-1+d_kappa)/double(d_sinkTotal[k]+d_kappa*d_sinkNum);
				weights[k]=weights[k] * (d_trajSet[i].topicStat(k)+0.5*d_topicNum);
			}

			sumWeight += weights[k];
		}


		for (k = 0; k < d_topicNum;k++){
			weights[k] = weights[k]/sumWeight;
            //printf("%f  ", weights[k]);
		}
        //printf("\n");
        
		newTopic = randmult(weights, d_topicNum);
        //printf("newTopic= %d  \n", newTopic);
        //system( "pause ");
		if (newTopic != oldTopic){
			d_trajSet[i].assignTopic(j, newTopic);
            increaseClassqq(newTopic, word);
			decreaseClassqq(oldTopic, word);
			d_topicTotal[oldTopic]--;
			d_topicTotal[newTopic]++;

            d_sourceTotal[oldTopic]--;
			d_sourceTotal[newTopic]++;
			d_sinkTotal[oldTopic]--;
			d_sinkTotal[newTopic]++;


            increaseClassSink(newTopic, sink-1);
			decreaseClassSink(oldTopic, sink-1);
			increaseClassSource(newTopic, source-1);
			decreaseClassSource(oldTopic, source-1);


		}
		topicStat[newTopic]++;
	}

}


void MRF_Traj::sampleOneTrajSS(size_t i){

    size_t len, j, k, m, newTopic, word,oldTopic,oldSource,oldSink, newSource, newSink;
	vector<double> weightSources;
	vector<double> weightSinks;
	weightSources.resize(d_sourceNum, 0);
	weightSinks.resize(d_sinkNum, 0);
    double sumWeight;

    len = d_trajSet[i].length();

	for (j = 0; j < len; j++){
		//if (i==2){
		//	printf("curSource=%d,curSink=%d,curSourGibbs=%d, curSinkGibbs=%d/n ",d_trajSet[i].source(j),d_trajSet[i].sink(j),d_trajSet[i].source_gibbs[j],d_trajSet[i].sink_gibbs[j]);
		//}
		if (d_trajSet[i].source_gibbs[j]==1)
		{
		oldSource=d_trajSet[i].source(j);
		oldTopic=d_trajSet[i].topic(j);
		sumWeight = 0;
         
		for (k = 0; k < d_sourceNum; k++){
			if (k == oldSource-1){
				
                weightSources[k]=double(classSource(oldTopic, k)-1+d_eta)/double(d_sourceTotal[oldTopic]-1+d_eta*d_sourceNum);
 

			}else
			{
				
                weightSources[k]=double(classSource(oldTopic, k)+d_eta)/double(d_sourceTotal[oldTopic]+d_eta*d_sourceNum);
				
			}
			sumWeight += weightSources[k];
		}

	    for (k = 0; k < d_sourceNum;k++){
			weightSources[k] = weightSources[k]/sumWeight;
			//if (i==2){
			//printf("%f  ", weightSources[k]);
			//}
		}
		//if (i==2){
		//printf("\n");
		//}
		newSource = randmult(weightSources, d_sourceNum);
		//if (i==2){
		//printf("newSource= %d  ", newSource);
		//}
		//system( "pause ");

        if (newSource != (oldSource-1)){
			d_trajSet[i].assignSource(j, newSource+1);
 
			//d_sourceTotal[oldSource-1]--;
			//d_sourceTotal[newSource]++;
      //      if (i==2){
		    //  printf("sourceTotal oldSource=%d,newSource=%d\n",oldSource,newSource);
		    //}
			increaseClassSource(oldTopic, newSource);
			decreaseClassSource(oldTopic, oldSource-1);
		}

  //      if (i==2){
  //      printf("sampling Source of current Word= %d is OK\n",j);
		//}
		}

		if (d_trajSet[i].sink_gibbs[j]==1)
		{
			
		oldSink=d_trajSet[i].sink(j);
		oldTopic=d_trajSet[i].topic(j);
		sumWeight = 0;
		for (k = 0; k < d_sinkNum; k++){
			if (k == oldSink-1){
				
                weightSinks[k]=double(classSink(oldTopic, k)-1+d_eta)/double(d_sinkTotal[oldTopic]-1+d_eta*d_sinkNum);
			}else
			{
                weightSinks[k]=double(classSink(oldTopic, k)+d_eta)/double(d_sinkTotal[oldTopic]+d_eta*d_sinkNum);
			}
			sumWeight += weightSinks[k];
		}

	    for (k = 0; k < d_sinkNum;k++){
			weightSinks[k] = weightSinks[k]/sumWeight;
			//printf("%f  ", weightSinks[k]);
		}
        //printf("\n");
		newSink = randmult(weightSinks, d_sinkNum);
		//printf("chosen newSink= %d\n",newSink);
        if (newSink != (oldSink-1)){
			d_trajSet[i].assignSink(j, newSink+1);
 


			increaseClassSink(oldTopic, newSink);
			decreaseClassSink(oldTopic, oldSink-1);
		}


		}


		


	}

}

void MRF_Traj::randSample(){

	size_t i,k,w;

	for (i = 0; i < d_trajSet.size(); i++){
		if (d_trajSet[i].d_index==1)
		{
			sampleOneTraj(i);
		}
        if ((i % 1000) == 0)
			printf("trajectory %d\n", i);
	}

    FILE* fileptr;
    printf("save middle-source data ... \n");
	fileptr = fopen("mid-classSource.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_sourceNum; w++){
			fprintf(fileptr, "%d	", classSource(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save middle-sink data ... \n");
	fileptr = fopen("mid-classSink.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_sourceNum; w++){
			fprintf(fileptr, "%d	", classSink(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save middle-class source and sink data ... \n");
	fileptr = fopen("middle-classqq.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_termNum; w++){
			fprintf(fileptr, "%d	", classqq(k,w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);


	for (i = 0; i < d_trajSet.size(); i++){
		if (d_trajSet[i].d_index==1)
		{
			sampleOneTrajSS(i);
		}
        if ((i % 1000) == 0)
			printf("trajectory Source and Sink %d\n", i);
	}

	
}

void MRF_Traj::loadData(){

	size_t d, i;
	int length, x, y, t, s1,s2,numTrk=0;
	Point p;
    Trajectory* traj;
    char filename[100];
	FILE *fileptr;

	printf("load data ...\n");
	for (d = 1; d < 17; d++){
		//sprintf(filename, "trks_grand_ss8.txt");
        sprintf(filename, "filttracks_grand_%d.txt", d);
        //sprintf(filename, "parkinglot_trk.txt");
		fileptr = fopen(filename, "r");
		while ((fscanf(fileptr, "%d ", &length) != EOF))
		{
           traj = new Trajectory(d_topicNum, d_sourceNum, d_sinkNum, length);
           fscanf(fileptr, "[%d,%d]", &s1, &s2);
		  /* traj->d_source=s1;
		   traj->d_sink=s2;*/
		   for (i = 0; i < length; i++){
              fscanf(fileptr, "(%d,%d,%d)", &x, &y, &t);
			  p.x = x; p.y = y; p.t = t;
			  
			  traj->setPoint(i, p, s1, s2);
			 
		   }
		   
		   fscanf(fileptr, "\n");
		   d_trajSet.push_back(*traj);
		   d_trajSet[numTrk].trkSource=s1;
		   d_trajSet[numTrk].trkSink=s2;
		   //printf("cur TrkSource is %d, curTrkSink is %d. \n", d_trajSet[numTrk].trkSource,d_trajSet[numTrk].trkSink);
		   //system( "pause ");
		   
		   numTrk++;
           delete traj;
		}
		fclose(fileptr);
	}
	printf("%d trajectories are loaded. \n", d_trajSet.size());
}

void MRF_Traj::loadTree(){

	size_t d, i;
	int length, x, y, t, s1,s2,numTrk=0,noSS=0;
	Point p;
    Trajectory* traj;
    char filename[100];
	FILE *fileptr;

	d_trajSet.clear();

	printf("load forest of trees ...\n");
	for (d = 1; d < 2; d++){
		sprintf(filename, "trk-tree.txt");
        //sprintf(filename, "trks_filter_%d.txt", d);
        //sprintf(filename, "parkinglot_trk.txt");
		fileptr = fopen(filename, "r");
		while ((fscanf(fileptr, "%d ", &length) != EOF))
		{
           traj = new Trajectory(d_topicNum, d_sourceNum, d_sinkNum, length);
           fscanf(fileptr, "[%d,%d]", &s1, &s2);

		   if (s1==0  || s2==0){
			   noSS++;
			   printf("error:tree_source and sink...%d.\n",noSS);
		   }


		   for (i = 0; i < length; i++){
              fscanf(fileptr, "(%d,%d,%d)", &x, &y, &t);
			  p.x = x; p.y = y; p.t = t;
			  
			  traj->setPoint(i, p, s1, s2);
			 
		   }
		   
		   fscanf(fileptr, "\n");
		   d_trajSet.push_back(*traj);
		   d_trajSet[numTrk].trkSource=s1;
		   d_trajSet[numTrk].trkSink=s2;
		   //printf("cur TrkSource is %d, curTrkSink is %d. \n", d_trajSet[numTrk].trkSource,d_trajSet[numTrk].trkSink);
		   //system( "pause ");
		   
		   numTrk++;
           delete traj;
		}
		fclose(fileptr);
	}
	printf("%d trajectories are loaded. \n", d_trajSet.size());
}

void MRF_Traj::initialization(){

	size_t i, k, word, topic,source,sink,w;

	printf("initialization...\n");
	for (i = 0; i < d_trajSet.size(); i++){
		d_trajSet[i].d_index=0;
		d_trajSet[i].convert2Words(d_imW, d_imH, d_sizeOfCell);
		
		d_trajSet[i].initializeTopicLabels();
        
	}
    printf("initialization_begin compute classSource and ClassSink...\n");
	d_classqq.resize( d_topicNum * d_termNum, 0);
	d_classSource.resize(d_topicNum * d_sourceNum, 0);
    d_classSink.resize(d_topicNum * d_sinkNum, 0);
	d_topicTotal.resize(d_topicNum, 0);

	d_sourceTotal.resize(d_topicNum, 0);
	d_sinkTotal.resize(d_topicNum, 0);
    
	for (i = 0; i < d_trajSet.size(); i++){
        printf("current is %d, source=%d,sink=%d \n",i,d_trajSet[i].trkSource,d_trajSet[i].trkSink);
		for (k = 0; k < d_trajSet[i].length(); k++){

			
			word = d_trajSet[i].word(k);
			topic = d_trajSet[i].topic(k);
			source=d_trajSet[i].source(k);
            sink=d_trajSet[i].sink(k);
			//printf("topic=%d,d_terNum=%d,d_sourceNum=%d,d_sinkNum=%d \n",topic,d_termNum,d_sourceNum,d_sinkNum);
			
			d_classqq[topic * d_termNum +word]++;
			d_classSource[topic * d_sourceNum +source-1]++;
			d_classSink[topic * d_sinkNum +sink-1]++;
			//

			d_topicTotal[topic]++;
			
			d_sourceTotal[topic]++;
			d_sinkTotal[topic]++;       

		}

	}

	FILE* fileptr;
	printf("save pre class SS data ... \n");
	fileptr = fopen("pre-classqq.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_termNum; w++){
			fprintf(fileptr, "%d	", classqq(k,w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

 //   FILE* fileptr;
 //   printf("save pre-source data ... \n");
	//fileptr = fopen("pre-classSource.txt","w");
	//for (k = 0; k < d_topicNum; k++){
	//	for (w = 0; w < d_sourceNum; w++){
	//		fprintf(fileptr, "%d	", classSource(k, w));
	//	}
	//	fprintf(fileptr, "\n");
	//}
	//fclose(fileptr);

	//printf("save pre-sink data ... \n");
	//fileptr = fopen("pre-classSink.txt","w");
	//for (k = 0; k < d_topicNum; k++){
	//	for (w = 0; w < d_sourceNum; w++){
	//		fprintf(fileptr, "%d	", classSink(k, w));
	//	}
	//	fprintf(fileptr, "\n");
	//}
	//fclose(fileptr);

	buildGraph();
	buildForest();
	
    
}

void MRF_Traj::initializationTree(){

	size_t i, k, word, topic,source,sink,w;

	printf("initialization...\n");
	for (i = 0; i < d_trajSet.size(); i++){
		d_trajSet[i].convert2Words(d_imW, d_imH, d_sizeOfCell);
		
		d_trajSet[i].initializeTopicLabels();
        
	}
    printf("initialization_begin compute classSource and ClassSink...\n");
	d_classqq.resize( d_topicNum * d_termNum, 0);
	d_classSource.resize(d_topicNum * d_sourceNum, 0);
    d_classSink.resize(d_topicNum * d_sinkNum, 0);
	d_topicTotal.resize(d_topicNum, 0);

	d_sourceTotal.resize(d_topicNum, 0);
	d_sinkTotal.resize(d_topicNum, 0);
    
	for (i = 0; i < d_trajSet.size(); i++){
		for (k = 0; k < d_trajSet[i].length(); k++){

			
			word = d_trajSet[i].word(k);
			topic = d_trajSet[i].topic(k);
			source=d_trajSet[i].source(k);
            sink=d_trajSet[i].sink(k);

		//	printf("topic=%d,d_terNum=%d,d_sourceNum=%d,d_sinkNum=%d \n",topic,d_termNum,d_sourceNum,d_sinkNum);
			//printf("source=%d,sink=%d \n",source,sink);

			d_classqq[topic * d_termNum +word]++;
			d_classSource[topic * d_sourceNum +source-1]++;
			d_classSink[topic * d_sinkNum +sink-1]++;
			

			d_topicTotal[topic]++;
			
			d_sourceTotal[topic]++;
			d_sinkTotal[topic]++;
		}
	}
    FILE* fileptr;
	printf("save pre class SS data ... \n");
	fileptr = fopen("pre-classqq.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_termNum; w++){
			fprintf(fileptr, "%d	", classqq(k,w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);




   
    printf("save pre-source data ... \n");
	fileptr = fopen("pre-classSource.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_sourceNum; w++){
			fprintf(fileptr, "%d	", classSource(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save pre-sink data ... \n");
	fileptr = fopen("pre-classSink.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_sourceNum; w++){
			fprintf(fileptr, "%d	", classSink(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save topic statistics for trajectory ... \n");
	fileptr = fopen("pre-trk-topic.txt","w");
	for (k = 0; k < d_trajSet.size(); k++){
        fprintf(fileptr, "%d ", k+1);

		for (w = 0; w < d_topicNum; w++){
			    fprintf(fileptr, "(%d)", d_trajSet[k].topicStat(w));
		}
		
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);


}

void MRF_Traj::iterRandSample(int numIter){

	size_t i;

	printf("Gibbs sampling ... \n");
	for (i = 0; i < numIter; i++){
		
		randSample();
		printf("iteration %d\n", i);
	}
}

void MRF_Traj::saveData(){

	FILE* fileptr;
	size_t k, w,maxn,index;

	printf("save topic data ... \n");
	fileptr = fopen("classqq.txt","w");
	for (k = 0; k < d_topicNum; k++){
		for (w = 0; w < d_termNum; w++){
			fprintf(fileptr, "%d	", classqq(k,w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save topic source data ... \n");
	fileptr = fopen("classSource.txt","w");
	for (k = 0; k < d_topicNum; k++){
		fprintf(fileptr, "(%d)", k+1);
		for (w = 0; w < d_sourceNum; w++){
			fprintf(fileptr, "%d	", classSource(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save topic sink data ... \n");
	fileptr = fopen("classSink.txt","w");
	for (k = 0; k < d_topicNum; k++){
        fprintf(fileptr, "(%d)", k+1);
		for (w = 0; w < d_sinkNum; w++){
			
			fprintf(fileptr, "%d	", classSink(k, w));
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);


	printf("save source and sink labels for trajectory ... \n");
	fileptr = fopen("trk-ss-label.txt","w");
	for (k = 0; k < d_trajSet.size(); k++){
        fprintf(fileptr, "%d ", d_trajSet[k].length());
		for (w = 0; w < d_trajSet[k].length(); w++){
			    fprintf(fileptr, "(%d,%d)", d_trajSet[k].source(w),d_trajSet[k].sink(w));
		}
		
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save included label for trajectory ... \n");
	fileptr = fopen("included-label.txt","w");
	for (k = 0; k < d_trajSet.size(); k++){
        fprintf(fileptr, "%d ", d_trajSet[k].d_index);				
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	
	printf("save topic statistics for trajectory ... \n");
	fileptr = fopen("trk-topic.txt","w");
	for (k = 0; k < d_trajSet.size(); k++){
        fprintf(fileptr, "%d ", k+1);

		for (w = 0; w < d_topicNum; w++){
			    fprintf(fileptr, "(%d)", d_trajSet[k].topicStat(w));
		}
		
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);

	printf("save trk cluster label ... \n");
	fileptr = fopen("trk-cluster.txt","w");
	for (k = 0; k < d_trajSet.size(); k++){
		maxn=0;
		if (d_trajSet[k].d_index==1){
			
			fprintf(fileptr, "%d", k+1);
			for (w = 0; w < d_topicNum; w++){
				if (d_trajSet[k].topicStat(w)>maxn){maxn=d_trajSet[k].topicStat(w);index=w;}
				 
			}
			fprintf(fileptr, "(%d)", index);
			fprintf(fileptr, "\n");
		}
	}
	fclose(fileptr);


}