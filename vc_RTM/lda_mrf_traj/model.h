#include <vector>
#include "cokus.h"


/* Define a point on the trajectory */
typedef struct{
   
	int x;
	int y;
    size_t t;

} Point;


typedef struct{

	double vx;
	double vy;

} Velocity;

/* 
 * Trajectory class
 *
 *  Variables:
 *
 *    d_topicNum : the number of topics
 *    d_points : points along the trajectory
 *    d_words : word values of points along the trajectory
 *    d_topics : topic assignment of each point along the trajectory
 *    d_topicsStat : the number of points on the trajectory assigned to each topic
 *    d_neightborIds : the ids of neighboring trajectories
 *    d_neighborWeights : the weights of between this trajectory and neighboring trajectories
 *
 *  Functions:
 *    
 *    Trajectory(size_t topicNum, size_t len);
 *       Constructor. Cteate the trajectory with length len and set the topic number.
 *
 *    size_t length();
 *       return the length of the trajectory.     
 *
 *    Point startPoint();
 *       return the starting point of the trajectory
 *
 *    Point endPoint();
 *       return the ending point of the trajectory.
 *
 *    Velocity startVelocity(size_t sizeOfBatch);
 *       compute the velocity at the starting point of the trajectory. sizeOfBatch is the number of points used to compute the average velocity
 *
 *    Velocity endVelocity(size_t sizeOfBatch);
 *       compute the velocity at the ending point of the trajectory. sizeOfBatch is the number of points used to compute the average velocity
 *
 *    size_t word(size_t i);
 *       return the word value of the ith point on the trajectory.
 *
 *    size_t topic(size_t i);
 *      return the topic assignment of the ith point on the trajectory.
 *
 *    std::vector<size_t> neighborIds();
 *      return the ids of neighboring trajectories.
 *
 *    std::vector<double> neighborWeights();
 *      return the weights of the neighboring trajectories 
 *
 *    void setPoint(size_t ii, Point p)
 *       set the iith point the trajectory as value p
 *
 *    void assignTopic(size_t ii, size_t newTopic);
 *       assign newTopic to the iith point on the trajectory. Update the topicsStat at the same time.
 *
 *    void initializeTopicLabels();
 *       randomly initialize topic assignment
 *
 *    size_t topicStat(size_t k);
 *       return the number of points assigned to topic k
 *    
 *    void convert2Words(size_t imW, size_t imH, size_t sizeOfCell);
 *       compute the word values of points
 *
 *    void addNeighbor(size_t id, double weight);
 *       add a neighboring trajectory
 */

class Trajectory{
private:
    size_t d_topicNum;
	size_t d_sourceNum;
	size_t d_sinkNum;
	std::vector<Point> d_points;
	std::vector<size_t> d_words;
	std::vector<size_t> d_topics;
	std::vector<size_t> d_topicsStat;
    std::vector<size_t> d_sourceStat;
    std::vector<size_t> d_sinkStat;


    std::vector<size_t> d_sources;
	std::vector<size_t> d_sinks;
	


public:
	Trajectory(size_t topicNum, size_t sourceNum, size_t sinkNum, size_t len);
	size_t length();
	Point startPoint();
	Point endPoint();
	Point getPoint(size_t i);
	Velocity startVelocity(size_t sizeOfBatch);
	Velocity endVelocity(size_t sizeOfBatch);
	size_t word(size_t i);
	size_t topic(size_t i);
    size_t source(size_t i);
	size_t sink(size_t i);
	size_t d_index;
	std::vector<size_t> neighborIds();
	std::vector<double> neighborWeights();

	std::vector<size_t> d_neighborIds;
	std::vector<size_t> d_neighborIdsSink;
    std::vector<double> d_neighborWeightsSink;
	std::vector<double> d_neighborWeights;

	std::vector<size_t> d_neighborIdsSource;
    std::vector<double> d_neighborWeightsSource;

   
	
	void setPoint(size_t ii, Point p, size_t s1,size_t s2);
	void assignTopic(size_t ii, size_t newTopic);
	void assignSource(size_t ii, size_t newSource);
	void assignSink(size_t ii, size_t newSink);


	void initializeTopicLabels();
	size_t topicStat(size_t k);
	void convert2Words(size_t imW, size_t imH, size_t sizeOfCell);
	void addNeighbor(size_t id, double weight);
    std::vector<size_t> source_gibbs;
	std::vector<size_t> sink_gibbs;
	size_t trkSource;
	size_t trkSink;
};

/*
 * MRF_Traj is the class modeling the spatial and temporal MRF for clustering the observations on the trajectories
 *
 *  Variables
 *
 *    d_beta : the Dirichlet prior for the distributions of topics over words
 *    d_sigma : the parameter controlling the weight of edge between two neighboring trajectories
 *    d_lambda : the parameter converting distance to similarity between two points
 *    d_angleThreshold : if the angle between the velocities of two trajectories is larger than a threshold, they
 *       are not considered as neighboring trajectories. Here we use correlation instead of computing the angle.
 *    d_spatialRadius : if the starting point and the ending point of two trajectories are too far in space, they are not 
 *       neighboring trajectories.
 *    d_temporalRadius : if the starting point and the ending point of two trajectories are too far in time, they are not
 *       neighboring trajectories.
 *    d_topicNum   :   the number of topics;
 *    d_imH : the height of the scene
 *    d_imW : the width of the scene
 *    d_sizeOfCell : the size of cell used to quantizing words
 *    d_sizeOfBatch : the number of points used to compute the average velocity
 *    d_trajSet : the set of trajectories
 *    d_classqq : the number of points with word value w assgined to topic k.
 *    d_topicTotal : the number of points assigned to each topic.
 *
 *  Function
 *
 *    size_t classqq(size_t topic, size_t word);
 *       return the number of points with value word and being assigned to topic
 *
 *    void increaseClassqq(size_t topic, size_t word);
 *       increase the number of points with value word and being assigned to topic
 *
 *    void decreaseClassqq(size_t topic, size_t word);
 *       decrease the number of points with value word and being assigned to topic
 *
 *    void buildGraph();
 *       build the MRF graph
 *
 *    void initialization();
 *       initialization
 *
 *    void sampleOneTraj(size_t i);
 *       sample topic assignment on one trajectory.
 *
 *    void randSample();
 *       sample all the trajectories in the data set.
 *
 *    void iterRandSample(int numIter);
 *       iterate the sampling for numIter times
 *
 *    void loadData();
 *       load trajectories.
 *
 *    void saveData();
 *       save topic models
 */

class MRF_Traj{
private:
	double d_beta, d_sigma,d_lambda, d_angleThreshold, d_eta, d_kappa;
	size_t d_topicNum, d_termNum, d_imW, d_imH, d_sizeOfCell, d_spatialRadius, d_temporalRadius, d_sizeOfBatch,d_sourceNum, d_sinkNum,d_level_tree;
	std::vector<Trajectory> d_trajSet;
	std::vector<Trajectory> d_treeSet;

	std::vector<std::vector<size_t>> treeIndexSet;
	std::vector<std::vector<double>> treeWeightSet;



	std::vector<size_t> d_classqq;
	std::vector<size_t> d_classSource;
	std::vector<size_t> d_classSink;



	std::vector<size_t> d_topicTotal;
    std::vector<size_t> d_sourceTotal;
	std::vector<size_t> d_sinkTotal;


	size_t classqq(size_t topic, size_t word);
    size_t classSource(size_t topic, size_t word);
    size_t classSink(size_t topic, size_t word);
	void increaseClassqq(size_t topic, size_t word);
	void decreaseClassqq(size_t topic, size_t word);
    void increaseClassSource(size_t topic, size_t word);
	void decreaseClassSource(size_t topic, size_t word);
    void increaseClassSink(size_t topic, size_t word);
	void decreaseClassSink(size_t topic, size_t word);
    int genSink(size_t i);
	int genSource(size_t i);
	int genSinkSource(size_t i);
	void seekSink(size_t i,double wei);
	void seekSource(size_t i,double wei);
    void buildForest();
	void buildGraph();
	void sampleOneTraj(size_t i);
	void randSample();
    void sampleOneTrajSS(size_t i);

public:
    MRF_Traj(double beta, double eta, double kappa, double sigma, double lambda, size_t topicNum, size_t sourceNum,size_t sinkNum, size_t imW, size_t imH, size_t sizeOfCell, size_t spatialRadius, size_t temporalRadius, double angleThreshold, size_t sizeOfBatch,size_t level_tree);
	void loadData();
	void loadTree();
	void initialization();
	void initializationTree();
	void iterRandSample(int numIter);
	void saveData();
	void addTree(std::vector<size_t> curPath);
	//std::vector<size_t> sink_index;
	//std::vector<size_t> source_index;
	
	std::vector<size_t> curPath;
	std::vector<double> curWeight;

};