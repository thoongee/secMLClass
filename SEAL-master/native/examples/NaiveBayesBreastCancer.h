#include <cmath>
#include <iostream>
#include <iostream>
#include <fstream>          // for ifstream
#include <cctype>           // for tolower(), isalpha()
#include <string>
#include <vector>
using namespace std;

#ifndef _NAIVEBAYESBREASTCANCER_H
#define _NAIVEBAYESBREASTCANCER_H

class NaiveBayesBreastCancer
{
public:
	NaiveBayesBreastCancer(vector<vector<int>> inputDataset);
	NaiveBayesBreastCancer(string filename);
	//k is the Laplace smoothing value. By default it is 1
	void findSystemCounts(int);
	void findProbsAndLogsFromCounts(int);
	void accuracyTest();
	void getDatasetMatrix(vector<vector<int>> &datasetMatrixOutput);
	void printConfusionMatrix();
	void getCounts(vector<vector<int>> &countCancerFOutput, vector<vector<int>> &countNonCancerFOutput, int &countCancerOutput, int &countNonCancerOutput);
	void getProbs(vector<vector<long double>> &probCancerFOutput, vector<vector<long double>> &probNonCancerFOutput, long double &probCancerOutput, long double &probNonCancerOutput);
	void getLogsOfProbs(vector<vector<long double>> &logProbCancerFOutput, vector<vector<long double>> &logProbNonCancerFOutput, long double &logProbCancerOutput, long double &logProbNonCancerOutput);

private:
	vector<vector<int>> dataset; //we keep the dataset here as a matrix of ints	
	
	vector<vector<int>> countCancerF, countNonCancerF;
	vector<vector<long double>> probCancerF, probNonCancerF;
	vector<vector<long double>> logProbCancerF, logProbNonCancerF;
	
	int countCancer, countNonCancer;
	long double probCancer = 0, probNonCancer = 0;
	long double logProbCancer = 0, logProbNonCancer = 0;
	vector<vector<int>> confusionMatrix{ {0,0},{0,0} };

};

#endif