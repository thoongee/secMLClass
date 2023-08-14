#include <time.h>                // for time()
#include "NaiveBayesBreastCancer.h"
#include <cmath>
#include <iostream>
#include <iostream>
#include <fstream>          // for ifstream
#include <cctype>           
#include <string>
#include <vector>
using std::vector;
using namespace std;

NaiveBayesBreastCancer::NaiveBayesBreastCancer(string filename)
{
	vector<int> temp; //needed to temporary store a line of the dataset to be added later on to the matrix
	ifstream inputFile(filename);

	if (inputFile.fail())
	{
		cout << "could not open file " << filename << endl;
		cin.get();
		cin.ignore();
		return;
	}

	temp.empty();
	string line;
	int featureCount = 0, recordCount = 0;
	while (getline(inputFile,line))
	{
		recordCount++;
		for (int i = 0; i < line.length(); i+=2)
		{
			if ((i!=line.length()-1)&&line.at(i + 1) == '0')
			{
				temp.push_back(10);
				i++;
			}
			else
			{
				temp.push_back(atoi(&line.at(i)));
			}
		}
		dataset.push_back(temp);
		temp.clear();
	}
	
	/*
	std::cout << "Printing the dataset matrix: " << endl;
	for (int i = 0; i < dataset.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < dataset[0].size(); j++)
			std::cout << dataset[i][j] << ",";
	}*/
	
	vector<int> tmpIntVector(10, 0);
	vector<long double> tmpDoubleVector(10, 0.0);
	for (int i = 0; i < 9; i++)
	{
		countCancerF.push_back(tmpIntVector);
		countNonCancerF.push_back(tmpIntVector);

		probCancerF.push_back(tmpDoubleVector);
		probNonCancerF.push_back(tmpDoubleVector);
		
		logProbCancerF.push_back(tmpDoubleVector);
		logProbNonCancerF.push_back(tmpDoubleVector);		
	}
	
	countCancer = 0;
	countNonCancer = 0;
		
	probCancer = 0.0;
	probNonCancer = 0.0;

	logProbCancer = 0.0;
	logProbNonCancer = 0.0;

	vector<int> tmpConfusionMtrx(2, 0);
	confusionMatrix.push_back(tmpConfusionMtrx);
	confusionMatrix.push_back(tmpConfusionMtrx);

	cout << endl << "The BREAST CANCER dataset matrix has " << dataset[0].size() << " columns (features + the class) and "
		<< dataset.size() << " records" << endl;
}

NaiveBayesBreastCancer::NaiveBayesBreastCancer(vector<vector<int>> inputDataset)
{
	dataset = inputDataset;
	vector<int> tmpIntVector(10, 0);
	vector<long double> tmpDoubleVector(10, 0.0);
	for (int i = 0; i < 9; i++)
	{
		countCancerF.push_back(tmpIntVector);
		countNonCancerF.push_back(tmpIntVector);

		probCancerF.push_back(tmpDoubleVector);
		probNonCancerF.push_back(tmpDoubleVector);

		logProbCancerF.push_back(tmpDoubleVector);
		logProbNonCancerF.push_back(tmpDoubleVector);
	}

	countCancer = 0;
	countNonCancer = 0;

	probCancer = 0.0;
	probNonCancer = 0.0;

	logProbCancer = 0.0;
	logProbNonCancer = 0.0;

	vector<int> tmpConfusionMtrx(2, 0);
	confusionMatrix.push_back(tmpConfusionMtrx);
	confusionMatrix.push_back(tmpConfusionMtrx);
}

void NaiveBayesBreastCancer::findSystemCounts(int k=0)
{
	//cout << "In findSystemCounts() Beg: " << endl;
	for (int i = 0; i < dataset.size(); i++)
	{
		//std::cout << "Record Nr:" << i << endl;
		if (dataset[i][dataset[0].size() - 1] == 1) //1 FOR NON CANCER, 2 FOR CANCER
		{
			//std::cout << "Processing a cancer record"<<endl;
			countNonCancer++;
			for (int j = 0; j < dataset[0].size()-1; j++)
			{
				countNonCancerF[j] [dataset[i][j] - 1]++;
			}
		}
		else
		{
			//std::cout << "Processing a non cancer record" << endl;
			countCancer++;
			for (int j = 0; j < dataset[0].size()-1; j++)
			{
				countCancerF[j][dataset[i][j] - 1]++;
			}
		}
	}
	if (k != 0)
	{
		countNonCancer += 10*k;
		countCancer += 10*k;
		for(int i = 0; i<countCancerF.size(); i++)
			for (int  j = 0; j < countCancerF[0].size(); j++)
			{
				countCancerF[i][j]+=k;
				countNonCancerF[i][j] += k;
			}
	}

	//cout << "The number of non cancer records is: " << countNonCancer << endl;
	//cout << "The number of cancer records is: " << countCancer << endl;
	//cout << "In findSystemCounts() End: " << endl;
	
}

void NaiveBayesBreastCancer::findProbsAndLogsFromCounts(int k = 1)
{
	//cout << "In findProbsAndLogsFromCounts() Beg: " << endl;
	probCancer = (long double)countCancer / (long double) (countCancer + countNonCancer);
	probNonCancer = 1 - (long double) probCancer;
	logProbCancer = log(probCancer);
	logProbNonCancer = log(probNonCancer);

	for (int i = 0; i < probCancerF.size(); i++)
	{
		for  (int j = 0; j < probCancerF[0].size(); j++)
		{
			probCancerF[i][j] = (long double)(countCancerF[i][j] + k) / (countCancer + 10 * k);
			logProbCancerF[i][j] = log(probCancerF[i][j]);
			probNonCancerF[i][j] = (long double)(countNonCancerF[i][j] + k) / (countNonCancer + 10 * k);
			logProbNonCancerF[i][j] = log(probNonCancerF[i][j]);
		}
	}
	//cout << "In findProbsAndLogsFromCounts() End: " << endl;

}

void NaiveBayesBreastCancer::accuracyTest()
{
	//1 FOR NON CANCER, 2 FOR CANCER
	cout << "In accuracyTest() Beg: " << endl;
	int realClass, predClass;
	long double logCancerPred = 0.0, logNonCancerPred=0.0;
	for (int i = 0; i < dataset.size(); i++)
	{
		realClass = dataset[i][dataset[0].size()-1];
		
		logCancerPred = logProbCancer;
		logNonCancerPred = logProbNonCancer;

		for (int j = 0; j < dataset[0].size() - 1; j++)
		{
			logCancerPred += logProbCancerF[j][dataset[i][j] - 1];
			logNonCancerPred += logProbNonCancerF[j][dataset[i][j] - 1];
		}

		if (logNonCancerPred > logCancerPred)
			predClass = 1;
		else
			predClass = 2;

		if (realClass == 1)
		{
			if (predClass == 1)
				confusionMatrix[0][0]++;
			else
				confusionMatrix[0][1]++;
		}
		else
		{
			if (predClass == 1)
				confusionMatrix[1][0]++;
			else
				confusionMatrix[1][1]++;
		}
	}
	cout << "In accuracyTest() End: " << endl;
}

void NaiveBayesBreastCancer::getDatasetMatrix(vector<vector<int>> &datasetMatrixOutput)
{
	datasetMatrixOutput = dataset;
}

void NaiveBayesBreastCancer::printConfusionMatrix()
{
	cout << "Printing the confusion matrix:" << endl;
	cout<<"\t\tpred_NonCancer    pred Cancer"<<endl;
	cout << "pred_NonCancer\t"<< confusionMatrix[0][0]<<"\t\t  "<< confusionMatrix[0][1] << endl;
	cout<<"pred Cancer\t" << confusionMatrix[1][0] << "\t\t  " << confusionMatrix[1][1] << endl;
	std::cout << "The accuracy of the Breast Cancer Dataset for the Naive Bayes case is "
		<< (double)100*(confusionMatrix[0][0] + confusionMatrix[1][1]) / (confusionMatrix[0][0] + confusionMatrix[0][1] + confusionMatrix[1][0] + confusionMatrix[1][1])
		<< "%" << endl << endl;

	/*for (int i = 0; i < confusionMatrix.size(); i++)
	{
		cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			cout << "pred_NonCancer"<< confusionMatrix[i][j] << "     ";
		}
	}*/
}

void NaiveBayesBreastCancer::getCounts(vector<vector<int>> &countCancerFOutput, vector<vector<int>> &countNonCancerFOutput, int &countCancerOutput, int &countNonCancerOutput)
{
	countCancerFOutput = countCancerF;
	countNonCancerFOutput = countNonCancerF;
	countCancerOutput = countCancer;
	countNonCancerOutput = countNonCancer;
}

void NaiveBayesBreastCancer::getProbs(vector<vector<long double>> &probCancerFOutput, vector<vector<long double>> &probNonCancerFOutput, long double &probCancerOutput, long double &probNonCancerOutput)
{
	probCancerFOutput = probCancerF;
	probNonCancerFOutput = probNonCancerF;
	probCancerOutput = probCancer;
	probNonCancerOutput = probNonCancer;
}

void NaiveBayesBreastCancer::getLogsOfProbs(vector<vector<long double>> &logProbCancerFOutput, vector<vector<long double>> &logProbNonCancerFOutput, long double &logProbCancerOutput, long double &logProbNonCancerOutput)
{
	logProbCancerFOutput = logProbCancerF;
	logProbNonCancerFOutput = logProbNonCancerF;
	logProbCancerOutput = logProbCancer;
	logProbNonCancerOutput = logProbNonCancer;
}

