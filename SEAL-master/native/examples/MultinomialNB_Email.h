#include <cmath>
#include <iostream>
#include <iostream>
#include <fstream>          // for ifstream
#include <cctype>           // for tolower(), isalpha()
#include <string>
#include <vector>
using namespace std;

#ifndef _MultinomialNB_Email_H
#define _MultinomialNB_Email_H

struct fileNameNode
{
	string filename;
	fileNameNode * next;
	fileNameNode(string s, fileNameNode * ptr)
	{
		filename = s;
		next = ptr;
	}
	fileNameNode()
	{
		filename = "";
		next = nullptr;
	}
};

class MultinomialNB_Email
{
	public:
		MultinomialNB_Email();
		MultinomialNB_Email(const fileNameNode & hamFilesVector, const fileNameNode &spamFilesVector, int mValue, long long int &timeElapsed);
		MultinomialNB_Email(string path, int mValue, long long int &timeElapsed);
		void train(long long int &timeElapsed);
		void classify(string path, long long int &timeElapsed);

		void getWordsVector(vector<string> &wordsGlobalVectorOutput);

		void getPrivateData(vector<string> &wordsGlobalVectorOutput, 
			vector<long long int> &hamWordsInDocumentsCountOutput, 
			vector<long long int> &spamWordsInDocumentsCountOutput, 
			int & hamMailsOutput, int & spamMailsOutput);

		void getExtendedPrivateData(const vector<string> &wordsGlobalAggregated,
			vector<long long int> &hamWordsInDocumentsCountGlobal,
			vector<long long int> &spamWordsInDocumentsCountGlobal,
			int & hamMailsOutput, int & spamMailsOutput);

		void getSelectedFeaturesParameters(vector<string> &selectedFeaturesOutput,
			vector<long double>& logOfProbWordIsHamOutput,
			vector<long double>& logOfProbWordIsSpamOutput,
			long double &logProbHamOutput, long double & logProbSpamOutput);

		void getSelectedFeaturesOnly(vector<string> &selectedFeaturesOutput);

		void printConfionMat();

		void getTVaccordingToSelectedFeatures(vector<string> selectedFeatures,
			vector<int64_t> &trainingVector);

	
	private:
		int m, hamMails = 0, spamMails = 0;
		vector<string> wordsGlobalVector;
		vector<string> selectedFeatures;

		vector<long long int> spamWordsInDocumentsCount;
		vector<long long int> hamWordsInDocumentsCount;

		vector<long long int> spamWordsFreq;
		vector<long long int> hamWordsFreq;

		long double logProbSpam, logProbHam;
		vector<long double> logOfProbWordIsHam, logOfProbWordIsSpam;		

		vector<vector<int>> confusionMatrix{ {0,0},{0,0} };
};

#endif