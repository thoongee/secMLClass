#include <time.h>                // for time()
#include "MultinomialNB_Email.h"
#include <cmath>
#include <iostream>
#include <iostream>
#include <fstream>          // for ifstream
#include <cctype>           
#include <string>
#include <vector>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <experimental/filesystem>
namespace fs = std::experimental;
using std::vector;
using namespace std;

MultinomialNB_Email::MultinomialNB_Email()
{
	spamMails = 0;
	hamMails = 0;
	wordsGlobalVector.push_back("");
	spamWordsInDocumentsCount.push_back(0);
	hamWordsInDocumentsCount.push_back(0);
	spamWordsFreq.push_back(0);
	hamWordsFreq.push_back(0);
}

MultinomialNB_Email::MultinomialNB_Email(const fileNameNode &hamFilesVector, const fileNameNode &spamFilesVector, int mValue, long long int &timeElapsed)
{
	m = mValue;
	
	wordsGlobalVector.push_back("");
	spamWordsInDocumentsCount.push_back(0);
	hamWordsInDocumentsCount.push_back(0);
	spamWordsFreq.push_back(0);
	hamWordsFreq.push_back(0);

	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	string filename, ANSIWord;
	spamMails = 0;
	
	time_start = chrono::high_resolution_clock::now();	
	hamMails = 0;
	fileNameNode * ptr = hamFilesVector.next;
	while (ptr)
	{
		filename = ptr->filename;
		//std::cout <<"opening file: "<< filename << std::endl;

		ifstream inputFile(filename);
		/*if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			return 0;
		}*/

		//time_start = chrono::high_resolution_clock::now();
		vector<string> tmpLocalWords;
		vector<long long int> tmpLocalWordsInDoc;
		vector<long long int> tmpLocalWordsFreq;
		bool existsInList = false;
		/*time_end = chrono::high_resolution_clock::now();
		timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		totalTimeVector += timeVector;*/

		while (inputFile >> ANSIWord)//(getline(inputFile, line))
		{
			//ToLower(ANSIWord);
			//StripPunc(ANSIWord);
			//StripWhite(ANSIWord);

			/*
				wchar_t* UnicodeTextBuffer = new wchar_t[ANSIWord.length() + 1];
				std::wmemset(UnicodeTextBuffer, 0, ANSIWord.length() + 1);
				std::mbstowcs(UnicodeTextBuffer, ANSIWord.c_str(), ANSIWord.length());
				wWord = UnicodeTextBuffer;
				StemEnglish(wWord);
				//now the variable "word" should equal "document"
				//std::wcout << L"\nDemonstrating the stemming of an ANSI string:\n";
				//std::wcout << L"(English) Original text:\t" << ANSIWord.c_str() << std::endl;
				//std::wcout << L"(English) Stemmed text:\t" << wWord.c_str() << std::endl;
				ANSIWord = std::string(wWord.begin(), wWord.end());
				//std::cout << "the stemmed ANSIWord=" << ANSIWord << endl;
				//std::system("pause");
				*/
			//wordCount++;
			//now we process the word localy and putting it's data in a local vector
			//time_start = chrono::high_resolution_clock::now();
			existsInList = false;
			for (int i = 0; i < tmpLocalWords.size(); i++)
			{
				if (ANSIWord == tmpLocalWords[i])
				{
					existsInList = true;
					//tmpLocalWordsInDoc[i] = 1;
					tmpLocalWordsFreq[i] += 1;
					break;
				}
			}
			if (!existsInList)
			{
				tmpLocalWords.push_back(ANSIWord);
				tmpLocalWordsInDoc.push_back(1);
				tmpLocalWordsFreq.push_back(1);
			}
			/*time_end = chrono::high_resolution_clock::now();
			timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
			totalTimeVector += timeVector;*/

			/*now we process the word localy and putting it's data in a local linked list
			time_start = chrono::high_resolution_clock::now();
			tmpLocalParamLL = localParamLL;
			while ((tmpLocalParamLL->next)&&(tmpLocalParamLL->next->theWord > ANSIWord))
			{
				tmpLocalParamLL = tmpLocalParamLL->next;
			}
			if ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord == ANSIWord))
			{
				tmpLocalParamLL->next->hamDocFreq++;
				tmpLocalParamLL->next->hamDocAppearances = 1;
			}
			else
			{
				tmpLocalParamLL->next = new wordInfo(ANSIWord, 0, 1, 0, 1, 0, tmpLocalParamLL->next);
			}
			time_end = chrono::high_resolution_clock::now();
			timeLL = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
			totalTimeLL += timeLL;*/
			/*existsInList = false;
			tmpLocal = localParam->next;
			while ((tmpLocal->next))
			{
				if(ANSIWord == tmp)
			}*/
		}
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "after local processing of the mail" << endl;
		//std::cout << "LocalParamLL=" << endl;
		//printingLLData(localParamLL);
		//std::cout << "systemParamLL=" << endl;
		//printingLLData(systemParamLL);
		//std::cout << "-----------------------" << endl << endl;

		//now we put the local parameters into the global one for the vector case
		//time_start = chrono::high_resolution_clock::now();
		existsInList = false;
		for (int i = 0; i < tmpLocalWords.size(); i++)
		{
			existsInList = false;
			for (int j = 0; j < wordsGlobalVector.size(); j++)
			{
				if (tmpLocalWords[i] == wordsGlobalVector[j])
				{
					existsInList = true;
					hamWordsInDocumentsCount[j]++;
					hamWordsFreq[j] += tmpLocalWordsFreq[i];
					break;
				}
			}
			if (!existsInList)
			{
				wordsGlobalVector.push_back(tmpLocalWords[i]);
				hamWordsInDocumentsCount.push_back(1);
				hamWordsFreq.push_back(tmpLocalWordsFreq[i]);

				spamWordsInDocumentsCount.push_back(0);
				spamWordsFreq.push_back(0);
			}
		}
		/*time_end = chrono::high_resolution_clock::now();
		timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		totalTimeVector += timeVector;*/

		//now we put the local parameters into the global one for the linked list case
		//time_start = chrono::high_resolution_clock::now();
		//tmpLocalParamLL = localParamLL->next;
		//while (tmpLocalParamLL)
		//{
		//	//std::cout << "Here a" << endl;
		//	tempWordInfoPtr = tmpLocalParamLL;			
		//	
		//	tmpSystemParamLL = systemParamLL;
		//	while ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord > tmpLocalParamLL->theWord))
		//	{
		//		tmpSystemParamLL = tmpSystemParamLL->next;
		//	}
		//	if ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord == tmpLocalParamLL->theWord))
		//	{
		//		tmpSystemParamLL->next->hamDocFreq+= tmpLocalParamLL->hamDocFreq;
		//		tmpSystemParamLL->next->hamDocAppearances++;
		//	}
		//	else
		//	{
		//		nrNodes++;
		//		tmpSystemParamLL->next = 
		//			new wordInfo(tmpLocalParamLL->theWord,0,1,0, tmpLocalParamLL->hamDocFreq,0, tmpSystemParamLL->next);
		//	}			
		//	//std::cout << "Deleting node with theWord=" << tmpLocalParamLL->theWord << endl;
		//	tmpLocalParamLL = tmpLocalParamLL->next;
		//	delete tempWordInfoPtr;
		//}
		//localParamLL = new wordInfo("", 0, 0, 0, 0, 0, nullptr);
		//time_end = chrono::high_resolution_clock::now();
		//timeLL = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		//totalTimeLL += timeLL;
		hamMails++;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "After putting the local processing into the global one" << endl;
		//std::cout << "localParamLL=" << endl;
		//printingLLData(localParamLL);
		//std::cout << "systemParamLL=" << endl;
		//printingLLData(systemParamLL);
		//std::cout << "hamMails=" << hamMails << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::system("pause");

		/*if ((hamMails % 100) == 0)
		{
			std::cout << "Processing the " << hamMails << "th ham record" << endl;
			std::system("pause");
		}*/
		ptr = ptr->next;
	}

	ptr = spamFilesVector.next;
	while (ptr)
	{
		filename = ptr->filename;
		//std::cout <<"opening file: "<< filename << std::endl;

		ifstream inputFile(filename);
		/*if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			return 0;
		}*/

		vector<string> tmpLocalWords;
		vector<long long int> tmpLocalWordsInDoc;
		vector<long long int> tmpLocalWordsFreq;
		bool existsInList = false;

		while (inputFile >> ANSIWord)
		{
			//ToLower(ANSIWord);
			//StripPunc(ANSIWord);
			//StripWhite(ANSIWord);

			/*
				wchar_t* UnicodeTextBuffer = new wchar_t[ANSIWord.length() + 1];
				std::wmemset(UnicodeTextBuffer, 0, ANSIWord.length() + 1);
				std::mbstowcs(UnicodeTextBuffer, ANSIWord.c_str(), ANSIWord.length());
				wWord = UnicodeTextBuffer;
				StemEnglish(wWord);
				//now the variable "word" should equal "document"
				//std::wcout << L"\nDemonstrating the stemming of an ANSI string:\n";
				//std::wcout << L"(English) Original text:\t" << ANSIWord.c_str() << std::endl;
				//std::wcout << L"(English) Stemmed text:\t" << wWord.c_str() << std::endl;
				ANSIWord = std::string(wWord.begin(), wWord.end());
				//std::cout << "the stemmed ANSIWord=" << ANSIWord << endl;
				//std::system("pause");
				*/

			//wordCount++;
			//now we put the local parameters into the global one for the vector case
			existsInList = false;
			for (int i = 0; i < tmpLocalWords.size(); i++)
			{
				if (ANSIWord == tmpLocalWords[i])
				{
					existsInList = true;
					//tmpLocalWordsInDoc[i] = 1;
					tmpLocalWordsFreq[i] += 1;
					break;
				}
			}
			if (!existsInList)
			{
				tmpLocalWords.push_back(ANSIWord);
				tmpLocalWordsInDoc.push_back(1);
				tmpLocalWordsFreq.push_back(1);
			}

			//now we process the word localy and putting it's data in a local linked list
			/*tmpLocalParamLL = localParamLL;
			while ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord > ANSIWord))
			{
				tmpLocalParamLL = tmpLocalParamLL->next;
			}
			if ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord == ANSIWord))
			{
				tmpLocalParamLL->next->spamDocFreq++;
				tmpLocalParamLL->next->spamDocAppearances = 1;
			}
			else
			{
				tmpLocalParamLL->next = new wordInfo(ANSIWord,1,0,1,0,0,tmpLocalParamLL->next);
			}*/
		}
		//now we put the local parameters into the global one for the vector case	
		existsInList = false;
		for (int i = 0; i < tmpLocalWords.size(); i++)
		{
			existsInList = false;
			for (int j = 0; j < wordsGlobalVector.size(); j++)
			{
				if (tmpLocalWords[i] == wordsGlobalVector[j])
				{
					existsInList = true;
					spamWordsInDocumentsCount[j]++;
					spamWordsFreq[j] += tmpLocalWordsFreq[i];
					break;
				}
			}
			if (!existsInList)
			{
				wordsGlobalVector.push_back(tmpLocalWords[i]);
				spamWordsInDocumentsCount.push_back(1);
				spamWordsFreq.push_back(tmpLocalWordsFreq[i]);

				hamWordsInDocumentsCount.push_back(0);
				hamWordsFreq.push_back(0);
			}
		}

		//now we put the local parameters into the global one for the linked list case
		/*tmpLocalParamLL = localParamLL->next;
		while (tmpLocalParamLL)
		{
			tempWordInfoPtr = tmpLocalParamLL;

			tmpSystemParamLL = systemParamLL;
			while ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord > tmpLocalParamLL->theWord))
			{
				tmpSystemParamLL = tmpSystemParamLL->next;
			}
			if ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord == tmpLocalParamLL->theWord))
			{
				tmpSystemParamLL->next->spamDocFreq += tmpLocalParamLL->spamDocFreq;
				tmpSystemParamLL->next->spamDocAppearances++;
			}
			else
			{
				nrNodes++;
				tmpSystemParamLL->next =
					new wordInfo(tmpLocalParamLL->theWord,1,0, tmpLocalParamLL->spamDocFreq,0,0, tmpSystemParamLL->next);
			}
			tmpLocalParamLL = tmpLocalParamLL->next;
			delete tempWordInfoPtr;
		}
		localParamLL = new wordInfo("", 0, 0, 0, 0, 0, nullptr);*/
		spamMails++;
		/*if ((spamMails % 100) == 0)
		{
			std::cout << "Processing the " << spamMails << "th spam record" << endl;
		}*/
		ptr = ptr->next;
	}
	
	time_end = chrono::high_resolution_clock::now();
	time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	timeElapsed = time_diff.count() / 1000;
	std::cout << "The initialization was done in [" << time_diff.count() / 1000 << " ms]" << endl;
	std::cout << "In total we have " << wordsGlobalVector.size() << " total words in this dataset"<< endl;
}

MultinomialNB_Email::MultinomialNB_Email(string path, int mValue, long long int &timeElapsed)
{
	m = mValue;
	/*vector<string> tmpSelectedFeatures(m + 1, "");
	selectedFeatures = tmpSelectedFeatures;
	vector<long long int> spamWordsNOTInDocumentsCount; spamWordsNOTInDocumentsCount.push_back(0);
	vector<long long int> hamWordsNOTInDocumentsCount; hamWordsNOTInDocumentsCount.push_back(0);
	vector<long double> tmpLogOfProbWordIsHam(m + 1, 0), tmpLogOfProbWordIsSpam(m + 1, 0);
	logOfProbWordIsHam = tmpLogOfProbWordIsHam; logOfProbWordIsSpam = tmpLogOfProbWordIsSpam;*/
		
	wordsGlobalVector.push_back("");
	spamWordsInDocumentsCount.push_back(0);
	hamWordsInDocumentsCount.push_back(0);
	spamWordsFreq.push_back(0);
	hamWordsFreq.push_back(0);

	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);	
	std::cout << "Doing the training now" << endl;
	
	string hamPath = path + "\\ham", filename, ANSIWord;
	time_start = chrono::high_resolution_clock::now();
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(hamPath))
	{
		filename = entry.path().string();
		//std::cout <<"opening file: "<< filename << std::endl;

		ifstream inputFile(filename);
		/*if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			return 0;
		}*/

		//time_start = chrono::high_resolution_clock::now();
		vector<string> tmpLocalWords;
		vector<long long int> tmpLocalWordsInDoc;
		vector<long long int> tmpLocalWordsFreq;
		bool existsInList = false;
		/*time_end = chrono::high_resolution_clock::now();
		timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		totalTimeVector += timeVector;*/

		while (inputFile >> ANSIWord)//(getline(inputFile, line))
		{
			//ToLower(ANSIWord);
			//StripPunc(ANSIWord);
			//StripWhite(ANSIWord);

			/*
				wchar_t* UnicodeTextBuffer = new wchar_t[ANSIWord.length() + 1];
				std::wmemset(UnicodeTextBuffer, 0, ANSIWord.length() + 1);
				std::mbstowcs(UnicodeTextBuffer, ANSIWord.c_str(), ANSIWord.length());
				wWord = UnicodeTextBuffer;
				StemEnglish(wWord);
				//now the variable "word" should equal "document"
				//std::wcout << L"\nDemonstrating the stemming of an ANSI string:\n";
				//std::wcout << L"(English) Original text:\t" << ANSIWord.c_str() << std::endl;
				//std::wcout << L"(English) Stemmed text:\t" << wWord.c_str() << std::endl;
				ANSIWord = std::string(wWord.begin(), wWord.end());
				//std::cout << "the stemmed ANSIWord=" << ANSIWord << endl;
				//std::system("pause");
				*/
			//wordCount++;
				//now we process the word localy and putting it's data in a local vector
				//time_start = chrono::high_resolution_clock::now();
			
			existsInList = false;
			for (int i = 0; i < tmpLocalWords.size(); i++)
			{
				if (ANSIWord == tmpLocalWords[i])
				{
					existsInList = true;
					//tmpLocalWordsInDoc[i] = 1;
					tmpLocalWordsFreq[i] += 1;
					break;
				}
			}
			if (!existsInList)
			{
				tmpLocalWords.push_back(ANSIWord);
				tmpLocalWordsInDoc.push_back(1);
				tmpLocalWordsFreq.push_back(1);
			}
			/*time_end = chrono::high_resolution_clock::now();
			timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
			totalTimeVector += timeVector;*/

			/*now we process the word localy and putting it's data in a local linked list
			time_start = chrono::high_resolution_clock::now();
			tmpLocalParamLL = localParamLL;
			while ((tmpLocalParamLL->next)&&(tmpLocalParamLL->next->theWord > ANSIWord))
			{
				tmpLocalParamLL = tmpLocalParamLL->next;
			}
			if ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord == ANSIWord))
			{
				tmpLocalParamLL->next->hamDocFreq++;
				tmpLocalParamLL->next->hamDocAppearances = 1;
			}
			else
			{
				tmpLocalParamLL->next = new wordInfo(ANSIWord, 0, 1, 0, 1, 0, tmpLocalParamLL->next);
			}
			time_end = chrono::high_resolution_clock::now();
			timeLL = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
			totalTimeLL += timeLL;*/
			/*existsInList = false;
			tmpLocal = localParam->next;
			while ((tmpLocal->next))
			{
				if(ANSIWord == tmp)
			}*/
		}
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "after local processing of the mail" << endl;
		//std::cout << "LocalParamLL=" << endl;
		//printingLLData(localParamLL);
		//std::cout << "systemParamLL=" << endl;
		//printingLLData(systemParamLL);
		//std::cout << "-----------------------" << endl << endl;

		//now we put the local parameters into the global one for the vector case
		//time_start = chrono::high_resolution_clock::now();
		existsInList = false;
		for (int i = 0; i < tmpLocalWords.size(); i++)
		{
			existsInList = false;
			for (int j = 0; j < wordsGlobalVector.size(); j++)
			{
				if (tmpLocalWords[i] == wordsGlobalVector[j])
				{
					existsInList = true;
					hamWordsInDocumentsCount[j]++;
					hamWordsFreq[j] += tmpLocalWordsFreq[i];
					break;
				}
			}
			if (!existsInList)
			{
				wordsGlobalVector.push_back(tmpLocalWords[i]);
				hamWordsInDocumentsCount.push_back(1);
				hamWordsFreq.push_back(tmpLocalWordsFreq[i]);

				spamWordsInDocumentsCount.push_back(0);
				spamWordsFreq.push_back(0);
			}
		}
		/*time_end = chrono::high_resolution_clock::now();
		timeVector = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		totalTimeVector += timeVector;*/

		//now we put the local parameters into the global one for the linked list case
		//time_start = chrono::high_resolution_clock::now();
		//tmpLocalParamLL = localParamLL->next;
		//while (tmpLocalParamLL)
		//{
		//	//std::cout << "Here a" << endl;
		//	tempWordInfoPtr = tmpLocalParamLL;			
		//	
		//	tmpSystemParamLL = systemParamLL;
		//	while ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord > tmpLocalParamLL->theWord))
		//	{
		//		tmpSystemParamLL = tmpSystemParamLL->next;
		//	}
		//	if ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord == tmpLocalParamLL->theWord))
		//	{
		//		tmpSystemParamLL->next->hamDocFreq+= tmpLocalParamLL->hamDocFreq;
		//		tmpSystemParamLL->next->hamDocAppearances++;
		//	}
		//	else
		//	{
		//		nrNodes++;
		//		tmpSystemParamLL->next = 
		//			new wordInfo(tmpLocalParamLL->theWord,0,1,0, tmpLocalParamLL->hamDocFreq,0, tmpSystemParamLL->next);
		//	}			
		//	//std::cout << "Deleting node with theWord=" << tmpLocalParamLL->theWord << endl;
		//	tmpLocalParamLL = tmpLocalParamLL->next;
		//	delete tempWordInfoPtr;
		//}
		//localParamLL = new wordInfo("", 0, 0, 0, 0, 0, nullptr);
		//time_end = chrono::high_resolution_clock::now();
		//timeLL = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
		//totalTimeLL += timeLL;
		hamMails++;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "After putting the local processing into the global one" << endl;
		//std::cout << "localParamLL=" << endl;
		//printingLLData(localParamLL);
		//std::cout << "systemParamLL=" << endl;
		//printingLLData(systemParamLL);
		//std::cout << "hamMails=" << hamMails << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::cout << "-----------------------" << endl << endl;
		//std::system("pause");
		
		/*if ((hamMails % 100) == 0)
		{
			std::cout << "Processing the " << hamMails << "th ham record" << endl;
			std::system("pause");
		}*/		
	}

	string spamPath = path + "\\spam";
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(spamPath))
	{
		filename = entry.path().string();
		//std::cout <<"opening file: "<< filename << std::endl;

		ifstream inputFile(filename);
		/*if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			return 0;
		}*/

		vector<string> tmpLocalWords;
		vector<long long int> tmpLocalWordsInDoc;
		vector<long long int> tmpLocalWordsFreq;
		bool existsInList = false;

		while (inputFile >> ANSIWord)
		{
			//ToLower(ANSIWord);
			//StripPunc(ANSIWord);
			//StripWhite(ANSIWord);

			/*
				wchar_t* UnicodeTextBuffer = new wchar_t[ANSIWord.length() + 1];
				std::wmemset(UnicodeTextBuffer, 0, ANSIWord.length() + 1);
				std::mbstowcs(UnicodeTextBuffer, ANSIWord.c_str(), ANSIWord.length());
				wWord = UnicodeTextBuffer;
				StemEnglish(wWord);
				//now the variable "word" should equal "document"
				//std::wcout << L"\nDemonstrating the stemming of an ANSI string:\n";
				//std::wcout << L"(English) Original text:\t" << ANSIWord.c_str() << std::endl;
				//std::wcout << L"(English) Stemmed text:\t" << wWord.c_str() << std::endl;
				ANSIWord = std::string(wWord.begin(), wWord.end());
				//std::cout << "the stemmed ANSIWord=" << ANSIWord << endl;
				//std::system("pause");
				*/

			//wordCount++;
			//now we put the local parameters into the global one for the vector case
			existsInList = false;
			for (int i = 0; i < tmpLocalWords.size(); i++)
			{
				if (ANSIWord == tmpLocalWords[i])
				{
					existsInList = true;
					//tmpLocalWordsInDoc[i] = 1;
					tmpLocalWordsFreq[i] += 1;
					break;
				}
			}
			if (!existsInList)
			{
				tmpLocalWords.push_back(ANSIWord);
				tmpLocalWordsInDoc.push_back(1);
				tmpLocalWordsFreq.push_back(1);
			}

			//now we process the word localy and putting it's data in a local linked list
			/*tmpLocalParamLL = localParamLL;
			while ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord > ANSIWord))
			{
				tmpLocalParamLL = tmpLocalParamLL->next;
			}
			if ((tmpLocalParamLL->next) && (tmpLocalParamLL->next->theWord == ANSIWord))
			{
				tmpLocalParamLL->next->spamDocFreq++;
				tmpLocalParamLL->next->spamDocAppearances = 1;
			}
			else
			{
				tmpLocalParamLL->next = new wordInfo(ANSIWord,1,0,1,0,0,tmpLocalParamLL->next);
			}*/
		}
		//now we put the local parameters into the global one for the vector case	
		existsInList = false;
		for (int i = 0; i < tmpLocalWords.size(); i++)
		{
			existsInList = false;
			for (int j = 0; j < wordsGlobalVector.size(); j++)
			{
				if (tmpLocalWords[i] == wordsGlobalVector[j])
				{
					existsInList = true;
					spamWordsInDocumentsCount[j]++;
					spamWordsFreq[j] += tmpLocalWordsFreq[i];
					break;
				}
			}
			if (!existsInList)
			{
				wordsGlobalVector.push_back(tmpLocalWords[i]);
				spamWordsInDocumentsCount.push_back(1);
				spamWordsFreq.push_back(tmpLocalWordsFreq[i]);

				hamWordsInDocumentsCount.push_back(0);
				hamWordsFreq.push_back(0);
			}
		}

		//now we put the local parameters into the global one for the linked list case
		/*tmpLocalParamLL = localParamLL->next;
		while (tmpLocalParamLL)
		{
			tempWordInfoPtr = tmpLocalParamLL;

			tmpSystemParamLL = systemParamLL;
			while ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord > tmpLocalParamLL->theWord))
			{
				tmpSystemParamLL = tmpSystemParamLL->next;
			}
			if ((tmpSystemParamLL->next) && (tmpSystemParamLL->next->theWord == tmpLocalParamLL->theWord))
			{
				tmpSystemParamLL->next->spamDocFreq += tmpLocalParamLL->spamDocFreq;
				tmpSystemParamLL->next->spamDocAppearances++;
			}
			else
			{
				nrNodes++;
				tmpSystemParamLL->next =
					new wordInfo(tmpLocalParamLL->theWord,1,0, tmpLocalParamLL->spamDocFreq,0,0, tmpSystemParamLL->next);
			}
			tmpLocalParamLL = tmpLocalParamLL->next;
			delete tempWordInfoPtr;
		}
		localParamLL = new wordInfo("", 0, 0, 0, 0, 0, nullptr);*/
		spamMails++;
		/*if ((spamMails % 100) == 0)
		{
			std::cout << "Processing the " << spamMails << "th spam record" << endl;
		}*/
	}
	
	/*
	for (int i = 1; i < wordsGlobalVector.size(); i++)
	{
		hamWordsNOTInDocumentsCount.push_back(hamMails - hamWordsInDocumentsCount[i]);
		spamWordsNOTInDocumentsCount.push_back(spamMails - spamWordsInDocumentsCount[i]);
	}
	int totalNrMails = spamMails + hamMails;
	long double probSpam = (long double)spamMails / (long double)(totalNrMails), probHam = 1 - probSpam;
	logProbHam = log(probHam);
	logProbSpam = log(probSpam);

	vector<long double> probWordInDoc(wordsGlobalVector.size(), 0),
		probWordNOTInDoc(wordsGlobalVector.size(), 0); //p(tk) and p(not(tk))
	vector<long double> probWordInHamDocs(wordsGlobalVector.size(), 0),
		probWordNOTInHamDocs(wordsGlobalVector.size(), 0),
		probWordInSpamDocs(wordsGlobalVector.size(), 0),
		probWordNOTInSpamDocs(wordsGlobalVector.size(), 0);
	vector<long double> informationGain(wordsGlobalVector.size(), -9999999999.9);

	//vector<long long double> probHamWordsNOTInDoc;
	for (int i = 1; i < wordsGlobalVector.size(); i++)
	{
		if ((hamWordsInDocumentsCount[i] + spamWordsInDocumentsCount[i]) >= 5)
		{
			probWordInDoc[i] = ((long double)hamWordsInDocumentsCount[i] + (long double)spamWordsInDocumentsCount[i]) / ((long double)totalNrMails);
			probWordNOTInDoc[i] = (long double)(1 - probWordInDoc[i]);
			//(hamWordsNOTInDocumentsCount[i] + spamWordsNOTInDocumentsCount[i]) / (spamMails + hamMails);

			probWordInHamDocs[i] = (long double)(hamWordsInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInHamDocs[i] = (long double)(hamWordsNOTInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordInSpamDocs[i] = (long double)(spamWordsInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInSpamDocs[i] = (long double)(spamWordsNOTInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());

			informationGain[i] = probWordInHamDocs[i] * log((probWordInHamDocs[i]) / (probWordInDoc[i] * probHam))
				+ probWordNOTInHamDocs[i] * log((probWordNOTInHamDocs[i]) / (probWordNOTInDoc[i] * probHam))
				+ probWordInSpamDocs[i] * log((probWordInSpamDocs[i]) / (probWordInDoc[i] * probSpam))
				+ probWordNOTInSpamDocs[i] * log((probWordNOTInSpamDocs[i]) / (probWordNOTInDoc[i] * probSpam));
		}
	}

	//long double probWordInDocLL, probWordNOTInDocLL, probWordInHamDocsLL, probWordNOTInHamDocsLL,
	//	probWordInSpamDocsLL, probWordNOTInSpamDocsLL;
	//tmpSystemParamLL = systemParamLL->next;
	//if ((tmpSystemParamLL->hamDocAppearances + tmpSystemParamLL->spamDocAppearances) > 5)
	//{
	//	probWordInDocLL = ((long double)tmpSystemParamLL->hamDocAppearances + (long double)tmpSystemParamLL->spamDocAppearances) / (long double)totalNrMails;
	//	probWordNOTInDocLL = 1 - probWordInDocLL;
	//	probWordInHamDocsLL = (long double)(tmpSystemParamLL->hamDocAppearances + 1) / (long double)(totalNrMails + nrNodes);
	//}

	string tmpString;
	long long int tmpInt;
	long double tmpDouble;
	for (int i = 1; i < wordsGlobalVector.size() - 1; i++)
	{
		for (int j = i; j < wordsGlobalVector.size() - 1; j++)
		{
			if (informationGain[i] < informationGain[j])
			{
				//swap all the stuff
				//wordsGlobalVector
				tmpString = wordsGlobalVector[i];
				wordsGlobalVector[i] = wordsGlobalVector[j];
				wordsGlobalVector[j] = tmpString;
				//spamWordsInDocumentsCount
				tmpDouble = spamWordsInDocumentsCount[i];
				spamWordsInDocumentsCount[i] = spamWordsInDocumentsCount[j];
				spamWordsInDocumentsCount[j] = tmpDouble;

				//hamWordsInDocumentsCount
				tmpDouble = hamWordsInDocumentsCount[i];
				hamWordsInDocumentsCount[i] = hamWordsInDocumentsCount[j];
				hamWordsInDocumentsCount[j] = tmpDouble;
				
				//mutual information
				tmpDouble = informationGain[i];
				informationGain[i] = informationGain[j];
				informationGain[j] = tmpDouble;
				//spamWordsFreq
				tmpInt = spamWordsFreq[i];
				spamWordsFreq[i] = spamWordsFreq[j];
				spamWordsFreq[j] = tmpInt;
				//hamWordsFreq
				tmpInt = hamWordsFreq[i];
				hamWordsFreq[i] = hamWordsFreq[j];
				hamWordsFreq[j] = tmpInt;
			}
		}
	}

	if (wordsGlobalVector.size() < m)
		m = wordsGlobalVector.size();
	

	//std::cout << "THe feautes with the top " << m << " highest mutual information are:" << endl;
	//std::cout << "Feature:\tinformationGain,\tlogOfProbWordIsHam,\tlogOfProbWordIsSpam,\t" << endl;
	for (int i = 1; i <= m; i++)
	{
		selectedFeatures[i] = wordsGlobalVector[i];
		logOfProbWordIsHam[i] = log(((long double)hamWordsFreq[i]) / (long double)(hamMails));
		logOfProbWordIsSpam[i] = log(((long double)spamWordsFreq[i]) / (long double)(spamMails));
		//std::cout << selectedFeatures[i] << ": " << informationGain[i] << ", " << logOfProbWordIsHam[i]
		//	<< ", " << logOfProbWordIsSpam[i] << endl;		
	}

	for (int i = 1; i < m; i++)
	{
		for (int j = i + 1; j <= m; j++)
		{
			if (selectedFeatures[i] < selectedFeatures[j])
			{
				tmpString = selectedFeatures[i];
				selectedFeatures[i] = selectedFeatures[j];
				selectedFeatures[j] = tmpString;

				tmpDouble = logOfProbWordIsHam[i];
				logOfProbWordIsHam[i] = logOfProbWordIsHam[j];
				logOfProbWordIsHam[j] = tmpDouble;

				tmpDouble = logOfProbWordIsSpam[i];
				logOfProbWordIsSpam[i] = logOfProbWordIsSpam[j];
				logOfProbWordIsSpam[j] = tmpDouble;
			}
		}		
	}*/

	time_end = chrono::high_resolution_clock::now();
	time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	timeElapsed = time_diff.count()/1000;
	std::cout << "In total we have " << wordsGlobalVector.size() << " total words in this dataset" << endl;
	std::cout << "The initialization was done in [" << time_diff.count() / 1000 << " ms]" << endl;
}

void MultinomialNB_Email::getPrivateData(vector<string> &wordsGlobalVectorOutput,
	vector<long long int> &hamWordsInDocumentsCountOutput,
	vector<long long int> &spamWordsInDocumentsCountOutput,	
	int & hamMailsOutput, int & spamMailsOutput)
{
	wordsGlobalVectorOutput = wordsGlobalVector;
	hamWordsInDocumentsCountOutput = hamWordsInDocumentsCount;
	spamWordsInDocumentsCountOutput = spamWordsInDocumentsCount;	
	hamMailsOutput = hamMails;
	spamMailsOutput = spamMails;
}

void MultinomialNB_Email::getWordsVector(vector<string> &wordsGlobalVectorOutput)
{
	wordsGlobalVectorOutput = wordsGlobalVector;
}


void MultinomialNB_Email::getExtendedPrivateData(const vector<string> &wordsGlobalAggregated,
	vector<long long int> &hamWordsInDocumentsCountGlobal,
	vector<long long int> &spamWordsInDocumentsCountGlobal,
	int & hamMailsOutput, int & spamMailsOutput)
{
	hamMailsOutput = hamMails;
	spamMailsOutput = spamMails;

	for (int i = 0; i < wordsGlobalAggregated.size(); i++)
	{
		for (int j = 0; j < wordsGlobalVector.size(); j++)
		{
			if (wordsGlobalAggregated[i] == wordsGlobalVector[j])
			{
				hamWordsInDocumentsCountGlobal[i] = hamWordsInDocumentsCount[j];
				spamWordsInDocumentsCountGlobal[i] = spamWordsInDocumentsCount[j];
				break;
			}			
		}
	}
}


void  MultinomialNB_Email::getSelectedFeaturesParameters(vector<string> &selectedFeaturesOutput,
	vector<long double>& logOfProbWordIsHamOutput,
	vector<long double>& logOfProbWordIsSpamOutput,
	long double &logProbHamOutput, long double & logProbSpamOutput)
{
	selectedFeaturesOutput = selectedFeatures;
	logOfProbWordIsHamOutput = logOfProbWordIsHam;
	logOfProbWordIsSpamOutput = logOfProbWordIsSpam;
	logProbHamOutput = logProbHam;
	logProbSpamOutput = logProbSpam;	
}

void MultinomialNB_Email::getSelectedFeaturesOnly(vector<string> &selectedFeaturesOutput)
{
	selectedFeaturesOutput = selectedFeatures;
}


void MultinomialNB_Email::classify(string path, long long int &timeElapsed)
{
	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

	int countQuery = 0;
	string hamPath = path + "\\ham", filename, word;
	long double logSpamQueryProb = 0.0, logHamQueryProb = 0.0;
	
	time_start = chrono::high_resolution_clock::now();
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(hamPath))
	{
		filename = entry.path().string();
		//std::cout <<"opening file: "<< filename << std::endl;
		logSpamQueryProb = logProbSpam;
		logHamQueryProb = logProbHam;

		ifstream inputFile(filename);
    	/*	if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			exit(0);
		}*/

		int begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
		while (inputFile >> word)
		{
			//wordCount++;
			//existsInList = false;
			begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
			while (!((begIndex + 1) == endIndex))
			{
				/*std::cout << "begIndex = " << begIndex << endl;
				std::cout << "midIndex = " << midIndex << endl;
				std::cout << "endIndex = " << endIndex << endl;*/

				if (word == selectedFeatures[midIndex])
				{
					logSpamQueryProb += logOfProbWordIsSpam[midIndex];
					logHamQueryProb += logOfProbWordIsHam[midIndex];
					break;
				}
				else if (word > selectedFeatures[midIndex])
				{
					endIndex = midIndex;
					midIndex = (begIndex + endIndex) / 2;
				}
				else
				{
					begIndex = midIndex;
					midIndex = (begIndex + endIndex) / 2;
				}
			}
			if (word == selectedFeatures[begIndex])
			{
				logSpamQueryProb += logOfProbWordIsSpam[begIndex];
				logHamQueryProb += logOfProbWordIsHam[begIndex];
			}
			else if (word == selectedFeatures[endIndex])
			{
				logSpamQueryProb += logOfProbWordIsSpam[endIndex];
				logHamQueryProb += logOfProbWordIsHam[endIndex];
			}
		}

		if (logHamQueryProb > logSpamQueryProb)
			confusionMatrix[0][0]++;
		else
			confusionMatrix[0][1]++;

		/*countQuery++;
		if ((countQuery % 100) == 0)
		{
			std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
			std::system("pause");
		}*/		
	}

	countQuery = 0;
	string spamPath = path + "\\spam";
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(spamPath))
	{
		filename = entry.path().string();
		//std::cout <<"opening file: "<< filename << std::endl;
		logSpamQueryProb = logProbSpam;
		logHamQueryProb = logProbHam;

		ifstream inputFile(filename);
		/*if (inputFile.fail())
		{
			cout << "could not open file " << filename << endl;
			cin.get();
			cin.ignore();
			return 0;
		}*/

		int begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
		while (inputFile >> word)
		{
			//wordCount++;
			//existsInList = false;
			begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
			while (!((begIndex + 1) == endIndex))
			{
				if (word == selectedFeatures[midIndex])
				{
					logSpamQueryProb += logOfProbWordIsSpam[midIndex];
					logHamQueryProb += logOfProbWordIsHam[midIndex];
					break;
				}
				else if (word > selectedFeatures[midIndex])
				{
					endIndex = midIndex;
					midIndex = (begIndex + endIndex) / 2;
				}
				else
				{
					begIndex = midIndex;
					midIndex = (begIndex + endIndex) / 2;
				}
			}
			if (word == selectedFeatures[begIndex])
			{
				logSpamQueryProb += logOfProbWordIsSpam[begIndex];
				logHamQueryProb += logOfProbWordIsHam[begIndex];
			}
			else if (word == selectedFeatures[endIndex])
			{
				logSpamQueryProb += logOfProbWordIsSpam[endIndex];
				logHamQueryProb += logOfProbWordIsHam[endIndex];
			}

		}

		if (logSpamQueryProb > logHamQueryProb)
			confusionMatrix[1][1]++;
		else
			confusionMatrix[1][0]++;

		/*countQuery++;
		if ((countQuery % 100) == 0)
		{
			std::cout << "Processing the " << countQuery << "th spam QUERY record" << endl;
			std::system("pause");
		}*/
		
	}
	time_end = chrono::high_resolution_clock::now();
	time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	timeElapsed = time_diff.count() / 1000;
	std::cout << "The classification was done in [" << time_diff.count() / 1000 << " ms]" << endl;
}

void MultinomialNB_Email::printConfionMat()
{
	std::cout << "Printing the confusion mat for the improved classification case:" << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl << "The accuracy of the improved classification is: " << (long double)100 * (confusionMatrix[0][0] + confusionMatrix[1][1]) / (long double)(hamMails + spamMails) << " %" << endl << endl;
}

void MultinomialNB_Email::train(long long int &timeElapsed)
{
	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		
	vector<string> tmpSelectedFeatures(m + 1, "");
	selectedFeatures = tmpSelectedFeatures;
	vector<long long int> spamWordsNOTInDocumentsCount; spamWordsNOTInDocumentsCount.push_back(0);
	vector<long long int> hamWordsNOTInDocumentsCount; hamWordsNOTInDocumentsCount.push_back(0);
	vector<long double> tmpLogOfProbWordIsHam(m + 1, 0), tmpLogOfProbWordIsSpam(m + 1, 0);
	logOfProbWordIsHam = tmpLogOfProbWordIsHam; logOfProbWordIsSpam = tmpLogOfProbWordIsSpam;


	time_start = chrono::high_resolution_clock::now();
	for (int i = 1; i < wordsGlobalVector.size(); i++)
	{
		hamWordsNOTInDocumentsCount.push_back(hamMails - hamWordsInDocumentsCount[i]);
		spamWordsNOTInDocumentsCount.push_back(spamMails - spamWordsInDocumentsCount[i]);
	}

	int totalNrMails = spamMails + hamMails;
	long double probSpam = (long double)(spamMails + m) / (long double)(totalNrMails + 2*m), probHam = 1 - probSpam;
	logProbHam = log(probHam);
	logProbSpam = log(probSpam);

	vector<long double> probWordInDoc(wordsGlobalVector.size(), 0),
		probWordNOTInDoc(wordsGlobalVector.size(), 0); //p(tk) and p(not(tk))
	vector<long double> probWordInHamDocs(wordsGlobalVector.size(), 0),
		probWordNOTInHamDocs(wordsGlobalVector.size(), 0),
		probWordInSpamDocs(wordsGlobalVector.size(), 0),
		probWordNOTInSpamDocs(wordsGlobalVector.size(), 0);
	vector<long double> informationGain(wordsGlobalVector.size(), -9999999999.9);

	int countOfLessThanFive = 0;
	//vector<long long double> probHamWordsNOTInDoc;
	for (int i = 0; i < wordsGlobalVector.size(); i++)
	{
		if ((hamWordsInDocumentsCount[i] + spamWordsInDocumentsCount[i]) >= 5)
		{
			/*probWordInDoc[i] = ((long double)hamWordsInDocumentsCount[i] + (long double)spamWordsInDocumentsCount[i]) / ((long double)totalNrMails);
			probWordNOTInDoc[i] = (long double)(1 - probWordInDoc[i]);
			//(hamWordsNOTInDocumentsCount[i] + spamWordsNOTInDocumentsCount[i]) / (spamMails + hamMails);

			probWordInHamDocs[i] = (long double)(hamWordsInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInHamDocs[i] = (long double)(hamWordsNOTInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordInSpamDocs[i] = (long double)(spamWordsInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInSpamDocs[i] = (long double)(spamWordsNOTInDocumentsCount[i] + 1) / (long double)(totalNrMails + wordsGlobalVector.size());
			
			informationGain[i] = probWordInHamDocs[i] * log((probWordInHamDocs[i]) / (probWordInDoc[i] * probHam))
				+ probWordNOTInHamDocs[i] * log((probWordNOTInHamDocs[i]) / (probWordNOTInDoc[i] * probHam))
				+ probWordInSpamDocs[i] * log((probWordInSpamDocs[i]) / (probWordInDoc[i] * probSpam))
				+ probWordNOTInSpamDocs[i] * log((probWordNOTInSpamDocs[i]) / (probWordNOTInDoc[i] * probSpam));
			*/
			
			////////////

			probWordInDoc[i] = ((long double)hamWordsInDocumentsCount[i] + (long double)spamWordsInDocumentsCount[i]); // ((long double)totalNrMails);
			probWordNOTInDoc[i] = (long double)(totalNrMails - probWordInDoc[i]);

			probWordInHamDocs[i] = (long double)(hamWordsInDocumentsCount[i] + 1); // (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInHamDocs[i] = (long double)(hamWordsNOTInDocumentsCount[i] + 1); //(long double)(totalNrMails + wordsGlobalVector.size());
			probWordInSpamDocs[i] = (long double)(spamWordsInDocumentsCount[i] + 1); // (long double)(totalNrMails + wordsGlobalVector.size());
			probWordNOTInSpamDocs[i] = (long double)(spamWordsNOTInDocumentsCount[i] + 1); // (long double)(totalNrMails + wordsGlobalVector.size());

			informationGain[i] = probWordInHamDocs[i] * log((probWordInHamDocs[i])  / (probWordInDoc[i] * probHam))
				+ probWordNOTInHamDocs[i] * log((probWordNOTInHamDocs[i]) / (probWordNOTInDoc[i] * probHam))
				+ probWordInSpamDocs[i] * log((probWordInSpamDocs[i]) / (probWordInDoc[i] * probSpam))
				+ probWordNOTInSpamDocs[i] * log((probWordNOTInSpamDocs[i]) / (probWordNOTInDoc[i] * probSpam));

		}
		else
		{
			countOfLessThanFive++;
		}
	}
	std::cout << "At the train f-ion() of MultinomialNB_Email class countOfLessThanFive=" << countOfLessThanFive << endl;

	/*long double probWordInDocLL, probWordNOTInDocLL, probWordInHamDocsLL, probWordNOTInHamDocsLL,
		probWordInSpamDocsLL, probWordNOTInSpamDocsLL;
	tmpSystemParamLL = systemParamLL->next;
	if ((tmpSystemParamLL->hamDocAppearances + tmpSystemParamLL->spamDocAppearances) > 5)
	{
		probWordInDocLL = ((long double)tmpSystemParamLL->hamDocAppearances + (long double)tmpSystemParamLL->spamDocAppearances) / (long double)totalNrMails;
		probWordNOTInDocLL = 1 - probWordInDocLL;

		probWordInHamDocsLL = (long double)(tmpSystemParamLL->hamDocAppearances + 1) / (long double)(totalNrMails + nrNodes);
	}*/

	string tmpString;
	long long int tmpInt;
	long double tmpDouble;
	for (int i = 1; i < wordsGlobalVector.size() - 1; i++)
	{
		for (int j = i; j < wordsGlobalVector.size() - 1; j++)
		{
			if (informationGain[i] < informationGain[j])
			{
				//swap all the stuff
				//wordsGlobalVector
				tmpString = wordsGlobalVector[i];
				wordsGlobalVector[i] = wordsGlobalVector[j];
				wordsGlobalVector[j] = tmpString;

				//spamWordsInDocumentsCount
				tmpDouble = spamWordsInDocumentsCount[i];
				spamWordsInDocumentsCount[i] = spamWordsInDocumentsCount[j];
				spamWordsInDocumentsCount[j] = tmpDouble;

				//hamWordsInDocumentsCount
				tmpDouble = hamWordsInDocumentsCount[i];
				hamWordsInDocumentsCount[i] = hamWordsInDocumentsCount[j];
				hamWordsInDocumentsCount[j] = tmpDouble;

				//mutual information
				tmpDouble = informationGain[i];
				informationGain[i] = informationGain[j];
				informationGain[j] = tmpDouble;
				
				//spamWordsFreq
				tmpInt = spamWordsFreq[i];
				spamWordsFreq[i] = spamWordsFreq[j];
				spamWordsFreq[j] = tmpInt;
				
				//hamWordsFreq
				tmpInt = hamWordsFreq[i];
				hamWordsFreq[i] = hamWordsFreq[j];
				hamWordsFreq[j] = tmpInt;
			}
		}
	}

	if (wordsGlobalVector.size() < m)
		m = wordsGlobalVector.size();


	//std::cout << "THe feautes with the top " << m << " highest mutual information are:" << endl;
	//std::cout << "Feature:\tinformationGain,\tlogOfProbWordIsHam,\tlogOfProbWordIsSpam,\t" << endl;
	for (int i = 1; i <= m; i++)
	{
		selectedFeatures[i] = wordsGlobalVector[i];
		//if (hamWordsFreq[i] == 0)
		//	std::cout << "hamWordsFreq[" << i << "]=" << hamWordsFreq[i] << endl;
		logOfProbWordIsHam[i] = log(((long double)hamWordsFreq[i]+1) / (long double)(hamMails + m));
		//if (spamWordsFreq[i] == 0)
		//	std::cout << "spamWordsFreq[" << i << "]=" << spamWordsFreq[i] << endl;
		logOfProbWordIsSpam[i] = log(((long double)spamWordsFreq[i]+1) / (long double)(spamMails+ m));
		//std::cout << selectedFeatures[i] << ": " << informationGain[i] << ", " << logOfProbWordIsHam[i]
		//	<< ", " << logOfProbWordIsSpam[i] << endl;		
	}

	/*
	std::cout << endl << endl << endl;
	std::system("pause");
	std::cout << "BEFORE sorting by words-inside CLASS" << endl;
	std::cout << "selectedFeature\t" << "information gain" << endl;
	for (int i = 1; i < m; i++)
	{
		std::cout << "selectedFeatures[" << i << "]=" << selectedFeatures[i] << "\t  " << informationGain[i] << endl;
	}
	std::cout << endl << endl << endl;
	std::system("pause");
	*/
	
	for (int i = 1; i < m; i++)
	{
		for (int j = i + 1; j <= m; j++)
		{
			if (selectedFeatures[i] < selectedFeatures[j])
			{
				tmpString = selectedFeatures[i];
				selectedFeatures[i] = selectedFeatures[j];
				selectedFeatures[j] = tmpString;

				tmpDouble = logOfProbWordIsHam[i];
				logOfProbWordIsHam[i] = logOfProbWordIsHam[j];
				logOfProbWordIsHam[j] = tmpDouble;

				tmpDouble = logOfProbWordIsSpam[i];
				logOfProbWordIsSpam[i] = logOfProbWordIsSpam[j];
				logOfProbWordIsSpam[j] = tmpDouble;

				// swap all the stuff
				//wordsGlobalVector
				tmpString = wordsGlobalVector[i];
				wordsGlobalVector[i] = wordsGlobalVector[j];
				wordsGlobalVector[j] = tmpString;

				//spamWordsInDocumentsCount
				tmpDouble = spamWordsInDocumentsCount[i];
				spamWordsInDocumentsCount[i] = spamWordsInDocumentsCount[j];
				spamWordsInDocumentsCount[j] = tmpDouble;

				//hamWordsInDocumentsCount
				tmpDouble = hamWordsInDocumentsCount[i];
				hamWordsInDocumentsCount[i] = hamWordsInDocumentsCount[j];
				hamWordsInDocumentsCount[j] = tmpDouble;

				//mutual information
				tmpDouble = informationGain[i];
				informationGain[i] = informationGain[j];
				informationGain[j] = tmpDouble;

				//spamWordsFreq
				tmpInt = spamWordsFreq[i];
				spamWordsFreq[i] = spamWordsFreq[j];
				spamWordsFreq[j] = tmpInt;

				//hamWordsFreq
				tmpInt = hamWordsFreq[i];
				hamWordsFreq[i] = hamWordsFreq[j];
				hamWordsFreq[j] = tmpInt;
			}
		}
		/*if ((i % 100) == 0)
		{
			std::cout << "Processing the " << i << "th selected feature" << endl;
		}*/
	}
	
	/*std::cout << endl << endl << endl;
	std::system("pause");
	std::cout << "AFTER sorting by words-inside CLASS" << endl;
	std::cout << "selectedFeature\t" << "information gain" << endl;
	for (int i = 1; i < m; i++)
	{
		std::cout << "selectedFeatures[" << i << "]=" << selectedFeatures[i] << "\t  " << informationGain[i] << endl;
	}
	std::cout << endl << endl << endl;
	std::system("pause");
	*/
	
	time_end = chrono::high_resolution_clock::now();
	time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	timeElapsed = time_diff.count() / 1000;
	std::cout << "The real trainig of the selected features was done in " << timeElapsed << " ms" << endl;
}

void MultinomialNB_Email::getTVaccordingToSelectedFeatures(vector<string> selectedFeatures,
	vector<int64_t> &trainingVector)
{
	//std::cout << "selectedFeatures.size()=" << selectedFeatures.size() << endl;
	//std::cout << "trainingVector.size()=" << trainingVector.size() << endl;
	//std::cout << "wordsGlobalVector.size()=" << wordsGlobalVector.size() << endl;
	//std::cout << "hamWordsFreq.size()=" << hamWordsFreq.size() << endl;
	//std::cout << "spamWordsFreq.size()=" << spamWordsFreq.size() << endl <<endl;

	for (int i = 0; i < selectedFeatures.size(); i++)
	{
		for (int j = 0; j < wordsGlobalVector.size(); j++)
		{
			if (selectedFeatures[i] == wordsGlobalVector[j])
			{
				trainingVector[i + 1] = (int64_t) hamWordsFreq[j];
				trainingVector[i + 1 + m + 1] = (int64_t)spamWordsFreq[j];
			}
		}
	}
	trainingVector[0] = hamMails;
	trainingVector[0 + m + 1] = spamMails;
	trainingVector[trainingVector.size() - 1] = hamMails + spamMails;
}