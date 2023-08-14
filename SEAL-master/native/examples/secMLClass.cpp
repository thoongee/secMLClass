// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include <ctime> // jina add for clock()
#include <time.h>                // for time()
#include <stdlib.h>              // for rand/srand
#include <cmath>
#include <iostream>
#include <fstream>          // for ifstream
#include <cctype>           // for tolower(), isalpha()
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <cstddef>
#include <memory>
#include <limits>
#include <typeinfo>
#include <cstdlib>
#include <string>
#include "randgen.h"
#include "strutils.h"
#include "NaiveBayesBreastCancer.h"
#include "MultinomialNB_Email.h"
#include "seal/seal.h"

#include <experimental/filesystem>
namespace fs = std::experimental;

using namespace std;
using namespace seal;

/*
Helper function: Prints the name of the example in a fancy banner.
*/
void print_example_banner(string title)
{
    if (!title.empty())
    {
        size_t title_length = title.length();
        size_t banner_length = title_length + 2 + 2 * 10;
        string banner_top(banner_length, '*');
        string banner_middle = string(10, '*') + " " + title + " " + string(10, '*');

        cout << endl
            << banner_top << endl
            << banner_middle << endl
            << banner_top << endl
            << endl;
    }
}

/*
Helper function: Prints the parameters in a SEALContext.
*/
void print_parameters(shared_ptr<SEALContext> context)
{
    // Verify parameters
    if (!context)
    {
        throw invalid_argument("context is not set");
    }
    auto &context_data = *context->context_data();

    /*
    Which scheme are we using?
    */
    string scheme_name;
    switch (context_data.parms().scheme())
    {
    case scheme_type::BFV:
        scheme_name = "BFV";
        break;
    case scheme_type::CKKS:
        scheme_name = "CKKS";
        break;
    default:
        throw invalid_argument("unsupported scheme");
    }

    cout << "/ Encryption parameters:" << endl;
    cout << "| scheme: " << scheme_name << endl;
    cout << "| poly_modulus_degree: " << 
        context_data.parms().poly_modulus_degree() << endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    cout << "| coeff_modulus size: " << context_data.
        total_coeff_modulus_bit_count() << " bits" << endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == scheme_type::BFV)
    {
        cout << "| plain_modulus: " << context_data.
            parms().plain_modulus().value() << endl;
    }

    cout << "\\ noise_standard_deviation: " << context_data.
        parms().noise_standard_deviation() << endl;
    cout << endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
ostream &operator <<(ostream &stream, parms_id_type parms_id)
{
    stream << hex << parms_id[0] << " " << parms_id[1] << " "
        << parms_id[2] << " " << parms_id[3] << dec;
    return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template<typename T>
void print_vector(vector<T> vec, size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);

    size_t slot_count = vec.size();

    cout << fixed << setprecision(prec) << endl;
    if(slot_count <= 2 * print_size)
    {
        cout << "    [";
        for (size_t i = 0; i < slot_count; i++)
        {
            cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(max(vec.size(), 2 * print_size));
        cout << "    [";
        for (size_t i = 0; i < print_size; i++)
        {
            cout << " " << vec[i] << ",";
        }
        if(vec.size() > 2 * print_size)
        {
            cout << " ...,";
        }
        for (size_t i = slot_count - print_size; i < slot_count; i++)
        {
            cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    cout << endl;

    /*
    Restore the old std::cout formatting.
    */
    cout.copyfmt(old_fmt);
}


void SVM_LR_ClientCentric(int polyModulus, const vector<int64_t>& svm_weights, const vector<int64_t> &logistic_weights,
	const vector<vector<int64_t>> &data, GaloisKeys &gal_keys, Decryptor &decryptor, const vector<int> &true_label,
	BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator, 
	int NrOfQueryInstances, int classNotation)
{
	vector<int64_t> ones_v(polyModulus, 0);
	Plaintext ones_p;
	//Ciphertext temp;

	for (int i = 0; i < polyModulus; i += 16)
		ones_v[i] = 1;
	batch_encoder.encode(ones_v, ones_p);

	Plaintext M_svm_p, M_LR_p; //SVM's and LR's trained models - M in plaintext
	Ciphertext M_svm_c, M_LR_c; //SVM and LR's trained models trained model - M in ciphertext
	batch_encoder.encode(svm_weights, M_svm_p);
	encryptor.encrypt(M_svm_p, M_svm_c);
	batch_encoder.encode(logistic_weights, M_LR_p);
	encryptor.encrypt(M_LR_p, M_LR_c);

	if (polyModulus == 16384)
	{
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);

		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
	}
	else if (polyModulus == 8192)
	{
		evaluator.mod_switch_to_next_inplace(M_svm_c);
		evaluator.mod_switch_to_next_inplace(M_svm_c);

		evaluator.mod_switch_to_next_inplace(M_LR_c);
		evaluator.mod_switch_to_next_inplace(M_LR_c);
	}

	vector<Plaintext> Sdata_p(data.size()); //the set of user query - S in plaintext 
	vector<Ciphertext> Sdata_c(Sdata_p.size());//the set of user query - S in ciphertext

	vector<Ciphertext> CMS_svm_c(Sdata_p.size());//the classific. of S according to M as a ciphertext
	vector<Plaintext> CMS_svm_p(CMS_svm_c.size());//the classific. of S according to M as a plaintext
	vector<vector<int64_t>> CMS_svm_v(CMS_svm_p.size(), vector<int64_t>(polyModulus, 0));

	vector<Ciphertext> CMS_LR_c(Sdata_p.size());//the classific. of S according to M as a ciphertext
	vector<Plaintext> CMS_LR_p(CMS_LR_c.size());//the classific. of S according to M as a plaintext
	vector<vector<int64_t>> CMS_LR_v(CMS_LR_p.size(), vector<int64_t>(polyModulus, 0));
	Ciphertext tmp_c;
	for (int i = 0; i < Sdata_p.size(); i++)
	{

		std::cout << "Sdata i=" << i << endl;
		batch_encoder.encode(data[i], Sdata_p[i]);
		//encryptor.encrypt(SdataBCancer_p[i], Sdata_c[i]);

		//SVM
		//std::cout << "Dealing with SVM now" << endl;
		evaluator.multiply_plain(M_svm_c, Sdata_p[i], CMS_svm_c[i]);

		evaluator.rotate_rows(CMS_svm_c[i], 1, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 2, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 4, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 8, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.multiply_plain_inplace(CMS_svm_c[i], ones_p);

		//std::cout << "The noise budget of CMS_svm_c[" << i << "] is " << decryptor.invariant_noise_budget(CMS_svm_c[i]) << endl;
		decryptor.decrypt(CMS_svm_c[i], CMS_svm_p[i]);
		batch_encoder.decode(CMS_svm_p[i], CMS_svm_v[i]);

		//LR
		//std::cout << "Dealing with LR now" << endl;
		evaluator.multiply_plain(M_LR_c, Sdata_p[i], CMS_LR_c[i]);

		evaluator.rotate_rows(CMS_LR_c[i], 1, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 2, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 4, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 8, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.multiply_plain_inplace(CMS_LR_c[i], ones_p);

		//std::cout << "The noise budget of CMS_LR_c[" << i << "] is " << decryptor.invariant_noise_budget(CMS_LR_c[i]) << endl;
		decryptor.decrypt(CMS_LR_c[i], CMS_LR_p[i]);
		batch_encoder.decode(CMS_LR_p[i], CMS_LR_v[i]);
	}
	vector<int> predEncLabel_svm(NrOfQueryInstances); //labels according to the secMLCLass for SVM
	vector<int> predEncLabel_LR(NrOfQueryInstances); //logistic regression
	int index = 0;
	int secMLClass_SVM_acc = 0, secMLClass_LR_acc = 0;
	for (int i = 0; i < CMS_svm_v.size(); i++)
	{
		for (int j = 0; j < CMS_svm_v[0].size(); j += 16)
		{
			//SVM
			if (CMS_svm_v[i][j] >= 0)
			{
				predEncLabel_svm[index] = classNotation;
			}
			else
			{
				predEncLabel_svm[index] = classNotation-1;
			}
			if (predEncLabel_svm[index] == true_label[index])
			{
				secMLClass_SVM_acc++;
			}

			//LR
			if (CMS_LR_v[i][j] >= 0)
			{
				predEncLabel_LR[index] = classNotation;
			}
			else
			{
				predEncLabel_LR[index] = classNotation-1;
			}
			if (predEncLabel_LR[index] == true_label[index])
			{
				secMLClass_LR_acc++;
			}
			index++;
			if (index == predEncLabel_svm.size())
				break;
		}
		if (index == predEncLabel_svm.size())
			break;
	}
	std::cout << endl << "secMLClass SVM Accuracy: " << (double)secMLClass_SVM_acc * 100 / NrOfQueryInstances << "%" << "\n";
	std::cout << "secMLClass LR Accuracy: " << (double)secMLClass_LR_acc * 100 / NrOfQueryInstances << "%" << "\n";
}

void SVM_LR_ServerCentric(int polyModulus, const vector<int64_t>& svm_weights, const vector<int64_t> &logistic_weights,
	vector<vector<int64_t>> &data, GaloisKeys &gal_keys, Decryptor &decryptor, const vector<int> &true_label,
	BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator, 
	int NrOfQueryInstances, int classNotation)
{
	vector<int64_t> ones_v(polyModulus, 0);
	Plaintext ones_p;
	//Ciphertext temp;

	for (int i = 0; i < polyModulus; i += 16)
		ones_v[i] = 1;
	batch_encoder.encode(ones_v, ones_p);

	Plaintext M_svm_p, M_LR_p; //SVM's and LR's trained models - M in plaintext
	Ciphertext M_svm_c, M_LR_c; //SVM and LR's trained models trained model - M in ciphertext
	batch_encoder.encode(svm_weights, M_svm_p);
	encryptor.encrypt(M_svm_p, M_svm_c);
	batch_encoder.encode(logistic_weights, M_LR_p);
	encryptor.encrypt(M_LR_p, M_LR_c);

	vector<Plaintext> Sdata_p(data.size()); //the set of user query - S in plaintext 
	vector<Ciphertext> Sdata_c(Sdata_p.size());//the set of user query - S in ciphertext

	vector<Ciphertext> CMS_svm_c(Sdata_p.size());//the classific. of S according to M as a ciphertext
	vector<Plaintext> CMS_svm_p(CMS_svm_c.size());//the classific. of S according to M as a plaintext
	vector<vector<int64_t>> CMS_svm_v(CMS_svm_p.size(), vector<int64_t>(polyModulus, 0));

	vector<Ciphertext> CMS_LR_c(Sdata_p.size());//the classific. of S according to M as a ciphertext
	vector<Plaintext> CMS_LR_p(CMS_LR_c.size());//the classific. of S according to M as a plaintext
	vector<vector<int64_t>> CMS_LR_v(CMS_LR_p.size(), vector<int64_t>(polyModulus, 0));
	Ciphertext tmp_c;

	for (int i = 0; i < Sdata_p.size(); i++)
	{
		std::cout << "Sdata i=" << i << endl;
		batch_encoder.encode(data[i], Sdata_p[i]);
		encryptor.encrypt(Sdata_p[i], Sdata_c[i]);

		if (polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
		}
		else if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
			evaluator.mod_switch_to_next_inplace(Sdata_c[i]);
		}

		//SVM
		//std::cout << "Dealing with SVM now" << endl;
		evaluator.multiply_plain(Sdata_c[i], M_svm_p, CMS_svm_c[i]);

		evaluator.rotate_rows(CMS_svm_c[i], 1, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 2, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 4, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.rotate_rows(CMS_svm_c[i], 8, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_svm_c[i], tmp_c);

		evaluator.multiply_plain_inplace(CMS_svm_c[i], ones_p);

		//std::cout << "The noise budget of CMS_svm_ADFA_c[" << i << "] is " << decryptor.invariant_noise_budget(CMS_svm_c[i]) << endl;
		decryptor.decrypt(CMS_svm_c[i], CMS_svm_p[i]);
		batch_encoder.decode(CMS_svm_p[i], CMS_svm_v[i]);

		//LR
		//std::cout << "Dealing with LR now" << endl;
		evaluator.multiply_plain(Sdata_c[i], M_LR_p, CMS_LR_c[i]);

		evaluator.rotate_rows(CMS_LR_c[i], 1, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 2, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 4, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.rotate_rows(CMS_LR_c[i], 8, gal_keys, tmp_c);
		evaluator.add_inplace(CMS_LR_c[i], tmp_c);

		evaluator.multiply_plain_inplace(CMS_LR_c[i], ones_p);

		//std::cout << "The noise budget of CMS_LR_ADFA_c[" << i << "] is " << decryptor.invariant_noise_budget(CMS_LR_c[i]) << endl;
		decryptor.decrypt(CMS_LR_c[i], CMS_LR_p[i]);
		batch_encoder.decode(CMS_LR_p[i], CMS_LR_v[i]);
	}
	vector<int> predEncLabel_svm(NrOfQueryInstances); //labels according to the secMLCLass for SVM
	vector<int> predEncLabel_LR(NrOfQueryInstances); //logistic regression
	int index = 0;
	int secMLClass_SVM_acc = 0, secMLClass_LR_acc = 0;
	for (int i = 0; i < CMS_svm_v.size(); i++)
	{
		for (int j = 0; j < CMS_svm_v[0].size(); j += 16)
		{
			//SVM
			if (CMS_svm_v[i][j] >= 0)
			{
				predEncLabel_svm[index] = classNotation;
			}
			else
			{
				predEncLabel_svm[index] = classNotation-1;
			}
			if (predEncLabel_svm[index] == true_label[index])
			{
				secMLClass_SVM_acc++;
			}

			//LR
			if (CMS_LR_v[i][j] >= 0)
			{
				predEncLabel_LR[index] = classNotation;
			}
			else
			{
				predEncLabel_LR[index] = classNotation-1;
			}
			if (predEncLabel_LR[index] == true_label[index])
			{
				secMLClass_LR_acc++;
			}
			index++;
			if (index == predEncLabel_svm.size())
				break;
		}
		if (index == predEncLabel_svm.size())
			break;
	}
	std::cout << endl << "secMLClass SVM Accuracy: " << (double)secMLClass_SVM_acc * 100 / NrOfQueryInstances << "%" << "\n";
	std::cout << "secMLClass LR Accuracy: " << (double)secMLClass_LR_acc * 100 / NrOfQueryInstances << "%" << "\n";

}


void NB_BreastC_ServerCentric(int polyModulus, int plainModulus, std::shared_ptr<SEALContext> context, GaloisKeys &gal_keys, Decryptor &decryptor,
	BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator,
	vector<int64_t> trainedModelVector, int batchQueryCount, const vector<vector<int>> &localDataset, ofstream &output)
{
	//int polyModulus = batchQueryCount * 256;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "BEGGINING OF SERVER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "BEGGINING OF SERVER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	output << endl << "***********************************************************************" << endl;

	///// maxDifferenceOfLogsOfProbs=4468, maxDiffIndex=286, constK = 127 /////
	int maxDifferenceOfLogsOfProbs = 4468;
	vector<int64_t> batchedRVector(polyModulus, 0);
	vector<int64_t> batchedHVector(polyModulus, 0);
	vector<int64_t> batchedTrainedModel(polyModulus, 0);
	RandGen rndGen;
	
	for (int i = 0; i < batchQueryCount; i++)
	{
		batchedRVector[128 * i] = rndGen.RandInt(1, plainModulus / (2 * (maxDifferenceOfLogsOfProbs + 1) - 1));

		batchedHVector[128 * i] = rndGen.RandInt(0, batchedRVector[128 * i]);
		for (int j = 0; j < 128; j++)
		{
			batchedTrainedModel[i * 128 + j] = trainedModelVector[j];
		}
	}

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "*****PHASE 1: LOCALLY CONSTRUCTING QUERY(IES) AT USER(S)*****" << endl;
	std::cout << "*********ENCRYPTING THEM AND SENDING THEM TO TACS************" << endl;
	std::cout << "**************************(USER 1)***************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;;
	output << "*****PHASE 1: LOCALLY CONSTRUCTING QUERY(IES) AT USER(S)*****" << endl;
	output << "*********ENCRYPTING THEM AND SENDING THEM TO TACS************" << endl;
	output << "**************************(USER 1)***************************" << endl;


	//auto context = SEALContext::Create(parms);
	//print_parameters(context);
	//auto qualifiers = context->context_data()->qualifiers();
	//std::cout << "Batching enabled inside the FUNCTION: " << boolalpha << qualifiers.using_batching << endl;
	//KeyGenerator keygen(context);
	//auto public_key = keygen.public_key();
	//auto secret_key = keygen.secret_key();
	//auto gal_keys = keygen.galois_keys(DefaultParams::dbc_max());
	//auto relin_keys = keygen.relin_keys(DefaultParams::dbc_max());
	//Encryptor encryptor(context, public_key);
	//Evaluator evaluator(context);
	//Decryptor decryptor(context, secret_key);
	//BatchEncoder batch_encoder(context);

	Plaintext batchedTrainedModelPlain, batchedRPlain, batchedHPlain;
	Ciphertext batchedTrainedModelCipher;

	std::cout << "HERE 1" << endl;
	batch_encoder.encode(batchedTrainedModel, batchedTrainedModelPlain);
	encryptor.encrypt(batchedTrainedModelPlain, batchedTrainedModelCipher);

	/*std::cout << "Size of batchedTrainedModelCipher=" << sizeof(batchedTrainedModelCipher) << endl;
	std::cout << "Size of int=" << sizeof(int) << endl;
	std::cout << "Size of batchedFinalLogOfProbsVector="<<sizeof(batchedTrainedModel) << endl;
	*/

	std::cout << "HERE 2" << endl;
	batch_encoder.encode(batchedRVector, batchedRPlain);
	batch_encoder.encode(batchedHVector, batchedHPlain);

	std::cout << "Now constructing query(ies) at USER 1, encoding and encrypting them ... (BLACK 5)" << endl;
	output << "Now constructing query(ies) at USER 1, encoding and encrypting them ... (BLACK 5)" << endl;

	auto context_data = context->context_data();
	
	if (batchQueryCount >= 64)
	{
		while (context_data->next_context_data()->next_context_data())
		{
			cout << "Chain index: " << context_data->chain_index() << endl;
			cout << "parms_id of batchedTrainedModelCipher: " << batchedTrainedModelCipher.parms_id() << endl;
			std::cout << "The size of the coefficient modulus is: "
				<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
			cout << "Noise budget at this level: "
				<< decryptor.invariant_noise_budget(batchedTrainedModelCipher) << " bits" << endl;
			cout << "\\" << endl;
			cout << " \\-->" << endl;
			evaluator.mod_switch_to_next_inplace(batchedTrainedModelCipher);
			context_data = context_data->next_context_data();
		}
	}
	cout << "Chain index: " << context_data->chain_index() << endl;
	cout << "parms_id of batchedTrainedModelCipher: " << batchedTrainedModelCipher.parms_id() << endl;
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	cout << "Noise budget at this level: "
		<< decryptor.invariant_noise_budget(batchedTrainedModelCipher) << " bits" << endl;
	
	std::cout << "USER 1 sends his data to TACS ... (WHITE 4)" << endl;
	output << "USER 1 sends his data to TACS ... (WHITE 4)" << endl;

	int bitTransmitedUSER_TO_TACS_KB = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedUSER_TO_TACS_KB << " KB \nto transmit from the USER 1 to TACS server ...(WHITE 4)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "In all we have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedUSER_TO_TACS_KB << " KB \nto transmit from the USER 1 to TACS server ...(WHITE 4)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;;
	std::cout << "*****PHASE 2: HOMOMORPHICALLY PROCESSING AND RANDOMIZING*****" << endl;
	std::cout << "*********THE QUERY AT TACS AND SENDING IT TO EDS*************" << endl;
	std::cout << "***************************(TACS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;;
	output << "*****PHASE 2: HOMOMORPHICALLY PROCESSING AND RANDOMIZING*****" << endl;
	output << "*********THE QUERY AT TACS AND SENDING IT TO EDS*************" << endl;
	output << "***************************(TACS)****************************" << endl;

	std::cout << "TACS homomorphically processes the query ... (BLACK 6)" << endl;
	output << "TACS homomorphically processes the query ... (BLACK 6)" << endl;


	Ciphertext tmpBatchedTMC = batchedTrainedModelCipher;
	while (context_data->next_context_data())
	{
		cout << "Chain index: " << context_data->chain_index() << endl;
		cout << "parms_id of batchedTrainedModelCipher: " << tmpBatchedTMC.parms_id() << endl;
		std::cout << "The size of the coefficient modulus is: "
			<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
		cout << "Noise budget at this level: "
			<< decryptor.invariant_noise_budget(tmpBatchedTMC) << " bits" << endl;
		cout << "\\" << endl;
		cout << " \\-->" << endl;
		evaluator.mod_switch_to_next_inplace(tmpBatchedTMC);
		context_data = context_data->next_context_data();
	}

	cout << "Chain index: " << context_data->chain_index() << endl;
	cout << "parms_id of batchedTrainedModelCipher: " << tmpBatchedTMC.parms_id() << endl;
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	cout << "Noise budget at this level: "
		<< decryptor.invariant_noise_budget(tmpBatchedTMC) << " bits" << endl;


	std::cout << "TACS sends his processed randomized data to EDS ... (WHITE 5)" << endl;
	output << "TACS sends his processed randomized data to EDS ... (WHITE 5)" << endl;

	int bitTransmitedTACS_TO_EDS_KB = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedTACS_TO_EDS_KB << " KB \nto transmit from TACS to EDS ...(WHITE 5)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedTACS_TO_EDS_KB << " KB \nto transmit from TACS to EDS ...(WHITE 5)" << endl;


	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	std::cout << "****************AND SENDING IT TO THE USER 1*****************" << endl;
	std::cout << "****************************(EDS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	output << "****************AND SENDING IT TO THE USER 1*****************" << endl;
	output << "****************************(EDS)****************************" << endl;

	std::cout << "EDS decrypts the final randomized classification ... (BLACK 7)" << endl;
	output << "EDS decrypts the final randomized classification ... (BLACK 7)" << endl;

	std::cout << "EDS sends the randomized result in plain to USER 1 ... (WHITE 6)" << endl;
	output << "EDS sends the randomized result in plain to USER 1 ... (WHITE 6)" << endl;

	std::cout << "We have batchQueryCount*log(Nr_classes)=" << batchQueryCount << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 1, which is negligable ... (WHITE 6)" << endl;
	output << "We have batchQueryCount*log(Nr_classes)=" << batchQueryCount << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 1, which is negligable ... (WHITE 6)" << endl;

	int totalKBClassification = bitTransmitedUSER_TO_TACS_KB + bitTransmitedTACS_TO_EDS_KB;

	std::cout << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the server centric classification" << endl
		<< "with batchQueryCount=" << batchQueryCount << " queries per batch" << endl;

	output << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the server centric classification" << endl
		<< "with batchQueryCount=" << batchQueryCount << " queries per batch" << endl;

	std::cout << endl << endl;
	std::cout << "*********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 4: THE USER 1 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	std::cout << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	std::cout << "****************************(USER 1)**************************" << endl;
	output << endl << endl;
	output << "*********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 4: THE USER 1 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	output << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	output << "****************************(USER 1)**************************" << endl;

	std::cout << "USER 1 de-randomizes the result to get the final classification ... (BLACK 8) " << endl;
	output << "USER 1 de-randomizes the result to get the final classification ... (BLACK 8) " << endl;

	std::cout << "For the Server Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of " << batchQueryCount << " batched queries" << endl;
	output << "For the Server Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of batchQueryCount=" << batchQueryCount << " batched queries" << endl;

	//std::system("pause");

	vector<vector<int>> confusionMatrix{ {0,0},{0,0} };

	chrono::high_resolution_clock::time_point time_start, time_end;
	long long int black5MicroSeconds = 0, black6MicroSeconds = 0, black7MicroSeconds = 0, black8MicroSeconds = 0;
	int queryCount = 0;

	for (int i = 0; i < (localDataset.size() - batchQueryCount); i += batchQueryCount)
	{
		std::cout << "Dealing with i=" << i << endl;
		queryCount++;
		vector <int> realClass(batchQueryCount, 0);
		vector <int> predictedClass(batchQueryCount, 0);
		vector<int64_t> queryVector(polyModulus, 0);
		Plaintext queryVectorPlain;
		Ciphertext queryVectorCipher;

		//cout<<"BEGGINING OF PHASE1 - BLACK 5"<<endl;

		time_start = chrono::high_resolution_clock::now();
		time_end = chrono::high_resolution_clock::now();
		auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

		for (int j = 0; j < batchQueryCount; j++)
		{
			time_start = chrono::high_resolution_clock::now();
			realClass[j] = localDataset[i + j][localDataset[0].size() - 1];
			queryVector[j * 128 + 0] = 1; 
			//queryVector[j * 256 + 128] = 1;

			for (int k = 0; k < localDataset[0].size() - 1; k++)
			{
				queryVector[j * 128 + 10 * k + localDataset[i + j][k]] = 1;
				//queryVector[j * 256 + 10 * k + localDataset[i + j][k] + 128] = 1;
			}
			time_end = chrono::high_resolution_clock::now();
			time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
			black5MicroSeconds += time_diff.count();
		}
		time_start = chrono::high_resolution_clock::now();
		batch_encoder.encode(queryVector, queryVectorPlain);
		encryptor.encrypt(queryVectorPlain, queryVectorCipher);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black5MicroSeconds += time_diff.count();

		evaluator.mod_switch_to_inplace(queryVectorCipher, batchedTrainedModelCipher.parms_id());

		//cout<<"END OF PHASE1 - BLACK 5"<<endl; 

		//cout<<"BEGGINING OF PHASE2 - BLACK 6"<<endl;

		Ciphertext tmpQueryResult;
		int step = 0;

		time_start = chrono::high_resolution_clock::now();
		evaluator.multiply_plain_inplace(queryVectorCipher, batchedTrainedModelPlain);
		//evaluator.relinearize_inplace(queryVectorCipher, relin_keys);

		for (int j = 0; j < log2(128); j++)
		{
			step = pow(2, j);
			tmpQueryResult = queryVectorCipher;
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(queryVectorCipher, tmpQueryResult);
		}
		//evaluator.rotate_rows(queryVectorCipher, 128, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(queryVectorCipher, tmpQueryResult);//[A]-[B],
		evaluator.multiply_plain_inplace(queryVectorCipher, batchedRPlain);//([A]-[B])R
		evaluator.add_plain_inplace(queryVectorCipher, batchedHPlain);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, batchedHPlain);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black6MicroSeconds += time_diff.count();

		//cout<<"END OF PHASE2 - BLACK 6"<<endl; 

		std::cout << "the noise budget before switching to next mod inplace is: " << decryptor.invariant_noise_budget(queryVectorCipher) << endl;
		evaluator.mod_switch_to_inplace(queryVectorCipher, tmpBatchedTMC.parms_id());

		Plaintext queryVectorResultPlain;
		vector<int64_t> queryVectorResult;
		std::cout << "the noise budget before decrypting is: " << decryptor.invariant_noise_budget(queryVectorCipher) << endl;

		//cout<<"BEGGINING OF PHASE3 - BLACK 7"<<endl;		
		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(queryVectorCipher, queryVectorResultPlain);
		batch_encoder.decode(queryVectorResultPlain, queryVectorResult);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black7MicroSeconds += time_diff.count();
		//cout<<"END OF PHASE3 - BLACK 7"<<endl;


		//cout<<"BEGGINING OF PHASE4 - BLACK 8"<<endl;
		time_start = chrono::high_resolution_clock::now();
		for (int queryCount = 0; queryCount < batchQueryCount; queryCount++)
		{
			if (queryVectorResult[128 * queryCount] > 0)
			{
				predictedClass[queryCount] = 1;
			}
			else
			{
				predictedClass[queryCount] = 2;
			}
		}
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black8MicroSeconds += time_diff.count();
		//cout<<"END OF PHASE4 - BLACK 8"<<endl;

		//Finding the accuracy (confusion matrix)
		for (int k = 0; k < predictedClass.size(); k++)
		{
			if (realClass[k] == 1)
			{
				if (predictedClass[k] == 1)
					confusionMatrix[0][0]++;
				else
					confusionMatrix[0][1]++;
			}
			else
			{
				if (predictedClass[k] == 1)
					confusionMatrix[1][0]++;
				else
					confusionMatrix[1][1]++;
			}
		}
	}

	std::cout << "Printing the confusion matrix:" << endl;
	std::cout << "\t\tpred_NonCancer    pred Cancer" << endl;
	std::cout << "real_NonCancer\t" << confusionMatrix[0][0] << "\t\t  " << confusionMatrix[0][1] << endl;
	std::cout << "real Cancer\t" << confusionMatrix[1][0] << "\t\t  " << confusionMatrix[1][1] << endl;
	output << "Printing the confusion matrix:" << endl;
	output << "\t\tpred_NonCancer    pred Cancer" << endl;
	output << "real_NonCancer\t" << confusionMatrix[0][0] << "\t\t  " << confusionMatrix[0][1] << endl;
	output << "real Cancer\t" << confusionMatrix[1][0] << "\t\t  " << confusionMatrix[1][1] << endl;

	// jina add
	long blackTotalMicroSeconds = 0;
	blackTotalMicroSeconds = black5MicroSeconds+black6MicroSeconds+black7MicroSeconds+black8MicroSeconds;

	std::cout << "black5MicroSeconds =" << black5MicroSeconds << endl;
	std::cout << "black6MicroSeconds =" << black6MicroSeconds << endl;
	std::cout << "black7MicroSeconds =" << black7MicroSeconds << endl;
	std::cout << "black8MicroSeconds =" << black8MicroSeconds << endl;
	output << "black5MicroSeconds =" << black5MicroSeconds << endl;
	output << "black6MicroSeconds =" << black6MicroSeconds << endl;
	output << "black7MicroSeconds =" << black7MicroSeconds << endl;
	output << "black8MicroSeconds =" << black8MicroSeconds << endl;

	std::cout << endl << "The average time for constructing and encrypting the query(ies) at USER 1 (BLACK 5) is "
		<< (long double)black5MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << "The average time for query(ies) processing at TACS (BLACK 6) is "
		<< (long double)black6MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << endl << "The average time for decrypting the query(ies) at EDS (BLACK 7) is "
		<< (long double)black7MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << "The average time for finding final classifications of query at USER 1 (BLACK 8) is "
		<< (long double)black8MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	
	output << endl << "The average time for constructing and encrypting the query(ies) at USER 1 (BLACK 5) is "
		<< (long double)black5MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << "The average time for query(ies) processing at TACS (BLACK 6) is "
		<< (long double)black6MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << endl << "The average time for decrypting the query(ies) at EDS (BLACK 7) is "
		<< (long double)black7MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << "The average time for finding final classifications of query at USER 1 (BLACK 8) is "
		<< (long double)black8MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;

	std::cout << "blackTotalMicroSeconds =" << blackTotalMicroSeconds << endl;
	output << "blackTotalMicroSeconds =" << blackTotalMicroSeconds << endl;
	std::cout << endl << "The average time for Total is "
		<< (long double)blackTotalMicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << "query count =" << queryCount << endl;
	output << "query count =" << queryCount << endl;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "END OF SERVER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "END OF SERVER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	output << endl << "***********************************************************************" << endl;
}

void NB_BreastC_ClientCentric(int polyModulus, int plainModulus, std::shared_ptr<SEALContext> context, GaloisKeys &gal_keys, Decryptor &decryptor,
	BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator, 
	vector<int64_t> trainedModelVector, int batchQueryCount, const vector<vector<int>> &localDataset, ofstream &output)
{
	//int polyModulus = batchQueryCount * 256;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "BEGGINING OF USER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "BEGGINING OF USER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	output << endl << "***********************************************************************" << endl;

	///// maxDifferenceOfLogsOfProbs=4468, maxDiffIndex=286, constK = 127 /////
	int maxDifferenceOfLogsOfProbs = 4468;
	vector<int64_t> batchedRVector(polyModulus, 0);
	vector<int64_t> batchedHVector(polyModulus, 0);
	vector<int64_t> batchedTrainedModel(polyModulus, 0);
	//vector<uint64_t> batchedSTCGuardVector(polyModulus, 0);
	RandGen rndGen;

	//uint64_t plainModulus;//needed after decryption comparison
	//if (polyModulus == 4096)
	//	plainModulus = 65537; //65537, 6144001
	//else
	//	plainModulus = 65536524289;
	//EncryptionParameters parms(scheme_type::BFV);
	//parms.set_poly_modulus_degree(polyModulus); //4096, 8192, 16384,
	//parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));
	//parms.set_plain_modulus(plainModulus);
	//for (int i = 0; i < 128; i++)
	//{
	//	batchedSTCGuardVector[i] = rndGen.RandInt(1, plainModulus / (maxDifferenceOfLogsOfProbs + 10));
	//	batchedSTCGuardVector[i + 128] = batchedSTCGuardVector[i];
	//}

	for (int i = 0; i < batchQueryCount; i++)
	{
		batchedRVector[128 * i] = rndGen.RandInt(1, plainModulus / (2 * (maxDifferenceOfLogsOfProbs + 1) - 1));

		batchedHVector[128 * i] = rndGen.RandInt(0, batchedRVector[128 * i]);
		for (int j = 0; j < 128; j++)
		{
			batchedTrainedModel[i * 128 + j] = trainedModelVector[j];
			//batchedSTCGuardVector[i * 256 + j] = batchedSTCGuardVector[j];
		}
	}

	//for (int i = 0; i < 128; i++)
	//{
	//	batchedSTCGuardVector[i] = rndGen.RandInt(1, plainModulus / (maxDifferenceOfLogsOfProbs + 10));
	//	batchedSTCGuardVector[i + 128] = batchedSTCGuardVector[i];
	//}

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "*****PHASE 1: TACS RANDOMIZES THE ENCYPTED TRAINED MODEL*****" << endl;
	std::cout << "******IN ORDER TO PROTECT FROM STC, SENDS IT TO USER 2*******" << endl;
	std::cout << "***************************(TACS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "*****PHASE 1: TACS RANDOMIZES THE ENCYPTED TRAINED MODEL*****" << endl;
	output << "******IN ORDER TO PROTECT FROM STC, SENDS IT TO USER 2*******" << endl;
	output << "***************************(TACS)****************************" << endl;

	//auto context = SEALContext::Create(parms);
	//print_parameters(context);
	//auto qualifiers = context->context_data()->qualifiers();
	//std::cout << "Batching enabled inside the FUNCTION: " << boolalpha << qualifiers.using_batching << endl;

	//KeyGenerator keygen(context);
	//auto public_key = keygen.public_key();
	//auto secret_key = keygen.secret_key();
	//auto gal_keys = keygen.galois_keys(DefaultParams::dbc_max());
	//auto relin_keys = keygen.relin_keys(DefaultParams::dbc_max());
	//Encryptor encryptor(context, public_key);
	//Evaluator evaluator(context);
	//Decryptor decryptor(context, secret_key);
	//BatchEncoder batch_encoder(context);
	
	Plaintext batchedTrainedModelPlain, batchedRPlain, batchedHPlain, batchedSTCguardPlain;
	Ciphertext batchedTrainedModelCipher;

	batch_encoder.encode(batchedTrainedModel, batchedTrainedModelPlain);
	encryptor.encrypt(batchedTrainedModelPlain, batchedTrainedModelCipher);


	/*std::cout << "Size of batchedTrainedModelCipher=" << sizeof(batchedTrainedModelCipher) << endl;
	std::cout << "Size of int=" << sizeof(int) << endl;
	std::cout << "Size of batchedFinalLogOfProbsVector="<<sizeof(batchedTrainedModel) << endl;
	*/
	//batch_encoder.encode(batchedSTCGuardVector, batchedSTCguardPlain);
	batch_encoder.encode(batchedRVector, batchedRPlain);
	batch_encoder.encode(batchedHVector, batchedHPlain);

	//std::cout << "Adding the STCguard to the trained model at TACS ... (BLACK 9)" << endl;
	//std::cout << "Adding the STCguard to the trained model at TACS ... (BLACK 9)" << endl;

	auto context_data = context->context_data();

	if (batchQueryCount >= 64)
	{
		while (context_data->next_context_data()->next_context_data())
		{
			cout << "Chain index: " << context_data->chain_index() << endl;
			cout << "parms_id of batchedTrainedModelCipherr: " << batchedTrainedModelCipher.parms_id() << endl;
			std::cout << "The size of the coefficient modulus is: "
				<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
			cout << "Noise budget at this level: "
				<< decryptor.invariant_noise_budget(batchedTrainedModelCipher) << " bits" << endl;
			cout << "\\" << endl;
			cout << " \\-->" << endl;
			evaluator.mod_switch_to_next_inplace(batchedTrainedModelCipher);
			context_data = context_data->next_context_data();
		}
	}
	cout << "Chain index: " << context_data->chain_index() << endl;
	cout << "parms_id of batchedTrainedModelCipher: " << batchedTrainedModelCipher.parms_id() << endl;
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	cout << "Noise budget at this level: "
		<< decryptor.invariant_noise_budget(batchedTrainedModelCipher) << " bits" << endl;


	std::cout << "TACS sends the randomized trained model to USER 2 ... (WHITE 7)" << endl;
	output << "TACS sends the randomized trained model to USER 2 ... (WHITE 7)" << endl;

	int bitTransmitedTACS_TO_USER_KB = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedTACS_TO_USER_KB << " KB \nto transmit from TACS to USER 2 ... (WHITE 7)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "In all we have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedTACS_TO_USER_KB << " KB \nto transmit from TACS to USER 2 ...(WHITE 7)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---*************************" << endl;
	std::cout << "**********PHASE 2: HOMOMORPHICALLY PROCESSING THE ************" << endl;
	std::cout << "*********THE QUERY AT USER 2 AND SENDING IT TO EDS************" << endl;
	std::cout << "***************************(USER 2)***************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---*************************" << endl;
	output << "**********PHASE 2: HOMOMORPHICALLY PROCESSING THE ************" << endl;
	output << "*********THE QUERY AT USER 2 AND SENDING IT TO EDS************" << endl;
	output << "***************************(USER 2)***************************" << endl;

	std::cout << "USER 2 homomorphically processes the query and randomizes it ... (BLACK 10)" << endl;
	output << "USER 2 homomorphically processes the query and randomizes it ... (BLACK 10)" << endl;

	std::cout << "USER 2 sends his processed randomized data to EDS ... (WHITE 8)" << endl;
	output << "USER 2 sends his processed randomized data to EDS ... (WHITE 8)" << endl;

	Ciphertext tmpBatchedFLOPC = batchedTrainedModelCipher;
	while (context_data->next_context_data())
	{
		cout << "Chain index: " << context_data->chain_index() << endl;
		cout << "parms_id of tmpBatchedFLOPC: " << tmpBatchedFLOPC.parms_id() << endl;
		std::cout << "The size of the coefficient modulus is: "
			<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
		cout << "Noise budget at this level: "
			<< decryptor.invariant_noise_budget(tmpBatchedFLOPC) << " bits" << endl;
		cout << "\\" << endl;
		cout << " \\-->" << endl;
		evaluator.mod_switch_to_next_inplace(tmpBatchedFLOPC);
		context_data = context_data->next_context_data();
	}

	cout << "Chain index: " << context_data->chain_index() << endl;
	cout << "parms_id of tmpBatchedFLOPC: " << tmpBatchedFLOPC.parms_id() << endl;
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	cout << "Noise budget at this level: "
		<< decryptor.invariant_noise_budget(tmpBatchedFLOPC) << " bits" << endl;


	int bitTransmitedUSER_TO_EDS_KB = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedUSER_TO_EDS_KB << " KB \nto transmit from USER 2 to EDS ...(WHITE 8)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< bitTransmitedUSER_TO_EDS_KB << " KB \nto transmit from USER 2 to EDS ...(WHITE 8)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	std::cout << "****************AND SENDING IT TO THE USER 2*****************" << endl;
	std::cout << "****************************(EDS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	output << "****************AND SENDING IT TO THE USER 2*****************" << endl;
	output << "****************************(EDS)****************************" << endl;

	std::cout << "EDS decrypts the final randomized classification ... (BLACK 11)" << endl;
	output << "EDS decrypts the final randomized classification ... (BLACK 11)" << endl;

	std::cout << "EDS sends the randomized result in plain to USER 2 ... (WHITE 9)" << endl;
	output << "EDS sends the randomized result in plain to USER 2 ... (WHITE 9)" << endl;

	std::cout << "We have batchQueryCount*log(Nr_classes)=" << batchQueryCount << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 2, which is negligable ... (WHITE 9)" << endl;
	output << "We have batchQueryCount*log(Nr_classes)=" << batchQueryCount << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 2, which is negligable ... (WHITE 9)" << endl;

	int totalKBClassification = bitTransmitedTACS_TO_USER_KB + bitTransmitedUSER_TO_EDS_KB;

	std::cout << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the user centric classification" << endl
		<< "with batchQueryCount=" << batchQueryCount << " queries per batch" << endl;

	output << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the user centric classification" << endl
		<< "with batchQueryCount=" << batchQueryCount << " queries per batch" << endl;

	std::cout << endl << endl;
	std::cout << "*********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 4: THE USER 2 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	std::cout << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	std::cout << "****************************(USER 2)**************************" << endl;
	output << endl << endl;
	output << "*********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 4: THE USER 2 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	output << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	output << "****************************(USER 2)**************************" << endl;

	std::cout << "USER 1 de-randomizes the result to get the final classification ... (BLACK 12) " << endl;
	output << "USER 1 de-randomizes the result to get the final classification ... (BLACK 12) " << endl;

	std::cout << "For the User Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of " << batchQueryCount << " batched queries" << endl;
	output << "For the User Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of batchQueryCount=" << batchQueryCount << " batched queries" << endl;

	vector<vector<int>> confusionMatrix{ {0,0},{0,0} };

	chrono::high_resolution_clock::time_point time_start, time_end;
	long long int black9MicroSeconds = 0, black10MicroSeconds = 0, black11MicroSeconds = 0, black12MicroSeconds = 0;
	int queryCount = 0;

	for (int i = 0; i < (localDataset.size() - batchQueryCount); i += batchQueryCount)
	{
		std::cout << "Dealing with i=" << i << endl;
		queryCount++;
		vector <int> realClass(batchQueryCount, 0);
		vector <int> predictedClass(batchQueryCount, 0);
		vector<int64_t> queryVector(polyModulus, 0);
		Plaintext queryVectorPlain;

		time_start = chrono::high_resolution_clock::now();
		time_end = chrono::high_resolution_clock::now();
		auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

		//cout<<" BEGGINNING OF PHASE1 - BLACK 9"<<endl;
		Ciphertext tmpBatchedTrainedModelCipher = batchedTrainedModelCipher;
		//time_start = chrono::high_resolution_clock::now();
		//evaluator.add_plain_inplace(tmpBatchedTrainedModelCipher, batchedSTCguardPlain);
		//time_end = chrono::high_resolution_clock::now();
		//time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		//black9MicroSeconds += time_diff.count();
		//cout<<" END OF PHASE1 - BLACK 9"<<endl;

		//cout<<" BEGGINNING OF PHASE2 - BLACK 10"<<endl;
		Ciphertext tmpQueryResult;
		int step = 0;
		for (int j = 0; j < batchQueryCount; j++)
		{
			time_start = chrono::high_resolution_clock::now();
			realClass[j] = localDataset[i + j][localDataset[0].size() - 1];
			queryVector[j * 128 + 0] = 1; 
			//queryVector[j * 256 + 128] = 1;

			for (int k = 0; k < localDataset[0].size() - 1; k++)
			{
				queryVector[j * 128 + 10 * k + localDataset[i + j][k]] = 1;
				//queryVector[j * 256 + 10 * k + localDataset[i + j][k] + 128] = 1;
			}
			time_end = chrono::high_resolution_clock::now();
			time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
			black10MicroSeconds += time_diff.count();
		}
		time_start = chrono::high_resolution_clock::now();
		batch_encoder.encode(queryVector, queryVectorPlain);
		evaluator.multiply_plain_inplace(tmpBatchedTrainedModelCipher, queryVectorPlain);
		//evaluator.relinearize_inplace(tmpBatchedTrainedModelCipher, relin_keys);

		for (int j = 0; j < log2(128); j++)
		{
			step = pow(2, j);
			tmpQueryResult = tmpBatchedTrainedModelCipher;
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(tmpBatchedTrainedModelCipher, tmpQueryResult);
		}
		//evaluator.rotate_rows(tmpBatchedTrainedModelCipher, 128, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(tmpBatchedTrainedModelCipher, tmpQueryResult);//[A]-[B],
		evaluator.multiply_plain_inplace(tmpBatchedTrainedModelCipher, batchedRPlain);//([A]-[B])R
		evaluator.add_plain_inplace(tmpBatchedTrainedModelCipher, batchedHPlain);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, batchedHPlain);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black10MicroSeconds += time_diff.count();


		std::cout << "the noise budget before switching to next mod inplace is: " << decryptor.invariant_noise_budget(tmpBatchedTrainedModelCipher) << endl;
		evaluator.mod_switch_to_inplace(tmpBatchedTrainedModelCipher, tmpBatchedFLOPC.parms_id());

		Plaintext queryVectorResultPlain;
		vector<int64_t> queryVectorResult;
		std::cout << "the noise budget before decrypting is: " << decryptor.invariant_noise_budget(tmpBatchedTrainedModelCipher) << endl;
		//cout<<" END OF PHASE2 - BLACK 10"<<endl;

		//cout<<" BEGGINNING OF PHASE 3 - BLACK 11"<<endl;
		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(tmpBatchedTrainedModelCipher, queryVectorResultPlain);
		batch_encoder.decode(queryVectorResultPlain, queryVectorResult);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black11MicroSeconds += time_diff.count();
		//cout<<" END OF PHASE 3 - BLACK 11"<<endl;

		//cout<<" BEGGINNING OF PHASE 4 - BLACK 12"<<endl;
		time_start = chrono::high_resolution_clock::now();
		for (int queryCount = 0; queryCount < batchQueryCount; queryCount++)
		{
			if (queryVectorResult[128 * queryCount] > 0)
			{
				predictedClass[queryCount] = 1;
			}
			else
			{
				predictedClass[queryCount] = 2;
			}
		}
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		black12MicroSeconds += time_diff.count();
		//cout<<"END OF PHASE4 - BLACK 8"<<endl;

		//Finding the accuracy (confusion matrix)
		for (int k = 0; k < predictedClass.size(); k++)
		{
			if (realClass[k] == 1)
			{
				if (predictedClass[k] == 1)
					confusionMatrix[0][0]++;
				else
					confusionMatrix[0][1]++;
			}
			else
			{
				if (predictedClass[k] == 1)
					confusionMatrix[1][0]++;
				else
					confusionMatrix[1][1]++;
			}
		}
	}

	std::cout << "Printing the confusion matrix:" << endl;
	std::cout << "\t\tpred_NonCancer    pred Cancer" << endl;
	std::cout << "real_NonCancer\t" << confusionMatrix[0][0] << "\t\t  " << confusionMatrix[0][1] << endl;
	std::cout << "real Cancer\t" << confusionMatrix[1][0] << "\t\t  " << confusionMatrix[1][1] << endl;
	output << "Printing the confusion matrix:" << endl;
	output << "\t\tpred_NonCancer    pred Cancer" << endl;
	output << "real_NonCancer\t" << confusionMatrix[0][0] << "\t\t  " << confusionMatrix[0][1] << endl;
	output << "real Cancer\t" << confusionMatrix[1][0] << "\t\t  " << confusionMatrix[1][1] << endl;

	std::cout << "black9MicroSeconds =" << black9MicroSeconds << endl;
	std::cout << "black10MicroSeconds =" << black10MicroSeconds << endl;
	std::cout << "black11MicroSecond =" << black11MicroSeconds << endl;
	std::cout << "black12MicroSeconds =" << black12MicroSeconds << endl;
	output << "black9MicroSeconds =" << black9MicroSeconds << endl;
	output << "black10MicroSeconds =" << black10MicroSeconds << endl;
	output << "black11MicroSecond =" << black11MicroSeconds << endl;
	output << "black12MicroSeconds =" << black12MicroSeconds << endl;

	std::cout << endl << "The average time for STC guard at TACS (BLACK 9) is "
		<< (long double)black9MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << "The average time for query(ies) processing at USER 2 (BLACK 10) is "
		<< (long double)black10MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << endl << "The average time for decrypting the query(ies) at EDS (BLACK 11) is "
		<< (long double)black11MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	std::cout << "The average time for finding final classifications of query at USER 2 (BLACK 12) is "
		<< (long double)black12MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << endl << "The average time for STC guard at TACS (BLACK 9) is "
		<< (long double)black9MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << "The average time for query(ies) processing at USER 2 (BLACK 10) is "
		<< (long double)black10MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << endl << "The average time for decrypting the query(ies) at EDS (BLACK 11) is "
		<< (long double)black11MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;
	output << "The average time for finding final classifications of query at USER 2 (BLACK 12) is "
		<< (long double)black12MicroSeconds / ((long double)queryCount*batchQueryCount) << " microseconds per query." << endl;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "END OF USER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "END OF USER CENTRIC Processing of " << localDataset.size() / batchQueryCount << " queries with batchQuerySize=" << batchQueryCount << " queries per batch" << endl;
	output << endl << "***********************************************************************" << endl;
}

void DeepLearning_BreastC_ServerCentric(int polyModulus, const vector<double> & DLdata,
	vector<vector<int>> DL_weights1, const vector<int> & DL_weights2, vector<int> DL_true_label)
{
	EncryptionParameters parms(scheme_type::BFV);
		
	int untilForRndNumber = 200, dbc = 40;
	uint64_t tmpPlainModulus = 40961;
	
	if (polyModulus == 4096)
	{
		tmpPlainModulus = 40961;
		untilForRndNumber = 200;
		dbc = 60;
	}
	else if (polyModulus == 8192)
	{
		tmpPlainModulus = 65536524289; //65536524289 -> 37 bits; 65537 -> 15 bits, 65537
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}
	else
	{
		tmpPlainModulus = 65536524289; //60 bits 576460752304439297, 557057 -> 20 bits
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}

	parms.set_poly_modulus_degree(polyModulus);//4096, 8192, 16384, 32768 
	parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));//4096, 8192, 16384, 32768 

	parms.set_plain_modulus(tmpPlainModulus);//40961 for 4096, and 65536524289 for 8192, 576460752304439297 for others

	auto context = SEALContext::Create(parms);

	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();
	std::cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(dbc);

	auto relin_keys = keygen.relin_keys(dbc);

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);

	vector<vector<int>> DL_weights1Transp(16, vector<int>(16, 0));
	for (int i = 0; i < DL_weights1.size(); i++)
	{
		for (int j = 0; j < DL_weights1[0].size(); j++)
		{
			DL_weights1Transp[i][j] = DL_weights1[j][i];
		}
	}
	DL_weights1 = DL_weights1Transp;

	vector<vector<int64_t>> SdataDL_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	std::cout << "SdataDL_v.size()=" << SdataDL_v.size() << ", SdataDL_v[0].size()=" << SdataDL_v[0].size() << endl;

	vector<int64_t> M_DL_weights1_v(polyModulus, 0);
	vector<int64_t> M_DL_weights2_v(polyModulus, 0);

	vector<Plaintext> SdataDL_p(SdataDL_v.size());
	Plaintext M_DL_weights1_p, M_DL_weights2_p;

	vector<Ciphertext> SdataDL_c(SdataDL_v.size());
	Ciphertext M_DL_weights1_c, M_DL_weights2_c;

	//packing the SdataDL_v and the corresponding plaintexts and ciphertexts now
	for (int i = 0; i < 699; i++)
	{
		vector<int64_t> tempData(16, 0);
		for (int j = 0; j < tempData.size(); j++)
		{
			tempData[j] = (int64_t)DLdata[i * 16 + j];
		}

		for (int k = 0; k < tempData.size(); k++)
		{			
			for (int j = 0; j < tempData.size(); j++)
			{
				SdataDL_v[(i * 256) / polyModulus][(i * 256 + k * 16 + j) % polyModulus] = tempData[j];
			}
		}
	}


	for (int i = 0; i < SdataDL_p.size(); i++)
	{
		//std::cout << "Printing vector SdataDL_v[" << i << "]:" << endl;
		//print_vector(SdataDL_v[i], 1536);
		batch_encoder.encode(SdataDL_v[i], SdataDL_p[i]);
		encryptor.encrypt(SdataDL_p[i], SdataDL_c[i]);
		if (polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(SdataDL_c[i]);
			evaluator.mod_switch_to_next_inplace(SdataDL_c[i]);
			evaluator.mod_switch_to_next_inplace(SdataDL_c[i]);
			evaluator.mod_switch_to_next_inplace(SdataDL_c[i]);			
		}		
	}

	//packing the M_DL_weights1_v and the corresponding plaintexts and ciphertexts
	int queriesPerCiphertext = polyModulus / 256;
	for (int i = 0; i < queriesPerCiphertext; i++)
	{
		for (int j = 0; j < DL_weights1.size(); j++)
		{
			for (int k = 0; k < DL_weights1[0].size(); k++)
			{
				M_DL_weights1_v[i * 256 + j * DL_weights1[0].size() + k] = (int64_t)DL_weights1[j][k];
			}
		}
	}

	//std::cout << "Printing vector M_DL_weights1_v:" << endl;
	//print_vector(M_DL_weights1_v, 1536);
	batch_encoder.encode(M_DL_weights1_v, M_DL_weights1_p);
	encryptor.encrypt(M_DL_weights1_p, M_DL_weights1_c);

	//packing the M_DL_weights2_v and the corresponding plaintexts and ciphertexts
	for (int i = 0; i < queriesPerCiphertext; i++)
	{
		for (int j = 0; j < DL_weights2.size(); j++)
		{
			M_DL_weights2_v[i * 256 + j * 16] = (int64_t)DL_weights2[j];
		}
	}
	//std::cout << "Printing vector M_DL_weights2_v:" << endl;
	//print_vector(M_DL_weights2_v, 1536);
	batch_encoder.encode(M_DL_weights2_v, M_DL_weights2_p);
	encryptor.encrypt(M_DL_weights2_p, M_DL_weights2_c);
	//std::system("pause");


	///////////////////////////////////////////////////////////////////   	  
	std::cout << "DOING CHECKING NOW BEFORE THE REAL STUFF" << endl;
	//vector<vector<int64_t>> SdataDL_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resMultW1_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resSqr_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resMultW2_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	//vector<vector<int64_t>> resMultW1_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));

	double maxDLRes = 0;
	for (int i = 0; i < SdataDL_v.size(); i++)
	{
		for (int j = 0; j < SdataDL_v[0].size(); j += 16)
		{
			for (int k = 0; k < 16; k++)
			{
				resMultW1_v[i][j] += SdataDL_v[i][j + k] * M_DL_weights1_v[j + k];
			}
			resSqr_v[i][j] = pow(resMultW1_v[i][j], 2);
			if (maxDLRes < abs(resSqr_v[i][j]))
				maxDLRes = resSqr_v[i][j];
		}
	}

	int index = 0;
	vector<int> testLabels(699);

	for (int i = 0; i < SdataDL_v.size(); i++)
	{
		for (int j = 0; j < SdataDL_v[0].size(); j += 256)
		{
			for (int k = 0; k < 256; k += 16)
			{
				resMultW2_v[i][j] += resSqr_v[i][j + k] * M_DL_weights2_v[j + k];
			}
			
			if (maxDLRes < abs(resMultW2_v[i][j]))
				maxDLRes = resMultW2_v[i][j];
			if (resMultW2_v[i][j] >= 0)
				testLabels[index] = 2;
			else
				testLabels[index] = 1;
			index++;
			if (index == 699)
				break;
		}
	}	

	double acc = 0;
	for (int i = 0; i < 699; i++)
	{
		if (testLabels[i] == DL_true_label[i])
		{
			acc++;
		}		
	}

	cout << "Accuracy before the real stuff: " << acc * 100 / 699 << "%" << endl;
	cout << "maxDLRes=" << maxDLRes << endl;
	//std::system("pause");
	////////////////////////////////////////////////////////////////////


	std::cout << "Doing the secMlClass now:" << endl;
	Ciphertext tmpQueryResult;

	vector<Plaintext> SM_result_p(SdataDL_p.size());
	vector<vector<int64_t>> SM_result_v(SdataDL_p.size(), vector<int64_t>(polyModulus, 0));

	for (int i = 0; i < SdataDL_p.size(); i++)
	{
		std::cout << "Doing the multiplication with the hidden layer for packed query " << i << endl;
		evaluator.multiply_plain_inplace(SdataDL_c[i], M_DL_weights1_p);
		for (int j = 0; j < log2(16); j++)
		{
			int step = pow(2, j);
			tmpQueryResult = SdataDL_c[i];
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(SdataDL_c[i], tmpQueryResult);
		}

		std::cout << "Doing the squaring now for packed query " << i << endl;
		
		///////////////
		evaluator.mod_switch_to_next_inplace(SdataDL_c[i]);
		////////////////////		
		
		evaluator.square_inplace(SdataDL_c[i]);
		evaluator.relinearize_inplace(SdataDL_c[i], relin_keys);
		
		std::cout << "Doing the multiplication with the output layer for packed query " << i << endl;
		evaluator.multiply_plain_inplace(SdataDL_c[i], M_DL_weights2_p);
		for (int j = 4; j < log2(256); j++)
		{
			int step = pow(2, j);
			tmpQueryResult = SdataDL_c[i];
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(SdataDL_c[i], tmpQueryResult);
		}
		std::cout << "The noise before decryption is: " << decryptor.invariant_noise_budget(SdataDL_c[i]) << " bits" << endl;

		decryptor.decrypt(SdataDL_c[i], SM_result_p[i]);
		batch_encoder.decode(SM_result_p[i], SM_result_v[i]);
	}

	vector <int64_t> resDLSecML(699);
	vector<int> predLabDL(699);
	acc = 0;
	for (int i = 0; i < 699; i++)
	{
		int64_t tmpResult = SM_result_v[(i * 256) / polyModulus][(i * 256) % polyModulus];
		//cout << "tmpResult[" << i << "]=" << tmpResult << endl;
		resDLSecML[i] = tmpResult;
		if (tmpResult >= 0)
		{
			predLabDL[i] = 2;
		}
		else
		{
			predLabDL[i] = 1;
		}

		if (predLabDL[i] == DL_true_label[i])
		{
			acc++;
		}
	}
	cout << "Accuracy: " << acc * 100 / 699 << "%" << endl;

}

void DeepLearning_BreastC_ClientCentric (int polyModulus, const vector<double> & DLdata, 
	 vector<vector<int>> DL_weights1, const vector<int> & DL_weights2, vector<int> DL_true_label)
{	
	EncryptionParameters parms(scheme_type::BFV);
		
	int untilForRndNumber = 200, dbc = 40;
	uint64_t tmpPlainModulus = 40961;
		
	if (polyModulus == 4096)
	{
		tmpPlainModulus = 40961;
		untilForRndNumber = 200;
		dbc = 60;
	}
	else if (polyModulus == 8192)
	{
		tmpPlainModulus = 65536524289; //65536524289 -> 37 bits; 65537 -> 15 bits, 65537
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}
	else
	{
		tmpPlainModulus = 65536524289; //60 bits 576460752304439297, 557057 -> 20 bits
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}

	parms.set_poly_modulus_degree(polyModulus);//4096, 8192, 16384, 32768 
	parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));//4096, 8192, 16384, 32768 

	parms.set_plain_modulus(tmpPlainModulus);//40961 for 4096, and 65536524289 for 8192, 576460752304439297 for others

	auto context = SEALContext::Create(parms);

	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();
	std::cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(dbc);

	auto relin_keys = keygen.relin_keys(dbc);

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);	
	
	vector<vector<int>> DL_weights1Transp(16, vector<int>(16, 0));
	for (int i = 0; i < DL_weights1.size(); i++)
	{
		for (int j = 0; j < DL_weights1[0].size(); j++)
		{
			DL_weights1Transp[i][j] = DL_weights1[j][i];
		}
	}
	DL_weights1 = DL_weights1Transp;

	vector<vector<int64_t>> SdataDL_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	std::cout << "SdataDL_v.size()=" << SdataDL_v.size() << ", SdataDL_v[0].size()=" << SdataDL_v[0].size() << endl;

	vector<int64_t> M_DL_weights1_v(polyModulus, 0);
	vector<int64_t> M_DL_weights2_v(polyModulus, 0);

	vector<Plaintext> SdataDL_p(SdataDL_v.size());
	Plaintext M_DL_weights1_p, M_DL_weights2_p;

	vector<Ciphertext> SdataDL_c(SdataDL_v.size());
	Ciphertext M_DL_weights1_c, M_DL_weights2_c;

	//packing the SdataDL_v and the corresponding plaintexts and ciphertexts now
	for (int i = 0; i < 699; i++)
	{
		vector<int64_t> tempData(16, 0);
		for (int j = 0; j < tempData.size(); j++)
		{
			tempData[j] = (int64_t)DLdata[i * 16 + j];
		}

		for (int k = 0; k < tempData.size(); k++)
		{
			//if(((i * 256) / polyModulus) >= SdataDL_v.size())
			//	cout << "(i * 256) / polyModulus) = " << (i * 256) / polyModulus << endl;			

			for (int j = 0; j < tempData.size(); j++)
			{
				//if (((i * 256 + k * 16 + j) % polyModulus) >= SdataDL_v[0].size())
				//	cout << "((i * 256 + k * 16 + j) % polyModulus) = "
				//	<< ((i * 256 + k * 16 + j) % polyModulus) << endl;
				SdataDL_v[(i * 256) / polyModulus][(i * 256 + k * 16 + j) % polyModulus] = tempData[j];
			}
		}
	}


	for (int i = 0; i < SdataDL_p.size(); i++)
	{
		//std::cout << "Printing vector SdataDL_v[" << i << "]:" << endl;
		//print_vector(SdataDL_v[i], 1536);
		batch_encoder.encode(SdataDL_v[i], SdataDL_p[i]);
		encryptor.encrypt(SdataDL_p[i], SdataDL_c[i]);
	}

	//packing the M_DL_weights1_v and the corresponding plaintexts and ciphertexts
	int queriesPerCiphertext = polyModulus / 256;

	for (int i = 0; i < queriesPerCiphertext; i++)
	{
		for (int j = 0; j < DL_weights1.size(); j++)
		{
			for (int k = 0; k < DL_weights1[0].size(); k++)
			{
				M_DL_weights1_v[i * 256 + j * DL_weights1[0].size() + k] = (int64_t)DL_weights1[j][k];
			}
		}
	}

	//std::cout << "Printing vector M_DL_weights1_v:" << endl;
	//print_vector(M_DL_weights1_v, 1536);
	batch_encoder.encode(M_DL_weights1_v, M_DL_weights1_p);
	encryptor.encrypt(M_DL_weights1_p, M_DL_weights1_c);
	
	//packing the M_DL_weights2_v and the corresponding plaintexts and ciphertexts
	for (int i = 0; i < queriesPerCiphertext; i++)
	{
		for (int j = 0; j < DL_weights2.size(); j++)
		{
			M_DL_weights2_v[i * 256 + j * 16] = (int64_t)DL_weights2[j];
		}
	}
	//std::cout << "Printing vector M_DL_weights2_v:" << endl;
	//print_vector(M_DL_weights2_v, 1536);
	batch_encoder.encode(M_DL_weights2_v, M_DL_weights2_p);
	encryptor.encrypt(M_DL_weights2_p, M_DL_weights2_c);

	if (polyModulus == 16384)
	{
		evaluator.mod_switch_to_next_inplace(M_DL_weights1_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights1_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights1_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights1_c);
		
		evaluator.mod_switch_to_next_inplace(M_DL_weights2_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights2_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights2_c);
		evaluator.mod_switch_to_next_inplace(M_DL_weights2_c);		
		
	}	
	//std::system("pause");

	
	///////////////////////////////////////////////////////////////////   	  
	std::cout << "DOING CHECKING NOW BEFORE THE REAL STUFF" << endl;
	//vector<vector<int64_t>> SdataDL_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resMultW1_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resSqr_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	vector<vector<int64_t>> resMultW2_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));
	//vector<vector<int64_t>> resMultW1_v(ceil(((double)699 * 256) / polyModulus), vector<int64_t>(polyModulus, 0));

	double maxDLRes = 0;
	for (int i = 0; i < SdataDL_v.size(); i++)
	{
		for (int j = 0; j < SdataDL_v[0].size(); j += 16)
		{
			for (int k = 0; k < 16; k++)
			{
				resMultW1_v[i][j] += SdataDL_v[i][j + k] * M_DL_weights1_v[j + k];
			}
			resSqr_v[i][j] = pow(resMultW1_v[i][j], 2);
			if (maxDLRes < abs(resSqr_v[i][j]))
				maxDLRes = resSqr_v[i][j];
		}
	}

	int index = 0;
	vector<int> testLabels(699);

	for (int i = 0; i < SdataDL_v.size(); i++)
	{
		for (int j = 0; j < SdataDL_v[0].size(); j += 256)
		{
			for (int k = 0; k < 256; k += 16)
			{
				resMultW2_v[i][j] += resSqr_v[i][j + k] * M_DL_weights2_v[j + k];
			}

			//for (int k = 0; k < 256; k+=16)
			//{
			//	resMultW2_v[i][j+k] = resSqr_v[i][j + k] * M_DL_weights2_v[j + k];
			//}

			if (maxDLRes < abs(resMultW2_v[i][j]))
				maxDLRes = resMultW2_v[i][j];
			if (resMultW2_v[i][j] >= 0)
				testLabels[index] = 2;
			else
				testLabels[index] = 1;
			index++;
			if (index == 699)
				break;
		}
	}

	/*for (int i = 0; i < SdataDL_v.size(); i++)
	{
		std::cout << "Printing resMultW1_v[" << i << "]:" << endl;
		print_vector(resMultW1_v[i], 1536);

		std::cout << "Printing resSqr_v[" << i << "]:" << endl;
		print_vector(resSqr_v[i], 1536);

		std::cout << "Printing resMultW2_v[" << i << "]:" << endl;
		print_vector(resMultW2_v[i], 1536);
	}*/

	double acc = 0;
	for (int i = 0; i < 699; i++)
	{
		if (testLabels[i] == DL_true_label[i])
		{
			acc++;
		}
		/*if (testLabels[i] != DL_true_label[i])
		{
			std::cout << "Results are NOT the same for i=" << i << endl;

			for (int j = 0; j < 256; j+=16)
			{
				std::cout << "resMultW1_v[" << i << "]=" << resMultW1_v[(i * 256) / polyModulus][(i * 256 + j) % polyModulus] << endl;
				std::cout << "resSqr_v[" << i << "]=" << resSqr_v[(i * 256) / polyModulus][(i * 256 + j) % polyModulus] << endl << endl;
			}

			std::cout << "resMultW2_v["<<i<<"]=" << resMultW2_v[(i * 256) / polyModulus][(i * 256) % polyModulus]<<endl;
			std::cout << "res_mlp[" << i << "]=" << res_mlp[i] << endl<<endl;
			std::system("pause");
		}
		else
		{
			std::cout << "Results are the same for i=" << i << endl;

			for (int j = 0; j < 256; j += 16)
			{
				std::cout << "resMultW1_v[" << i << "]=" << resMultW1_v[(i * 256) / polyModulus][(i * 256 + j) % polyModulus] << endl;
				std::cout << "resSqr_v[" << i << "]=" << resSqr_v[(i * 256) / polyModulus][(i * 256 + j) % polyModulus] << endl << endl;
			}
			std::cout << "resMultW2_v[" << i << "]=" << resMultW2_v[(i * 256) / polyModulus][(i * 256) % polyModulus] << endl;
			std::cout << "res_mlp[" << i << "]=" << res_mlp[i] << endl << endl;
			std::system("pause");
		}*/
	}

	cout << "Accuracy before the real stuff: " << acc * 100 / 699 << "%" << endl;
	cout << "maxDLRes=" << maxDLRes << endl;
	//std::system("pause");
	////////////////////////////////////////////////////////////////////
	   	  

	std::cout << "Doing the secMlClass now:" << endl;
	Ciphertext tmpQueryResult;

	vector<Plaintext> SM_result_p(SdataDL_p.size());
	vector<vector<int64_t>> SM_result_v(SdataDL_p.size(), vector<int64_t>(polyModulus, 0));

	evaluator.mod_switch_to_next_inplace(M_DL_weights2_c);
	Ciphertext tmpResult;
	for (int i = 0; i < SdataDL_p.size(); i++)
	{
		std::cout << "Doing the multiplication with the hidden layer for packed query " << i << endl;
		tmpResult = M_DL_weights1_c;
		evaluator.multiply_plain_inplace(tmpResult, SdataDL_p[i]);
		for (int j = 0; j < log2(16); j++)
		{
			int step = pow(2, j);
			tmpQueryResult = tmpResult;
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(tmpResult, tmpQueryResult);
		}

		std::cout << "Doing the squaring now for packed query " << i << endl;

		//////////////////
		evaluator.mod_switch_to_next_inplace(tmpResult);
		
		//////////////////

		evaluator.square_inplace(tmpResult);
		evaluator.relinearize_inplace(tmpResult, relin_keys);

		std::cout << "Doing the multiplication with the output layer for packed query " << i << endl;
		evaluator.multiply_inplace(tmpResult, M_DL_weights2_c);
		evaluator.relinearize_inplace(tmpResult, relin_keys);
		for (int j = 4; j < log2(256); j++)
		{
			int step = pow(2, j);
			tmpQueryResult = tmpResult;
			evaluator.rotate_rows_inplace(tmpQueryResult, step, gal_keys);
			evaluator.add_inplace(tmpResult, tmpQueryResult);
		}
		std::cout << "The noise before decryption is: " << decryptor.invariant_noise_budget(tmpResult) << " bits" << endl;

		decryptor.decrypt(tmpResult, SM_result_p[i]);
		batch_encoder.decode(SM_result_p[i], SM_result_v[i]);
	}

	vector <int64_t> resDLSecML(699);
	vector<int> predLabDL(699);
	acc = 0;
	for (int i = 0; i < 699; i++)
	{
		int64_t tmpResult = SM_result_v[(i * 256) / polyModulus][(i * 256) % polyModulus];
		//cout << "tmpResult[" << i << "]=" << tmpResult << endl;
		resDLSecML[i] = tmpResult;
		if (tmpResult >= 0)
		{
			predLabDL[i] = 2;
		}
		else
		{
			predLabDL[i] = 1;
		}

		if (predLabDL[i] == DL_true_label[i])
		{
			acc++;
		}
	}
	cout << "Accuracy: " << acc * 100 / 699 << "%" << endl;

}

void PPServCenMNClass(string path, const vector<string>& finalSelectedFeaturesDecrypted,
	vector<int64_t> trainedModel_v, int polyModulus, int m, ofstream &output)
{

	int countQuery = 0, queriesPerCphrtxt = polyModulus / ((m + 1));

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "BEGGINING OF SERVER CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "BEGGINING OF SERVER CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	output << endl << "***********************************************************************" << endl;

	long long int totalTimeBlack12 = 0, totalTimeBlack13 = 0, totalTimeBlack14 = 0, totalTimeBlack15 = 0,
		sizeWhite12 = 0, sizeWhite13 = 0, sizeWhite14 = 0, sizeWhite15 = 0;

	vector<vector<int>> confusionMatrix{ {0,0},{0,0} };
	string word;
	string hamPath = path + "\\ham", filename;
	vector<string> hamFiles, spamFiles;

	vector<int64_t> R_v(polyModulus, 0);
	vector<int64_t> h_v(polyModulus, 0);
	RandGen rndGen;

	uint64_t plainModulus;//needed after decryption comparison
	if (polyModulus == 4096)
		plainModulus = 65537; //65537, 6144001
	if (polyModulus == 8192)
		plainModulus = 65536524289; //65536524289, 3288203001857,
	if (polyModulus == 16384)
		plainModulus = 65536524289; //65536524289, 3288203001857, 576460752302473217

	for (int i = 0; i < queriesPerCphrtxt; i++)
	{
		R_v[(m + 1) * i] = rndGen.RandInt(1, 4096);
		h_v[(m + 1) * i] = rndGen.RandInt(-R_v[(m + 1) * i], R_v[(m + 1) * i]);

		//R_v[2*(m + 1) * i + m + 1] = rndGen.RandInt(1, 1024);
		//h_v[2*(m + 1) * i + m + 1] = rndGen.RandInt(0, R_v[2*(m + 1) * i + m + 1]);		
	}

	EncryptionParameters parms(scheme_type::BFV);

	parms.set_poly_modulus_degree(polyModulus); //4096, 8192, 16384,
	parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));
	parms.set_plain_modulus(plainModulus);

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "*****PHASE 1: LOCALLY CONSTRUCTING QUERY(IES) AT USER(S)*****" << endl;
	std::cout << "*********ENCRYPTING THEM AND SENDING THEM TO TACS************" << endl;
	std::cout << "**************************(USER 1)***************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;;
	output << "*****PHASE 1: LOCALLY CONSTRUCTING QUERY(IES) AT USER(S)*****" << endl;
	output << "*********ENCRYPTING THEM AND SENDING THEM TO TACS************" << endl;
	output << "**************************(USER 1)***************************" << endl;

	auto context = SEALContext::Create(parms);
	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();

	std::cout << "Batching enabled inside the FUNCTION: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(DefaultParams::dbc_max());

	auto relin_keys = keygen.relin_keys(DefaultParams::dbc_max());

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);

	Plaintext trainedModel_p, R_p, h_p;
	Ciphertext trainedModel_c;

	batch_encoder.encode(trainedModel_v, trainedModel_p);
	encryptor.encrypt(trainedModel_p, trainedModel_c);

	/*std::cout << "Size of trainedModel_c=" << sizeof(trainedModel_c) << endl;
	std::cout << "Size of int=" << sizeof(int) << endl;
	std::cout << "Size of batchedFinalLogOfProbsVector="<<sizeof(batchedTrainedModel) << endl;
	*/

	batch_encoder.encode(R_v, R_p);
	batch_encoder.encode(h_v, h_p);

	std::cout << "Now constructing query(ies) at USER 1, encoding and encrypting them ... (BLACK 12)" << endl;
	output << "Now constructing query(ies) at USER 1, encoding and encrypting them ... (BLACK 12)" << endl;

	auto context_data = context->context_data();
	
	if (polyModulus == 8192)
	{
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();

		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		//evaluator.mod_switch_to_next_inplace(trainedModel_c);
	}
	else //if(polyModulus == 16384)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();

		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		//evaluator.mod_switch_to_next_inplace(trainedModel_c);
	}


	std::cout << "USER 1 sends his data to TACS ... (WHITE 12)" << endl;
	output << "USER 1 sends his data to TACS ... (WHITE 12)" << endl;

	sizeWhite12 = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= sizeWhite12="
		<< sizeWhite12 << " KB \nto transmit from the USER 1 to TACS server ...(WHITE 12)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "In all we have 2*coeff_modulus_bit_count()*poly_modulus_degree()= sizeWhite12"
		<< sizeWhite12 << " KB \nto transmit from the USER 1 to TACS server ...(WHITE 12)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;;
	std::cout << "*****PHASE 2: HOMOMORPHICALLY PROCESSING AND RANDOMIZING*****" << endl;
	std::cout << "*********THE QUERY AT TACS AND SENDING IT TO EDS*************" << endl;
	std::cout << "***************************(TACS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;;
	output << "*****PHASE 2: HOMOMORPHICALLY PROCESSING AND RANDOMIZING*****" << endl;
	output << "*********THE QUERY AT TACS AND SENDING IT TO EDS*************" << endl;
	output << "***************************(TACS)****************************" << endl;

	std::cout << "TACS homomorphically processes the query ... (BLACK 13)" << endl;
	output << "TACS homomorphically processes the query ... (BLACK 13)" << endl;

	Ciphertext tmpTV_c = trainedModel_c;
	if (polyModulus == 8192)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
	}
	else //if(polyModulus == 16384)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();
	}
 // IMPROVEDCLASSIFICATION 

	std::cout << "TACS sends his processed randomized data to EDS ... (WHITE 13)" << endl;
	output << "TACS sends his processed randomized data to EDS ... (WHITE 13)" << endl;

	sizeWhite13 = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< sizeWhite13 << " KB \nto transmit from TACS to EDS ...(WHITE 13)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()= "
		<< sizeWhite13 << " KB \nto transmit from TACS to EDS ...(WHITE 13)" << endl;


	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	std::cout << "****************AND SENDING IT TO THE USER 1*****************" << endl;
	std::cout << "****************************(EDS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	output << "****************AND SENDING IT TO THE USER 1*****************" << endl;
	output << "****************************(EDS)****************************" << endl;

	std::cout << "EDS decrypts the final randomized classification ... (BLACK 14)" << endl;
	output << "EDS decrypts the final randomized classification ... (BLACK 14)" << endl;

	std::cout << "EDS sends the randomized result in plain to USER 1 ... (WHITE 14)" << endl;
	output << "EDS sends the randomized result in plain to USER 1 ... (WHITE 14)" << endl;

	std::cout << "We havequeriesPerCphrtxt*log2(Nr_classes)=" << queriesPerCphrtxt << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 1, which is negligable ... (WHITE 14)" << endl;
	output << "We havequeriesPerCphrtxt*log2(Nr_classes)=" << queriesPerCphrtxt << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 1, which is negligable ... (WHITE 14)" << endl;

	int totalKBClassification = sizeWhite12 + sizeWhite13;

	std::cout << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the server centric classification" << endl
		<< "with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries per batch" << endl;

	output << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the server centric classification" << endl
		<< "with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries per batch" << endl;

	std::cout << endl << endl;
	std::cout << "*********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 4: THE USER 1 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	std::cout << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	std::cout << "****************************(USER 1)**************************" << endl;
	output << endl << endl;
	output << "*********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 4: THE USER 1 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	output << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	output << "****************************(USER 1)**************************" << endl;

	std::cout << "USER 1 de-randomizes the result to get the final classification ... (BLACK 8) " << endl;
	output << "USER 1 de-randomizes the result to get the final classification ... (BLACK 8) " << endl;

	std::cout << "For the Server Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of " << queriesPerCphrtxt << " batched queries" << endl;
	output << "For the Server Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query ofqueriesPerCphrtxt=" << queriesPerCphrtxt << " batched queries" << endl;

	//ACTUAL PP SERVER CENTRIC CLASSIFICATION	  	  
	vector <int> powersOfTwo(log2(m + 1), 0);
	for (int i = 0; i < powersOfTwo.size(); i++)
	{
		powersOfTwo[i] = pow(2, i);
	}

	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

	Plaintext queryVector_p;
	Ciphertext queryVector_c;

	for (const auto & entry : fs::filesystem::recursive_directory_iterator(hamPath))
	{
		filename = entry.path().string();
		hamFiles.push_back(filename);
	}

	for (int j = 0; j < hamFiles.size(); j += queriesPerCphrtxt)
	{
		if ((j % 100) == 0)
			std::cout << "Processing the " << j << "th ham query" << endl;

		vector<int64_t> queryVector_v(polyModulus, 0);
		//std::cout << "USER 1 Constructing, encoding and encrypting " << queriesPerCphrtxt << " in one ciphertext ...(BLACK 12)" << endl;

		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = hamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			/*	if (inputFile.fail())
			{
				cout << "could not open file " << filename << endl;
				cin.get();
				cin.ignore();
				exit(0);
			}*/

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * 2 * (m + 1) + m + 1] = 1;

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

					if (word == finalSelectedFeaturesDecrypted[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex + 1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > finalSelectedFeaturesDecrypted[midIndex])
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
				if (word == finalSelectedFeaturesDecrypted[begIndex])
				{
					queryVector_v[i * (m + 1) + begIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == finalSelectedFeaturesDecrypted[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			//vector<int64_t> logHamQueryProb(queriesPerCphrtxt, 0),
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			//for (int k = 0; k <= m; k++)
			//{
			//	queryVector_v[i * 2 * (m + 1) + k] *= trainedModel_v[i * 2 * (m + 1) + k];
			//	logHamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k];
			//	queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[k + m + 1];
			//	logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			//}

			//if (logHamQueryProb[i] > logSpamQueryProb[i])
			//	confusionMatrix[0][0]++;
			//else
			//	confusionMatrix[0][1]++;

			/*countQuery++;
			if ((countQuery % 100) == 0)
			{
				std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
				std::system("pause");
			}*/
		}
		batch_encoder.encode(queryVector_v, queryVector_p);
		encryptor.encrypt(queryVector_p, queryVector_c);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack12 += time_diff.count() / 1000; //long long

		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		//std::system("pause");
		//std::cout << "TACS processing the query (multiply, additions, secComp) ... (BLACK 13)" << endl;

		Ciphertext tmpQueryResult;
		//int step = 0;
		vector <int> powersOfTwo(log2(m + 1), 0);
		for (int i = 0; i < powersOfTwo.size(); i++)
		{
			powersOfTwo[i] = pow(2, i);
		}

		time_start = chrono::high_resolution_clock::now();
		evaluator.multiply_plain_inplace(queryVector_c, trainedModel_p);
		//evaluator.relinearize_inplace(queryVector_c, relin_keys);
		for (int i = 0; i < log2(m + 1); i++)
		{
			//step = powersOfTwo[i];
			tmpQueryResult = queryVector_c;
			evaluator.rotate_rows_inplace(tmpQueryResult, powersOfTwo[i], gal_keys);
			evaluator.add_inplace(queryVector_c, tmpQueryResult);
		}
		//evaluator.rotate_rows(queryVector_c, m + 1, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(queryVector_c, tmpQueryResult);//[A]-[B],
		
		//evaluator.multiply_plain_inplace(queryVector_c, R_p);//([A]-[B])R
		//evaluator.add_plain_inplace(queryVector_c, h_p);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, h_p);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack13 += time_diff.count() / 1000; //long long
		//std::cout << "Here 1A: The noise budget of queryVector_c is " << decryptor.invariant_noise_budget(queryVector_c) << " bits" << endl;
		//std::system("pause");

		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}

		//std::cout << "Here 2: The noise budget of queryVector_c is " << decryptor.invariant_noise_budget(queryVector_c) << " bits" << endl;
		//std::system("pause");
		//std::cout << "EDS decrypting the randomized query ... (BLACK 15)" << endl;
		Plaintext queryResult_p;
		vector<int64_t> queryResult_v;
		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(queryVector_c, queryResult_p);
		batch_encoder.decode(queryResult_p, queryResult_v);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack14 += time_diff.count() / 1000; //long long

		//std::cout << "User 1 getting the final classification of the query(ies) ... (BLACK 15)" << endl;
		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			if (queryResult_v[i * (m + 1)] > 0)
				confusionMatrix[0][0]++;
			else
				confusionMatrix[0][1]++;
		}
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack15 += time_diff.count() / 1000; //long long
	}
	std::cout << endl << "Printing the confusion matrix after PP classification in the end (ham part only)" << endl << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl;

	countQuery = 0;
	string spamPath = path + "\\spam";
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(spamPath))
	{
		filename = entry.path().string();
		spamFiles.push_back(filename);
	}

	for (int j = 0; j < spamFiles.size(); j += queriesPerCphrtxt)
	{
		if ((j % 100) == 0)
			std::cout << "Processing the " << j << "th spam query" << endl;

		vector<int64_t> queryVector_v(polyModulus, 0);
		//std::cout << "USER 1 Constructing, encoding and encrypting " << queriesPerCphrtxt << " in one ciphertext ...(BLACK 12)" << endl;
		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = spamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			/*	if (inputFile.fail())
			{
				cout << "could not open file " << filename << endl;
				cin.get();
				cin.ignore();
				exit(0);
			}*/

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * (m + 1) + m + 1] = 1;

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

					if (word == finalSelectedFeaturesDecrypted[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex + 1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > finalSelectedFeaturesDecrypted[midIndex])
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
				if (word == finalSelectedFeaturesDecrypted[begIndex])
				{
					queryVector_v[i * (m + 1) + begIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == finalSelectedFeaturesDecrypted[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex + 1]++;
					//queryVector_v[i * 2*(m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			//vector<int64_t> logHamQueryProb(queriesPerCphrtxt, 0),
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			//for (int k = 0; k <= m; k++)
			//{
			//	queryVector_v[i * 2 * (m + 1) + k] *= trainedModel_v[i * 2 * (m + 1) + k];
			//	logHamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k];
			//	queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[i * 2 * (m + 1) + k + m + 1];
			//	logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			//}

			//if (logSpamQueryProb[i] > logHamQueryProb[i])
			//	confusionMatrix[1][1]++;
			//else
			//	confusionMatrix[1][0]++;

			/*countQuery++;
			if ((countQuery % 100) == 0)
			{
				std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
				std::system("pause");
			}*/
		}
		batch_encoder.encode(queryVector_v, queryVector_p);
		encryptor.encrypt(queryVector_p, queryVector_c);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack12 += time_diff.count() / 1000; //long long

		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}

		//std::cout << "The noise budget of queryVector_c is " << decryptor.invariant_noise_budget(queryVector_c) << " bits" << endl;
		//std::system("pause");
		//std::cout << "TACS processing the query (multiply, additions, secComp) ... (BLACK 13)" << endl;

		Ciphertext tmpQueryResult;
		//int step = 0;
		time_start = chrono::high_resolution_clock::now();
		evaluator.multiply_plain_inplace(queryVector_c, trainedModel_p);
		//evaluator.relinearize_inplace(queryVector_c, relin_keys);
		for (int i = 0; i < log2(m + 1); i++)
		{
			//step = pow(2, i);
			tmpQueryResult = queryVector_c;
			evaluator.rotate_rows_inplace(tmpQueryResult, powersOfTwo[i], gal_keys);
			evaluator.add_inplace(queryVector_c, tmpQueryResult);
		}
		//evaluator.rotate_rows(queryVector_c, m + 1, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(queryVector_c, tmpQueryResult);//[A]-[B],
		
		//evaluator.multiply_plain_inplace(queryVector_c, R_p);//([A]-[B])R
		//evaluator.add_plain_inplace(queryVector_c, h_p);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, h_p);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack13 += time_diff.count() / 1000; //long long

		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			evaluator.mod_switch_to_next_inplace(queryVector_c);
			//evaluator.mod_switch_to_next_inplace(queryVector_c);
		}
		//std::cout << "The noise budget of queryVector_c is " << decryptor.invariant_noise_budget(queryVector_c) << " bits" << endl;
		//std::system("pause");
		//std::cout << "EDS decrypting the query ...(BLACK 14)" << endl;

		Plaintext queryResult_p;
		vector<int64_t> queryResult_v;

		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(queryVector_c, queryResult_p);
		batch_encoder.decode(queryResult_p, queryResult_v);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack14 += time_diff.count() / 1000; //long long

		//std::cout << "User 1 getting the final classification of the query(ies) ... (BLACK 15)" << endl;
		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			if (queryResult_v[i * (m + 1)] < 0)
				confusionMatrix[1][1]++;
			else
				confusionMatrix[1][0]++;
		}
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack15 += time_diff.count() / 1000; //long long
	}

	std::cout << endl << "Printing the confusion matrix after PP classification in the end" << endl << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl << "The accuracy of the PP classification is: " << (long double)100 * (confusionMatrix[0][0] + confusionMatrix[1][1]) / (long double)(confusionMatrix[0][0] + confusionMatrix[1][1] + confusionMatrix[0][1] + confusionMatrix[1][0]) << " %" << endl << endl;

	std::cout << "totalTimeBlack12=" << totalTimeBlack12 << "ms, timePerQueryBlack12="
		<< (long double)totalTimeBlack12 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 12)" << endl;
	std::cout << "totalTimeBlack13=" << totalTimeBlack13 << "ms, timePerQueryBlack13="
		<< (long double)totalTimeBlack13 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 13)" << endl;
	std::cout << "totalTimeBlack14=" << totalTimeBlack14 << "ms, timePerQueryBlack14="
		<< (long double)totalTimeBlack14 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 14)" << endl;
	std::cout << "totalTimeBlack15=" << totalTimeBlack15 << "ms, timePerQueryBlack15="
		<< (long double)totalTimeBlack15 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 15)" << endl;

	output << "totalTimeBlack12=" << totalTimeBlack12 << "ms, timePerQueryBlack12="
		<< (long double)totalTimeBlack12 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 12)" << endl;
	output << "totalTimeBlack13=" << totalTimeBlack13 << "ms, timePerQueryBlack13="
		<< (long double)totalTimeBlack13 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 13)" << endl;
	output << "totalTimeBlack14=" << totalTimeBlack14 << "ms, timePerQueryBlack14="
		<< (long double)totalTimeBlack14 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 14)" << endl;
	output << "totalTimeBlack15=" << totalTimeBlack15 << "ms, timePerQueryBlack15="
		<< (long double)totalTimeBlack15 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 15)" << endl;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "END OF SERVER CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "END OF SERVER CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	output << endl << "***********************************************************************" << endl;

}

void PPClientCenMNClass(string path, const vector<string>& finalSelectedFeaturesDecrypted,
	vector<int64_t> trainedModel_v, int polyModulus, int m, ofstream &output)
{
	long long int totalTimeBlack16 = 0, totalTimeBlack17 = 0, totalTimeBlack18 = 0, totalTimeBlack19 = 0,
		sizeWhite16 = 0, sizeWhite17 = 0, sizeWhite18 = 0, sizeWhite19 = 0;

	int countQuery = 0, queriesPerCphrtxt = polyModulus / ((m + 1));

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "BEGGINING OF CLIENT CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "BEGGINING OF CLIENT CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	output << endl << "***********************************************************************" << endl;

	vector<vector<int>> confusionMatrix{ {0,0},{0,0} };
	string word;
	string hamPath = path + "\\ham", filename;
	vector<string> hamFiles, spamFiles;

	vector<int64_t> R_v(polyModulus, 0);
	vector<int64_t> h_v(polyModulus, 0);
	vector<int64_t> STCguard_v(polyModulus, 0);
	RandGen rndGen;

	uint64_t plainModulus;//needed after decryption comparison
	if (polyModulus == 4096)
		plainModulus = 65537; //65537, 6144001
	if (polyModulus == 8192)
		plainModulus = 65536524289; //65536524289, 3288203001857,
	if (polyModulus == 16384)
		plainModulus = 65536524289; //65536524289, 3288203001857, 576460752302473217

	for (int i = 0; i < queriesPerCphrtxt; i++)
	{
		R_v[(m + 1) * i] = rndGen.RandInt(1, 4096);
		h_v[(m + 1) * i] = rndGen.RandInt(-R_v[(m + 1) * i], R_v[(m + 1) * i]);

		for (int j = 0; j < m + 1; j++)
		{
			STCguard_v[i * (m + 1) + j] = rndGen.RandInt(-4096, 4096);
			//STCguard_v[i * (m + 1) + j + m + 1] = STCguard_v[i * (m + 1) + j];
		}
	}

	EncryptionParameters parms(scheme_type::BFV);

	parms.set_poly_modulus_degree(polyModulus); //4096, 8192, 16384,
	parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));
	parms.set_plain_modulus(plainModulus);

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "*****PHASE 1: TACS RANDOMIZES THE ENCYPTED TRAINED MODEL*****" << endl;
	std::cout << "******IN ORDER TO PROTECT FROM STC, SENDS IT TO USER 2*******" << endl;
	std::cout << "***************************(TACS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "*****PHASE 1: TACS RANDOMIZES THE ENCYPTED TRAINED MODEL*****" << endl;
	output << "******IN ORDER TO PROTECT FROM STC, SENDS IT TO USER 2*******" << endl;
	output << "***************************(TACS)****************************" << endl;

	auto context = SEALContext::Create(parms);
	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();

	std::cout << "Batching enabled inside the FUNCTION: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(DefaultParams::dbc_max());

	auto relin_keys = keygen.relin_keys(DefaultParams::dbc_max());

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);

	Plaintext trainedModel_p, R_p, h_p, STCguard_p;
	Ciphertext trainedModel_c;

	batch_encoder.encode(trainedModel_v, trainedModel_p);
	encryptor.encrypt(trainedModel_p, trainedModel_c);

	batch_encoder.encode(STCguard_v, STCguard_p);
	batch_encoder.encode(R_v, R_p);
	batch_encoder.encode(h_v, h_p);

	std::cout << "Adding the STCguard to the trained model at TACS ... (BLACK 16)" << endl;
	std::cout << "Adding the STCguard to the trained model at TACS ... (BLACK 16)" << endl;

	auto context_data = context->context_data();


	if (polyModulus == 8192)
	{
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();

		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		//evaluator.mod_switch_to_next_inplace(trainedModel_c);
	}
	else //if(polyModulus == 16384)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();

		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		evaluator.mod_switch_to_next_inplace(trainedModel_c);
		//evaluator.mod_switch_to_next_inplace(trainedModel_c);
	}


	std::cout << "TACS sends the randomized trained model to USER 2 ... (WHITE 16)" << endl;
	output << "TACS sends the randomized trained model to USER 2 ... (WHITE 16)" << endl;

	sizeWhite16 = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()=\nsizeWhite16="
		<< sizeWhite16 << " KB \nto transmit from TACS to USER 2 ... (WHITE 16)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "In all we have 2*coeff_modulus_bit_count()*poly_modulus_degree()=\nsizeWhite16"
		<< sizeWhite16 << " KB \nto transmit from TACS to USER 2 ...(WHITE 16)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---*************************" << endl;
	std::cout << "**********PHASE 2: HOMOMORPHICALLY PROCESSING THE ************" << endl;
	std::cout << "*********THE QUERY AT USER 2 AND SENDING IT TO EDS************" << endl;
	std::cout << "***************************(USER 2)***************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---*************************" << endl;
	output << "**********PHASE 2: HOMOMORPHICALLY PROCESSING THE ************" << endl;
	output << "*********THE QUERY AT USER 2 AND SENDING IT TO EDS************" << endl;
	output << "***************************(USER 2)***************************" << endl;

	std::cout << "USER 2 homomorphically processes the query and randomizes it ... (BLACK 17)" << endl;
	output << "USER 2 homomorphically processes the query and randomizes it ... (BLACK 17)" << endl;

	std::cout << "USER 2 sends his processed randomized data to EDS ... (WHITE 17)" << endl;
	output << "USER 2 sends his processed randomized data to EDS ... (WHITE 17)" << endl;


	Ciphertext tmpTV_c = trainedModel_c;
	if (polyModulus == 8192)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
	}
	else //if(polyModulus == 16384)
	{
		context_data = context_data->next_context_data();
		context_data = context_data->next_context_data();
		//context_data = context_data->next_context_data();
	}


	sizeWhite17 = (2 * context_data->total_coeff_modulus_bit_count()*(int)context_data->parms().poly_modulus_degree()) / (8 * 1024);
	std::cout << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	std::cout << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	std::cout << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()=\nsizeWhite17 "
		<< sizeWhite17 << " KB \nto transmit from USER 2 to EDS ...(WHITE 17)" << endl;

	output << "The size of the coefficient modulus is: "
		<< context_data->total_coeff_modulus_bit_count() << " bits" << endl;
	output << "The size (degree) of the polynomial is: "
		<< context_data->parms().poly_modulus_degree() << endl;
	output << "We have 2*coeff_modulus_bit_count()*poly_modulus_degree()=\nsizeWhite17 "
		<< sizeWhite17 << " KB \nto transmit from USER 2 to EDS ...(WHITE 17)" << endl;

	std::cout << endl << endl;
	std::cout << "********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	std::cout << "****************AND SENDING IT TO THE USER 2*****************" << endl;
	std::cout << "****************************(EDS)****************************" << endl;
	output << endl << endl;
	output << "********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 3: DECRYPTING AND DECODING THE RANDOMIZED QUERY****" << endl;
	output << "****************AND SENDING IT TO THE USER 2*****************" << endl;
	output << "****************************(EDS)****************************" << endl;

	std::cout << "EDS decrypts the final randomized classification ... (BLACK 18)" << endl;
	output << "EDS decrypts the final randomized classification ... (BLACK 18)" << endl;

	std::cout << "EDS sends the randomized result in plain to USER 2 ... (WHITE 18)" << endl;
	output << "EDS sends the randomized result in plain to USER 2 ... (WHITE 18)" << endl;

	std::cout << "We havequeriesPerCphrtxt*log2(Nr_classes)=" << queriesPerCphrtxt << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 2, which is negligable ... (WHITE 18)" << endl;
	output << "We havequeriesPerCphrtxt*log2(Nr_classes)=" << queriesPerCphrtxt << " bits transmitted (1 bit per query in our case) "
		<< "\nfrom EDS to USER 2, which is negligable ... (WHITE 18)" << endl;

	int totalKBClassification = sizeWhite16 + sizeWhite17;

	std::cout << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the user centric classification" << endl
		<< "with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries per batch" << endl;

	output << "In total we have " << totalKBClassification << " KB (KBytes) transmitted for the user centric classification" << endl
		<< "with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries per batch" << endl;

	std::cout << endl << endl;
	std::cout << "*********************---CLASSIFYING---************************" << endl;
	std::cout << "****PHASE 4: THE USER 2 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	std::cout << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	std::cout << "****************************(USER 2)**************************" << endl;
	output << endl << endl;
	output << "*********************---CLASSIFYING---************************" << endl;
	output << "****PHASE 4: THE USER 2 DE-RANDOMIZES THE RANDOMIZED QUERY****" << endl;
	output << "************RESULT TO GET THE FINAL CLASSIFICATION************" << endl;
	output << "****************************(USER 2)**************************" << endl;

	std::cout << "USER 1 de-randomizes the result to get the final classification ... (BLACK 19) " << endl;
	output << "USER 1 de-randomizes the result to get the final classification ... (BLACK 19) " << endl;

	std::cout << "For the User Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query of " << queriesPerCphrtxt << " batched queries" << endl;
	output << "For the User Centric classification a total of \n" << totalKBClassification
		<< " KB of data are transmmited per query ofqueriesPerCphrtxt=" << queriesPerCphrtxt << " batched queries" << endl;

	//ACTUAL PP SERVER CENTRIC CLASSIFICATION
	chrono::high_resolution_clock::time_point time_start, time_end;
	time_start = chrono::high_resolution_clock::now();
	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

	std::cout << "Adding the STC guard to the trained model" << endl;
	Ciphertext tmpSTCTrainedModel = trainedModel_c, tmpQueryResult, queryVector_c;

	Plaintext queryVector_p;

	vector <int> powersOfTwo(log2(m + 1), 0);
	for (int i = 0; i < powersOfTwo.size(); i++)
	{
		powersOfTwo[i] = pow(2, i);
	}

	for (const auto & entry : fs::filesystem::recursive_directory_iterator(hamPath))
	{
		filename = entry.path().string();
		hamFiles.push_back(filename);
	}

	for (int j = 0; j < hamFiles.size(); j += queriesPerCphrtxt)
	{
		if ((j % 100) == 0)
			std::cout << "Processing the " << j << "th ham query" << endl;

		//std::cout << "Adding the STC guard to the trained model and sending it to USER 2 ... (BLACK 16)" << endl;
		tmpSTCTrainedModel = trainedModel_c;
		time_start = chrono::high_resolution_clock::now();
		//evaluator.add_plain_inplace(tmpSTCTrainedModel, STCguard_p);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack16 += time_diff.count() / 1000; //long long

		//std::cout << "HERE 1: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");
		//std::cout << "User 2 constructs, encodes the queryVector and multiplies it with the tmpSTCTrainedModel" << endl;

		vector<int64_t> queryVector_v(polyModulus, 0);
		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = hamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			/*	if (inputFile.fail())
			{
				cout << "could not open file " << filename << endl;
				cin.get();
				cin.ignore();
				exit(0);
			}*/

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * 2 * (m + 1) + m + 1] = 1;

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

					if (word == finalSelectedFeaturesDecrypted[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex + 1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > finalSelectedFeaturesDecrypted[midIndex])
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
				if (word == finalSelectedFeaturesDecrypted[begIndex])
				{
					queryVector_v[i * (m + 1) + begIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == finalSelectedFeaturesDecrypted[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			//vector<int64_t> logHamQueryProb(queriesPerCphrtxt, 0),
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			//for (int k = 0; k <= m; k++)
			//{
			//	queryVector_v[i * 2 * (m + 1) + k] *= trainedModel_v[i * 2 * (m + 1) + k];
			//	logHamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k];
			//	queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[k + m + 1];
			//	logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			//}

			//if (logHamQueryProb[i] > logSpamQueryProb[i])
			//	confusionMatrix[0][0]++;
			//else
			//	confusionMatrix[0][1]++;

			/*countQuery++;
			if ((countQuery % 100) == 0)
			{
				std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
				std::system("pause");
			}*/
		}
		batch_encoder.encode(queryVector_v, queryVector_p);
		//encryptor.encrypt(queryVector_p, queryVector_c);

		evaluator.multiply_plain_inplace(tmpSTCTrainedModel, queryVector_p);
		//evaluator.relinearize_inplace(tmpSTCTrainedModel, relin_keys);
		for (int i = 0; i < log2(m + 1); i++)
		{
			tmpQueryResult = tmpSTCTrainedModel;
			evaluator.rotate_rows_inplace(tmpQueryResult, powersOfTwo[i], gal_keys);
			evaluator.add_inplace(tmpSTCTrainedModel, tmpQueryResult);
		}
		//evaluator.rotate_rows(tmpSTCTrainedModel, m + 1, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(tmpSTCTrainedModel, tmpQueryResult);//[A]-[B],
		
		//evaluator.multiply_plain_inplace(tmpSTCTrainedModel, R_p);//([A]-[B])R
		//evaluator.add_plain_inplace(tmpSTCTrainedModel, h_p);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, h_p);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack17 += time_diff.count() / 1000; //long long

		//std::cout << "HERE 2: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");


		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			//evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
		}

		//std::cout << "HERE 3: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");
		//std::cout << "EDS decrypting the query ... (BLACK 18)" << endl;
		Plaintext queryResult_p;
		vector<int64_t> queryResult_v;

		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(tmpSTCTrainedModel, queryResult_p);
		batch_encoder.decode(queryResult_p, queryResult_v);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack18 += time_diff.count() / 1000; //long long

		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			if (queryResult_v[i * (m + 1)] > 0)
				confusionMatrix[0][0]++;
			else
				confusionMatrix[0][1]++;
		}
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack19 += time_diff.count() / 1000; //long long
	}
	std::cout << endl << "Printing the confusion matrix after PP classification in the end (ham part only)" << endl << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl;

	countQuery = 0;
	string spamPath = path + "\\spam";
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(spamPath))
	{
		filename = entry.path().string();
		spamFiles.push_back(filename);
	}

	for (int j = 0; j < spamFiles.size(); j += queriesPerCphrtxt)
	{
		if ((j % 100) == 0)
			std::cout << "Processing the " << j << "th spam query" << endl;

		//std::cout << "Adding the STC guard to the trained model and sending it to USER 2 ... (BLACK 16)" << endl;
		tmpSTCTrainedModel = trainedModel_c;
		time_start = chrono::high_resolution_clock::now();
		//evaluator.add_plain_inplace(tmpSTCTrainedModel, STCguard_p);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack16 += time_diff.count() / 1000; //long long

		//std::cout << "HERE 1: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");
		//std::cout << "User 2 constructs, encodes the queryVector and multiplies it with the tmpSTCTrainedModel" << endl;

		vector<int64_t> queryVector_v(polyModulus, 0);
		time_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = spamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			/*	if (inputFile.fail())
			{
				cout << "could not open file " << filename << endl;
				cin.get();
				cin.ignore();
				exit(0);
			}*/

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * 2 * (m + 1) + m + 1] = 1;

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

					if (word == finalSelectedFeaturesDecrypted[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex + 1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > finalSelectedFeaturesDecrypted[midIndex])
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
				if (word == finalSelectedFeaturesDecrypted[begIndex])
				{
					queryVector_v[i * (m + 1) + begIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == finalSelectedFeaturesDecrypted[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			//vector<int64_t> logHamQueryProb(queriesPerCphrtxt, 0),
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			//for (int k = 0; k <= m; k++)
			//{
			//	queryVector_v[i * 2 * (m + 1) + k] *= trainedModel_v[i * 2 * (m + 1) + k];
			//	logHamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k];
			//	queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[k + m + 1];
			//	logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			//}

			//if (logHamQueryProb[i] > logSpamQueryProb[i])
			//	confusionMatrix[0][0]++;
			//else
			//	confusionMatrix[0][1]++;

			/*countQuery++;
			if ((countQuery % 100) == 0)
			{
				std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
				std::system("pause");
			}*/
		}
		batch_encoder.encode(queryVector_v, queryVector_p);
		//encryptor.encrypt(queryVector_p, queryVector_c);

		evaluator.multiply_plain_inplace(tmpSTCTrainedModel, queryVector_p);
		//evaluator.relinearize_inplace(tmpSTCTrainedModel, relin_keys);
		for (int i = 0; i < log2(m + 1); i++)
		{
			tmpQueryResult = tmpSTCTrainedModel;
			evaluator.rotate_rows_inplace(tmpQueryResult, powersOfTwo[i], gal_keys);
			evaluator.add_inplace(tmpSTCTrainedModel, tmpQueryResult);
		}
		//evaluator.rotate_rows(tmpSTCTrainedModel, m + 1, gal_keys, tmpQueryResult);
		//evaluator.sub_inplace(tmpSTCTrainedModel, tmpQueryResult);//[A]-[B],
		
		//evaluator.multiply_plain_inplace(tmpSTCTrainedModel, R_p);//([A]-[B])R
		//evaluator.add_plain_inplace(tmpSTCTrainedModel, h_p);//([A]-[B])R + H
		//evaluator.sub_plain_inplace(queryVectorCipher, h_p);//([A]-[B])R + H
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack17 += time_diff.count() / 1000; //long long

		//std::cout << "HERE 2: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");


		if (polyModulus == 8192)
		{
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
		}
		else //if(polyModulus == 16384)
		{
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
			//evaluator.mod_switch_to_next_inplace(tmpSTCTrainedModel);
		}

		//std::cout << "HERE 3: the noise budget of tmpSTCTrainedModel is "
		//	<< decryptor.invariant_noise_budget(tmpSTCTrainedModel) << " bits" << endl;
		//std::system("pause");
		//std::cout << "EDS decrypting the query ... (BLACK 18)" << endl;
		Plaintext queryResult_p;
		vector<int64_t> queryResult_v;

		time_start = chrono::high_resolution_clock::now();
		decryptor.decrypt(tmpSTCTrainedModel, queryResult_p);
		batch_encoder.decode(queryResult_p, queryResult_v);
		time_end = chrono::high_resolution_clock::now();
		time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		totalTimeBlack18 += time_diff.count() / 1000; //long long

		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			if (queryResult_v[i * (m + 1)] < 0)
				confusionMatrix[1][1]++;
			else
				confusionMatrix[1][0]++;
		}
	}

	std::cout << endl << "Printing the confusion matrix after PP classification in the end" << endl << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl << "The accuracy of the PP classification is: " << (long double)100 * (confusionMatrix[0][0] + confusionMatrix[1][1]) / (long double)(confusionMatrix[0][0] + confusionMatrix[1][1] + confusionMatrix[0][1] + confusionMatrix[1][0]) << " %" << endl << endl;

	std::cout << "totalTimeBlack16=" << totalTimeBlack16 << "ms, timePerQueryBlack16="
		<< (long double)totalTimeBlack16 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 12)" << endl;
	std::cout << "totalTimeBlack17=" << totalTimeBlack17 << "ms, timePerQueryBlack17="
		<< (long double)totalTimeBlack17 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 13)" << endl;
	std::cout << "totalTimeBlack18=" << totalTimeBlack18 << "ms, timePerQueryBlack18="
		<< (long double)totalTimeBlack18 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 14)" << endl;
	std::cout << "totalTimeBlack19=" << totalTimeBlack19 << "ms, timePerQueryBlack15="
		<< (long double)totalTimeBlack19 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 15)" << endl;

	output << "totalTimeBlack16=" << totalTimeBlack16 << "ms, timePerQueryBlack16="
		<< (long double)totalTimeBlack16 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 12)" << endl;
	output << "totalTimeBlack17=" << totalTimeBlack17 << "ms, timePerQueryBlack17="
		<< (long double)totalTimeBlack17 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 13)" << endl;
	output << "totalTimeBlack18=" << totalTimeBlack18 << "ms, timePerQueryBlack18="
		<< (long double)totalTimeBlack18 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 14)" << endl;
	output << "totalTimeBlack19=" << totalTimeBlack19 << "ms, timePerQueryBlack19="
		<< (long double)totalTimeBlack19 / (hamFiles.size() + spamFiles.size()) << " ms per one query ...(BLACK 15)" << endl;

	std::cout << endl << "***********************************************************************" << endl;
	std::cout << "END OF CLIENT CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	std::cout << endl << "***********************************************************************" << endl;
	output << endl << "***********************************************************************" << endl;
	output << "END OF CLIENT CENTRIC Processing of queries with queriesPerCphrtxt=" << queriesPerCphrtxt << " queries" << endl;
	output << endl << "***********************************************************************" << endl;
}


int main()
{
#ifdef SEAL_VERSION
    cout << "Microsoft SEAL version: " << SEAL_VERSION << endl;
#endif
	
	EncryptionParameters parms(scheme_type::BFV);
	// EncryptionParameters parms(scheme_type::CKKS);

	int polyModulus = 16384; //4096, 8192, 16384, 32768, put here the value for the polynomial Modulus
	//int polyModulus = 32768;
	int untilForRndNumber = 200, dbc = 40;
	uint64_t tmpPlainModulus = 40961;

	std::cout << "Enter the PolyModulus value (4096 or 8192 or 16384):" << endl;
	//cin >> tmpolyModulus;
	//polyModulus = 16384;
	while ((polyModulus != 8192) && (polyModulus != 16384) && (polyModulus != 4096))
	{
		std::cout << "Enter the PolyModulus value again (4096 or 8192 or 16384):" << endl;
		cin >> polyModulus;
	}

	if (polyModulus == 4096)
	{
		tmpPlainModulus = 40961;
		untilForRndNumber = 200;
		dbc = 60;
	}
	else if (polyModulus == 8192)
	{
		tmpPlainModulus = 65537; //65536524289 -> 37 bits; 65537 -> 15 bits, 65537
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}
	else
	{
		tmpPlainModulus = 65537; //60 bits 576460752304439297, 557057 -> 20 bits
		untilForRndNumber = pow(2, 10);
		dbc = 60;
	}

	parms.set_poly_modulus_degree(polyModulus);//4096, 8192, 16384, 32768 
	parms.set_coeff_modulus(DefaultParams::coeff_modulus_128(polyModulus));//4096, 8192, 16384, 32768 

	parms.set_plain_modulus(tmpPlainModulus);//40961 for 4096, and 65536524289 for 8192, 576460752304439297 for others

	auto context = SEALContext::Create(parms);

	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();
	std::cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(dbc);

	auto relin_keys = keygen.relin_keys(dbc);

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);

	int fileNr = 2;
	string fileNameOutput = "secMLClass,FileNr=" + itoa(fileNr)
		+ ",polyModulusSWITCH=" + itoa(polyModulus) + ".txt";
	ofstream output(fileNameOutput);
	
	std::cout << endl << endl << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "A) DEALING WITH THE WISCONSIN BREAST CANCER DATASET NOW" << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "********************************************************" << endl << endl;

	std::cout << "*******************************************" << endl;
	std::cout << "The NAIVE BAYES Classification - PLAIN CASE" << endl;
	std::cout << "*******************************************" << endl;
	
	// NaiveBayesBreastCancer nbBreast("BreastCancer.txt");

	// clock_t nb_start, nb_finish;
	// double duration;
	// nb_start = clock();

	NaiveBayesBreastCancer nbBreast("/home/byhan/heaan/heaan_stats-dev-1/heaan_stats/compare-papers/SEAL-master/native/examples/BreastCancer.txt");

	nbBreast.findSystemCounts(1);
	nbBreast.findProbsAndLogsFromCounts(0);
	nbBreast.accuracyTest();

	// nb_finish = clock();
	// duration = nb_finish - nb_start;
    // cout << duration << "ms" << endl;

	nbBreast.printConfusionMatrix();
	vector<vector<int>> localDataset;
	nbBreast.getDatasetMatrix(localDataset);

	vector<vector<long double>> globalOutputProbNonCancerF, globalOutputProbCancerF;
	vector<int64_t> globalTrainingModelVector; globalTrainingModelVector.resize(polyModulus, 0);
	long double globalProbNonCancer, globalProbCancer;
	nbBreast.getLogsOfProbs(globalOutputProbCancerF, globalOutputProbNonCancerF, globalProbCancer, globalProbNonCancer);

	int index = 0, constK=127;
	globalTrainingModelVector[index++] = (constK)*(globalProbNonCancer- globalProbCancer);
	for (int i = 0; i < globalOutputProbNonCancerF.size(); i++)
	{
		for (int j = 0; j < globalOutputProbNonCancerF[0].size(); j++)
		{
			globalTrainingModelVector[index++] =
				(constK)*(globalOutputProbNonCancerF[i][j]- globalOutputProbCancerF[i][j]);
		}
	}
	
	//index = 128;
	//globalTrainingModelVector[index++] = (-constK)*globalProbCancer;
	//for (int i = 0; i < globalOutputProbCancerF.size(); i++)
	//{
	//	for (int j = 0; j < globalOutputProbCancerF[0].size(); j++)
	//	{
	//		globalTrainingModelVector[index++] = (-constK)*globalOutputProbCancerF[i][j];
	//	}
	//}

	//for (int i = 0; i < 256; i++)
	//{
		//std::cout << globalTrainingModelVector[i] << endl;
	//}

	// 아... 딱 저 ServerCentric, ClientCentric이 암호화부분이네..
	int batchQuerySize = polyModulus / 128;
	cout << endl << endl << "------------------------------------------" << endl;
	std::cout << "polyModulus=" << polyModulus << ", batchQuerySize=" << batchQuerySize << endl;
	cout << "------------------------------------------" << endl;
	//output << endl << endl << "------------------------------------------" << endl;
	//output << "polyModulus=" << polyModulus << ", batchQuerySize=" << batchQuerySize << endl;
	//output << "------------------------------------------" << endl;
	
	NB_BreastC_ServerCentric(polyModulus,tmpPlainModulus, context, gal_keys, decryptor,
		batch_encoder,encryptor,evaluator, globalTrainingModelVector, batchQuerySize, localDataset, output);
	NB_BreastC_ClientCentric(polyModulus, tmpPlainModulus, context, gal_keys, decryptor,
		batch_encoder, encryptor, evaluator, globalTrainingModelVector, batchQuerySize, localDataset, output);
	//std::system("pause");
	//딱 여기까지만 보면 됌
	// system("pause");	
	// return 0;

	std::cout << endl;
	std::cout << "**********************************************" << endl;
	std::cout << "The SVM and LR ML Classifications - PLAIN CASE" << endl;
	std::cout << "**********************************************" << endl;

	vector<double> tmp_weights(10);
	vector<double> tmp_weights2(10);
	vector<int64_t> svm_weights(polyModulus,0);
	vector<int64_t> logistic_weights(polyModulus,0);

	vector<vector<int64_t>> dataBCancer(ceil((double)(699*16)/polyModulus), vector<int64_t> (polyModulus, 0));
	//std::cout << dataBCancer.size()<<endl<< dataBCancer[0].size()<<endl;
	
	vector<int> true_label(699);
	vector<int64_t> res_svm(699);
	vector<int64_t> res_logistic(699);
	vector<int> predicted_label_svm(699);
	vector<int> predicted_label_logistic(699);
	// SVM weights found by WEKA, weights multiplied with 100 and rounded
	tmp_weights[0] = -434; //-4.3373;
	tmp_weights[1] = 25; //0.2506;
	tmp_weights[2] = 0;  //-0.0045;
	tmp_weights[3] = 19; //0.1854;
	tmp_weights[4] = 7;  //0.0668;
	tmp_weights[5] = 5;  //0.047;
	tmp_weights[6] = 10;//.1977;
	tmp_weights[7] = 19; //0.1875;
	tmp_weights[8] = 8;  //0.0837;
	tmp_weights[9] = 18; //0.184;
	// Logistic Regression weights found by WEKA, weights multiplied with 100 and rounded
	tmp_weights2[0] = -967;//-9.6728;
	tmp_weights2[1] = 53;//0.5313;
	tmp_weights2[2] = 1;//0.0069;
	tmp_weights2[3] = 33;//0.3301;
	tmp_weights2[4] = 24;//0.2393;
	tmp_weights2[5] = 7;//0.0676;
	tmp_weights2[6] = 41;//0.4068;
	tmp_weights2[7] = 41;//0.4093;
	tmp_weights2[8] = 15;//0.1463;
	tmp_weights2[9] = 55;//0.5488;
	// Weights of SVM and Logistic Regression stacked to vector with size 16384
	int count1 = 0;
	int count2 = 0;
	while (count1 < polyModulus)
	{
		svm_weights[count1] = 1;
		count2 = 0;
		while (count2 < 10)
		{
			svm_weights[count1] = tmp_weights[count2];
			logistic_weights[count1] = tmp_weights2[count2];
			count1 += 1;
			count2 += 1;
		}
		count2 = 0;
		//Stack 6 zeros to make the weights power of 2 (i.e. 16)
		while (count2 < 6)
		{
			svm_weights[count1] = 0;
			logistic_weights[count1] = 0;
			count1 += 1;
			count2 += 1;
		}
	}

	// Data Reading
	// string filename = "BreastCancer.txt"; // Name of the data txt file
	string filename = "/home/byhan/heaan/heaan_stats-dev-1/heaan_stats/compare-papers/SEAL-master/native/examples/BreastCancer.txt";
	string line = "";
	
	ifstream file(filename);

	if (file.fail())
	{
		cout << "could not open file " << filename << endl;
		cin.get();
		cin.ignore();
		return 0;
	}
	
	count1 = 0;
	count2 = 0;
	int count_tt = 0;
	int data_counter = -1;
	double rate = 0;
	char trash;
	while (getline(file, line)) // Taking line from data txt file
	{
		istringstream strfile2(line);
		dataBCancer[count1/polyModulus][count1%polyModulus] = 1;
		count2 = 0;
		count1 += 1;
		data_counter += 1;
		while (strfile2 >> rate)
		{
			if (count2 != 9)
			{
				dataBCancer[count1 / polyModulus][count1%polyModulus] = rate;
				count2 += 1;
				count1 += 1;
				//cout << rate << " ";
			}
			else
			{
				true_label[data_counter] = rate;
				count2 += 1;
				//cout << "\n";
			}

			if (count2 != 10)
			{
				strfile2 >> trash;
			}
		}
		//Stack 6 zeros to make the features power of 2 (i.e. 16)
		for (int i = 0; i < 6; i++)
		{
			dataBCancer[count1 / polyModulus][count1%polyModulus] = 0;
			count1 += 1;
		}

	}
	// Results Multiplication
	count1 = 0;
	double maxValue = 0;
	double tmp_svm = 0;
	double tmp_logistic = 0;
	for (int i = 0; i < 699; i++)
	{
		count2 = 0;
		tmp_svm = 0;
		tmp_logistic = 0;
		while (count2 < 16)
		{
			tmp_svm += svm_weights[count1%polyModulus] * dataBCancer[count1 / polyModulus][count1%polyModulus];
			tmp_logistic += logistic_weights[count1%polyModulus] * dataBCancer[count1/polyModulus][count1%polyModulus];
			count2 += 1;
			count1 += 1;
		}

		if (maxValue < tmp_svm)
		{
			maxValue = tmp_svm;
		}
		if (maxValue < tmp_logistic)
		{
			maxValue = tmp_logistic;
		}

		res_svm[i] = tmp_svm;
		res_logistic[i] = tmp_logistic;
		//cout << tmp_logistic << "\n";
		// SVM final prediction
		if (tmp_svm >= 0)
		{
			predicted_label_svm[i] = 2;
		}
		else
		{
			predicted_label_svm[i] = 1;
		}
		// Logistic regression final prediction
		if (tmp_logistic >= 0)
		{
			predicted_label_logistic[i] = 2;
		}
		else
		{
			predicted_label_logistic[i] = 1;
		}

	}

	int SVM_accuracy = 0, LR_accuracy = 0;
	for (int i = 0; i < 699; i++)
	{
		if (true_label[i] == predicted_label_svm[i])
		{
			SVM_accuracy++;
			//std::cout << "true_label[i]=" << true_label[i] << endl;
			//std::cout << "predicted_label_svm[i]=" << predicted_label_svm[i] << endl;
		}
			
		if (true_label[i] == predicted_label_logistic[i])
		{
			LR_accuracy++;
			//std::cout << "true_label[i]=" << true_label[i] << endl;
			//std::cout << "predicted_label_logistic[i]=" << predicted_label_logistic[i] << endl;
		}
	}
	//std::cout << SVM_accuracy <<", "<< LR_accuracy << endl;
	std::cout << "The SVM accuracy is: " << (double)(SVM_accuracy * 100) / 699 << endl;
	std::cout << "The LR accuracy is: " << (double)(LR_accuracy * 100) / 699 << endl;
	//std::cout << "maxValue=" << maxValue << endl;
	
	
	std::cout << endl;
	std::cout << "************************************************************************" << endl;
	std::cout << "The SVM and LR secMLClass (secure ML Classification)-CLIENT CENTRIC CASE" << endl;
	std::cout << "************************************************************************" << endl;

	std::cout << "Calling the client centric function for both SVM and LR" << endl;
	//system("pause");
	SVM_LR_ClientCentric(polyModulus, svm_weights, logistic_weights, dataBCancer, gal_keys, decryptor, true_label,
		batch_encoder, encryptor, evaluator, 699, 2);
	
	std::cout << endl;
	std::cout << "************************************************************************" << endl;
	std::cout << "The SVM and LR secMLClass (secure ML Classification)-SERVER CENTRIC CASE" << endl;
	std::cout << "************************************************************************" << endl;

	std::cout << "Calling the server centric function for both SVM and LR" << endl;
	//system("pause");
	SVM_LR_ServerCentric(polyModulus, svm_weights, logistic_weights, dataBCancer, gal_keys, decryptor, true_label,
		batch_encoder, encryptor, evaluator, 699, 2);
	   	  

	std::cout << endl;
	std::cout << "**********************************************" << endl;
	std::cout << "The Deep Learning Classifications - PLAIN CASE" << endl;
	std::cout << "**********************************************" << endl;

	vector<vector<int>> DL_weights1(16, vector<int>(16, 0));
	vector<int> DL_weights2(16);
	//vector<vector<int>> mlp_weights1(16384, vector<int>(16));
	//vector<int> mlp_weights2(16384);
	
	vector<double> DLdata(16384, 0);
	vector<int> DL_true_label(699);

	vector<double> res_mlp(699);
	vector<int> predicted_label_mlp(699);

	// MLP last layer weights found by Keras Library Python, weights multiplied with 1000 and rounded
	/*DL_weights2[0] = 668;
	//DL_weights2[1] = 54;
	//DL_weights2[2] = -281;
	//DL_weights2[3] = -382;
	//DL_weights2[4] = 571;
	//DL_weights2[5] = -874;
	//DL_weights2[6] = 238;
	//DL_weights2[7] = 448;
	//DL_weights2[8] = -401;
	//DL_weights2[9] = 433;
	//DL_weights2[10] = 69;
	//DL_weights2[11] = -736;
	//DL_weights2[12] = 517;
	//DL_weights2[13] = -8;
	//DL_weights2[14] = -117;
	//DL_weights2[15] = -306;*/
	
	DL_weights2[0] = 66;
	DL_weights2[1] = 5;
	DL_weights2[2] = -28;
	DL_weights2[3] = -38;
	DL_weights2[4] = 57;
	DL_weights2[5] = -87;
	DL_weights2[6] = 23;
	DL_weights2[7] = 44;
	DL_weights2[8] = -40;
	DL_weights2[9] = 43;
	DL_weights2[10] = 6;
	DL_weights2[11] = -73;
	DL_weights2[12] = 51;
	DL_weights2[13] = 0;
	DL_weights2[14] = -11;
	DL_weights2[15] = -30;
	
	// Data Reading
	//int count1, count2;
	ifstream fileWeights1;
	string filenameWeights1 = "weights1.txt"; // Name of the data txt file
	line = "";
	fileWeights1.open(filenameWeights1.c_str());
	count1 = 0;
	count2 = 0;
	int rate1 = 0;
	//trash;
	while (getline(fileWeights1, line)) // Taking line from data txt file
	{
		istringstream strfile2(line);
		count2 = 0;
		while (strfile2 >> rate1)
		{
			DL_weights1[count1][count2] = rate1/10;
			count2 += 1;

			if (count2 != 16)
			{
				strfile2 >> trash;
			}
		}
		count1 += 1;
	}
	ifstream file2;
	// filename = "BreastCancer.txt"; // Name of the data txt file
	filename = "/home/byhan/heaan/heaan_stats-dev-1/heaan_stats/compare-papers/SEAL-master/native/examples/BreastCancer.txt";
	line = "";
	file2.open(filename.c_str());
	count1 = 0;
	count2 = 0;
	count_tt = 0;
	data_counter = -1;
	rate = 0;
	//char trash;
	while (getline(file2, line)) // Taking line from data txt file
	{
		istringstream strfile2(line);
		DLdata[count1] = 1;
		count2 = 0;
		count1 += 1;
		data_counter += 1;
		while (strfile2 >> rate)
		{
			if (count2 != 9)
			{
				DLdata[count1] = rate;
				count2 += 1;
				count1 += 1;
				//cout << rate << " ";
			}
			else
			{
				DL_true_label[data_counter] = rate;
				count2 += 1;
				//cout << "\n";
			}

			if (count2 != 10)
			{
				strfile2 >> trash;
			}
		}
		//Stack 6 zeros to make the features power of 2 (i.e. 16)
		for (int i = 0; i < 6; i++)
		{
			DLdata[count1] = 0;
			count1 += 1;
		}
	}

	long double maxDLRes = 0;
	// Results Multiplication
	count1 = 0;
	vector<double> tmp_mlp(16, 0);
	vector<double> initial(16, 0);
	double tmp_final_mlp = 0;
	//vector<double> data(16384, 0);
	for (int i = 0; i < 699; i++)
	{
		count2 = 0;
		tmp_final_mlp = 0;
		tmp_mlp = initial;
		//tmp_mlp = 0;		
		while (count2 < 16)
		{
			for (int k = 0; k < 16; k++)
			{
				tmp_mlp[k] += DLdata[count1] * DL_weights1[count2][k];
			}
			if (abs(tmp_mlp[count2]) > maxDLRes)
				maxDLRes = abs(tmp_mlp[count2]);
			count2 += 1;
			count1 += 1;
		}
		//std::cout << "For i = " << i << endl;
		for (int k = 0; k < 16; k++)
		//	std::cout << "tmp_mlp[" << k << "]=" << tmp_mlp[k] << endl;
		for (int j = 0; j < 16; j++)
		{
			tmp_final_mlp += (tmp_mlp[j] * tmp_mlp[j])* DL_weights2[j]; //(tmp_mlp[j]* tmp_mlp[j] ) part is square activation

		}
		
		res_mlp[i] = tmp_final_mlp;
		//std::cout << "res_mlp["<<i<<"]="<< res_mlp[i] << endl;
		if (abs(res_mlp[i]) > maxDLRes)
			maxDLRes = abs(tmp_mlp[count2]);
		//res_logistic[i] = tmp_logistic;
		//cout << res_mlp[i] << "\n";
		// SVM final prediction
		if (tmp_final_mlp >= 0)
		{
			predicted_label_mlp[i] = 2;
		}
		else
		{
			predicted_label_mlp[i] = 1;
		}
		//system("pause");
		// Logistic regression final prediction
	}
	double acc = 0;
	for (int i = 0; i < 699; i++)
	{
		if (predicted_label_mlp[i] == DL_true_label[i])
		{
			acc++;
		}
	}
	cout << "Accuracy: " << acc*100 / 699<<"%"<<endl;
	cout << "maxDLRes=" << maxDLRes << endl;

	
	std::cout << endl;
	std::cout << "******************************************************************************" << endl;
	std::cout << "The Deep Learning Classifications - secMLClass (secure ML Classification) CASE" << endl;
	std::cout << "******************************************************************************" << endl;
	
	std::cout << endl;
	std::cout << "********************************************************************" << endl;
	std::cout << "The Deep Learning Classifications - secMLClass - SERVER CENTRIC CASE" << endl;
	std::cout << "********************************************************************" << endl;
		
	DeepLearning_BreastC_ServerCentric(polyModulus, DLdata, DL_weights1, DL_weights2, DL_true_label);
	

	std::cout << endl;
	std::cout << "********************************************************************" << endl;
	std::cout << "The Deep Learning Classifications - secMLClass - CLIENT CENTRIC CASE" << endl;
	std::cout << "********************************************************************" << endl;
		
	DeepLearning_BreastC_ClientCentric(polyModulus, DLdata, DL_weights1, DL_weights2, DL_true_label);
	


	std::cout << endl << endl << endl;
	std::cout << "********************************************" << endl;
	std::cout << "********************************************" << endl;
	std::cout << "B) DEALING WITH THE E-MAIL ENRON DATASET NOW" << endl;
	std::cout << "********************************************" << endl;
	std::cout << "********************************************" << endl << endl;

	std::cout << endl;
	std::cout << "*******************************************************" << endl;
	std::cout << "The MULTINOMIAL NAIVE BAYES Classification - PLAIN CASE" << endl;
	std::cout << "*******************************************************" << endl;
	
	string path = "enron5";
	int m = 2047;
	long long int timeElpased = 0;
	MultinomialNB_Email mNEmail(path, m, timeElpased), tmpmNEmail;

	tmpmNEmail = mNEmail;
	tmpmNEmail.train(timeElpased);
	tmpmNEmail.classify(path, timeElpased);
	tmpmNEmail.printConfionMat();

	int64_t K2 = 8;
	vector<string> selectedFeaturesOutput;
	vector<int64_t> TVGlobal_v(polyModulus, 0);
	vector<int64_t> trainedModel_v(polyModulus, 0),
		trainedModelGlobal_v(polyModulus, 0);

	tmpmNEmail.getSelectedFeaturesOnly(selectedFeaturesOutput);
	tmpmNEmail.getTVaccordingToSelectedFeatures(selectedFeaturesOutput, TVGlobal_v);
	TVGlobal_v[0] += m;
	TVGlobal_v[m + 1] += m;
	TVGlobal_v[TVGlobal_v.size() - 1] += 2 * m;
	for (int i = 1; i <= m; i++)
	{
		TVGlobal_v[i]++;
		TVGlobal_v[i + m + 1]++;
	}

	trainedModel_v[0] = K2 * log(((long double)TVGlobal_v[0]) / (long double)TVGlobal_v[TVGlobal_v.size() - 1])
		- K2 * log((long double)TVGlobal_v[m + 1] / (long double)TVGlobal_v[TVGlobal_v.size() - 1]);
	//trainedModel_v[m + 1] = K2 * log((long double)TVGlobal_v[m + 1] / (long double)TVGlobal_v[TVGlobal_v.size() - 1]);
	//cout << "trainedModel_v[0]="<< trainedModel_v[0] << endl;

	for (int i = 1; i <= m; i++)
	{
		trainedModel_v[i] = K2 * log((long double)TVGlobal_v[i] / (long double)TVGlobal_v[0])
			- K2 * log((long double)TVGlobal_v[i + m + 1] / (long double)TVGlobal_v[m + 1]);
		//trainedModel_v[i + m + 1] = K2 * log((long double)TVGlobal_v[i + m + 1] / (long double)TVGlobal_v[m + 1]);
		//cout << "trainedModel_v["<<i<<"]=" << trainedModel_v[i] << endl;
	}

	//tmpmNEmail.getSelectedFeaturesParameters(selectedFeaturesOutput, logOfProbWordIsHamOutput,
	//	logOfProbWordIsSpamOutput, logProbHamOutput, logProbSpamOutput);
	//std::cout << "logOfProbWordIsHamOutput.size()=" << logOfProbWordIsHamOutput.size() << endl;
	//std::cout << "logOfProbWordIsSpamOutput.size()=" << logOfProbWordIsSpamOutput.size() << endl;

	//trainedModel_v[0] = (int64_t) (K2 * logProbHamOutput);
	//trainedModel_v[m + 1] = (int64_t)(K2 * logProbSpamOutput);
	//std::cout << "trainedModel_v[0]=" << trainedModel_v[0] << endl;
	//std::cout << "trainedModel_v["<< m + 1<<"]=" << trainedModel_v[m + 1] << endl;


	//for (int i = 1; i <= m; i++)
	//{
	//	trainedModel_v[i] = (int64_t)(K2 * logOfProbWordIsHamOutput[i]);
	//	std::cout << "trainedModel_v[" << i << "]=" << trainedModel_v[i] << endl;
	//	trainedModel_v[i + m + 1] = (int64_t)(K2 * logOfProbWordIsSpamOutput[i]);
	//	std::cout << "trainedModel_v[" << i + m + 1 << "]=" << trainedModel_v[i + m + 1] << endl;
	//}

	int queriesPerCphrtxt = polyModulus / ((m + 1));
	for (int i = 1; i < queriesPerCphrtxt; i++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			trainedModel_v[i * (m + 1) + j] = trainedModel_v[j];
			//trainedModel_v[i * 2 * (m + 1) + j + m + 1] = trainedModel_v[j + m + 1];
		}
	}

	vector<vector<int>> confusionMatrix(2, vector<int>(2, 0));
	std::cout << "DOING THE PLAIN CLASSIFICATION NOW" << endl;
	int countQuery = 0;

	string hamPath = path + "\\ham", word;

	vector<string> hamFiles, spamFiles;

	for (const auto & entry : fs::filesystem::recursive_directory_iterator(hamPath))
	{
		filename = entry.path().string();
		hamFiles.push_back(filename);
	}

	for (int j = 0; j < hamFiles.size(); j += queriesPerCphrtxt)
	{
		//std::cout << "Here 0" << endl;
		vector<int64_t> queryVector_v(polyModulus, 0);
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = hamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			//if (inputFile.fail())
			//{
			//	cout << "could not open file " << filename << endl;
			//	cin.get();
			//	cin.ignore();
			//	exit(0);
			//}

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * 2 * (m + 1) + m + 1] = 1;

			//std::cout << "Here 1" << endl;

			int begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
			while (inputFile >> word)
			{
				//wordCount++;
				//existsInList = false;
				begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
				while (!((begIndex + 1) == endIndex))
				{
					if (word == selectedFeaturesOutput[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex +1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > selectedFeaturesOutput[midIndex])
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
				if (word == selectedFeaturesOutput[begIndex])
				{
					queryVector_v[i *(m + 1) + begIndex+1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == selectedFeaturesOutput[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex+1]++;
					//queryVector_v[i * 2 * (m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			vector<int64_t> logQueryProb(queriesPerCphrtxt, 0);
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			for (int k = 0; k <= m; k++)
			{
				queryVector_v[i *(m + 1) + k] *= trainedModel_v[i * (m + 1) + k];
				logQueryProb[i] += queryVector_v[i * (m + 1) + k];
				//queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[i * 2 * (m + 1) + k + m + 1];
				//logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			}
			//std::cout << "Here 2" << endl;
			if (logQueryProb[i] > 0)
				confusionMatrix[0][0]++;
			else
				confusionMatrix[0][1]++;

			//countQuery++;
			//if ((countQuery % 100) == 0)
			//{
			//	std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
			//	std::system("pause");
			//}
		}
	}

	countQuery = 0;
	string spamPath = path + "\\spam";
	for (const auto & entry : fs::filesystem::recursive_directory_iterator(spamPath))
	{
		filename = entry.path().string();
		spamFiles.push_back(filename);
	}

	for (int j = 0; j < spamFiles.size(); j += queriesPerCphrtxt)
	{
		vector<int64_t> queryVector_v(polyModulus, 0);
		for (int i = 0; i < queriesPerCphrtxt; i++)
		{
			filename = spamFiles[j + i];
			//std::cout <<"opening file: "<< filename << std::endl;

			ifstream inputFile(filename);
			//if (inputFile.fail())
			//{
			//	cout << "could not open file " << filename << endl;
			//	cin.get();
			//	cin.ignore();
			//	exit(0);
			//}

			queryVector_v[i * (m + 1) + 0] = 1;
			//queryVector_v[i * 2 * (m + 1) + m + 1] = 1;

			int begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
			while (inputFile >> word)
			{
				//wordCount++;
				//existsInList = false;
				begIndex = 1, endIndex = m, midIndex = (begIndex + endIndex) / 2;
				while (!((begIndex + 1) == endIndex))
				{
					if (word == selectedFeaturesOutput[midIndex])
					{
						queryVector_v[i * (m + 1) + midIndex + 1]++;
						//queryVector_v[i * 2 * (m + 1) + midIndex + 1 + m + 1]++;
						break;
					}
					else if (word > selectedFeaturesOutput[midIndex])
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
				if (word == selectedFeaturesOutput[begIndex])
				{
					queryVector_v[i * (m + 1) + begIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + begIndex + 1 + m + 1]++;
				}
				else if (word == selectedFeaturesOutput[endIndex])
				{
					queryVector_v[i * (m + 1) + endIndex + 1]++;
					//queryVector_v[i * 2 * (m + 1) + endIndex + 1 + m + 1]++;
				}
			}

			vector<int64_t> logQueryProb(queriesPerCphrtxt, 0);
			//	logSpamQueryProb(queriesPerCphrtxt, 0);

			for (int k = 0; k <= m; k++)
			{
				queryVector_v[i * (m + 1) + k] *= trainedModel_v[i * (m + 1) + k];
				logQueryProb[i] += queryVector_v[i * (m + 1) + k];
				//queryVector_v[i * 2 * (m + 1) + k + m + 1] *= trainedModel_v[i * 2 * (m + 1) + k + m + 1];
				//logSpamQueryProb[i] += queryVector_v[i * 2 * (m + 1) + k + m + 1];
			}

			if (logQueryProb[i] < 0)
				confusionMatrix[1][1]++;
			else
				confusionMatrix[1][0]++;

			//countQuery++;
			//if ((countQuery % 100) == 0)
			//{
			//	std::cout << "Processing the " << countQuery << "th ham QUERY record" << endl;
			//	std::system("pause");
			//}
		}
	}

	std::cout << endl << "Printing the confusion matrix after PP classification in the end" << endl << endl;
	for (int i = 0; i < confusionMatrix.size(); i++)
	{
		std::cout << endl;
		for (int j = 0; j < confusionMatrix[0].size(); j++)
		{
			std::cout << confusionMatrix[i][j] << "\t";
		}
	}
	std::cout << endl << "The accuracy of the Plain classification is: " << (long double)100 * (confusionMatrix[0][0] + confusionMatrix[1][1]) / (long double)((confusionMatrix[0][0] + confusionMatrix[1][1] + confusionMatrix[0][1] + confusionMatrix[1][0])) << " %" << endl << endl;


	std::cout << endl;
	std::cout << "**************************************************************************************" << endl;
	std::cout << "The MULTINOMIAL NAIVE BAYES  secMLClass (secure ML Classification)-CLIENT CENTRIC CASE" << endl;
	std::cout << "**************************************************************************************" << endl;

	PPServCenMNClass(path, selectedFeaturesOutput,
		trainedModel_v, polyModulus, m, output);
	//std::system("pause");

	std::cout << endl;
	std::cout << "**************************************************************************************" << endl;
	std::cout << "The MULTINOMIAL NAIVE BAYES  secMLClass (secure ML Classification)-SERVER CENTRIC CASE" << endl;
	std::cout << "**************************************************************************************" << endl;
	 	  
	PPClientCenMNClass(path, selectedFeaturesOutput,
		trainedModel_v, polyModulus, m, output);
	//std::system("pause");

	std::cout << endl << endl << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "C) DEALING WITH THE ADFA-NB15 CYBERSECURITY DATASET NOW" << endl;
	std::cout << "********************************************************" << endl;
	std::cout << "********************************************************" << endl << endl;

	std::cout << endl;
	std::cout << "*********************************************" << endl;
	std::cout << "The SVM and LR ML Classification - PLAIN CASE" << endl;
	std::cout << "*********************************************" << endl;
	
	vector<double> tmp_weightsADFA(16);
	vector<double> tmp_weightsADFA2(16);
	vector<int64_t> svm_weightsADFA(polyModulus);
	vector<int64_t> logistic_weightsADFA(polyModulus);
	vector<vector<int64_t>> dataADFA(ceil((double)(82332 * 16) / polyModulus), vector<int64_t>(polyModulus, 0));

	vector<double> res_svmADFA(82332);// 82332 is the number of instances (rows) in the ADFA-NB15 CYBERSECURITY DATASET
	vector<double> res_logisticADFA(82332);

	vector<int> true_labelADFA(82332);
	vector<int> predicted_label_svmADFA(82332);
	vector<int> predicted_label_logisticADFA(82332);
	// SVM weights, multiplied with 1000 and rounded
	tmp_weightsADFA[0] = -1140; //-1.1404;
	tmp_weightsADFA[1] = 5; //0.0048;
	tmp_weightsADFA[2] = -25;  //-0.0246;
	tmp_weightsADFA[3] = 24; //0.0240;
	tmp_weightsADFA[4] = 0;  //0;
	tmp_weightsADFA[5] = 0;  //0;
	tmp_weightsADFA[6] = 1070;//1.0702;
	tmp_weightsADFA[7] = 0; //0;
	tmp_weightsADFA[8] = 0;  //0;
	tmp_weightsADFA[9] = 0; //0;
	tmp_weightsADFA[10] = 0;//0;
	tmp_weightsADFA[11] = 0;//0;
	tmp_weightsADFA[12] = 0;//0;
	tmp_weightsADFA[13] = 0;//0;
	tmp_weightsADFA[14] = 0;//0;
	tmp_weightsADFA[15] = -2000; //-2
	// Logistic Regression weights, multiplied with 10 and rounded
	tmp_weightsADFA2[0] = -233 ;//9.6728;
	tmp_weightsADFA2[1] = 2;//-0.5313;
	tmp_weightsADFA2[2] = -7;//-0.0069;
	tmp_weightsADFA2[3] = 5;//-0.3301;
	tmp_weightsADFA2[4] = 122;//-0.2393;
	tmp_weightsADFA2[5] = -8;//-0.0676;
	tmp_weightsADFA2[6] = 248;//-0.4068;
	tmp_weightsADFA2[7] = -2;//-0.4093;
	tmp_weightsADFA2[8] = 1;//-0.1463;
	tmp_weightsADFA2[9] = 37;//-0.5488;
	tmp_weightsADFA2[10] = 21;
	tmp_weightsADFA2[11] = 59;
	tmp_weightsADFA2[12] = 59;
	tmp_weightsADFA2[13] = -11;
	tmp_weightsADFA2[14] = -23;
	tmp_weightsADFA2[15] = -883;
	// Weights of SVM and Logistic Regression stacked to vector with size 16384
	count1 = 0;
	count2 = 0;
	while (count1 < polyModulus)
	{
		svm_weightsADFA[count1] = 1;
		count2 = 0;
		while (count2 < 16)
		{
			svm_weightsADFA[count1] = tmp_weightsADFA[count2];
			logistic_weightsADFA[count1] = tmp_weightsADFA2[count2];
			count1 += 1;
			count2 += 1;
		}
	}
	
	// Data Reading
	ifstream fileADFA;
	filename = "UNSW-testing.txt"; // Name of the data txt file
	line = "";
	fileADFA.open(filename.c_str());
	count1 = 0;
	count2 = 0;
	count_tt = 0;
	data_counter = -1;
	rate = 0;
	trash;
	maxValue = 0;
	int vec_counter = 0;
	while (getline(fileADFA, line)) // Taking line from data txt file
	{
		istringstream strfile2(line);

		if (count1 == polyModulus)
		{
			count1 = 0;
			vec_counter++;
		}
		dataADFA[vec_counter][count1] = 1;
		count2 = 0;
		count1 += 1;
		data_counter += 1;
		while (strfile2 >> rate)
		{
			if (count2 != 15)
			{
				dataADFA[vec_counter][count1] = rate;
				count2 += 1;
				count1 += 1;
				//cout << rate << " ";
			}
			else
			{
				true_labelADFA[data_counter] = rate;
				count2 += 1;
				//cout << "\n";
			}

			if (count2 != 16)
			{
				strfile2 >> trash;
			}
		}		

	}
	
	// Results Multiplication
	
	tmp_svm = 0;
	tmp_logistic = 0;
	double svm_acc = 0;
	double logistic_acc = 0;
	maxValue = 0;
	index = 0;
	for (int i = 0; i < dataADFA.size(); i++)
	{
		for (int j = 0; j < dataADFA[0].size(); j+=16)
		{
			count2 = 0;
			tmp_svm = 0;
			tmp_logistic = 0;

			while (count2 < 16)
			{
				tmp_svm += svm_weightsADFA[j+count2] * dataADFA[i][j + count2];
				tmp_logistic += logistic_weightsADFA[j + count2] * dataADFA[i][j + count2];
				count2 += 1;			
			}
			res_svmADFA[index] = tmp_svm;
			res_logisticADFA[index] = tmp_logistic;
			
			if (maxValue < res_svmADFA[index])
			{
				maxValue = res_svmADFA[index];
			}
			if (maxValue < res_logisticADFA[index])
			{
				maxValue = res_logisticADFA[index];
			}
			//cout << tmp_svm << "\n";
			// SVM final prediction
			if (tmp_svm >= 0)
			{
				predicted_label_svmADFA[index] = 1;
			}
			else
			{
				predicted_label_svmADFA[index] = 0;
			}
			// Logistic regression final prediction
			if (tmp_logistic >= 0)
			{
				predicted_label_logisticADFA[index] = 1;
			}
			else
			{
				predicted_label_logisticADFA[index] = 0;
			}
			if (predicted_label_svmADFA[index] == true_labelADFA[index])
			{
				svm_acc++;
			}
			if (predicted_label_logisticADFA[index] == true_labelADFA[index])
			{
				logistic_acc++;
			}
			index++;
			if (index == 82332)
				break;
		}		
	}
		
	std::cout << "SVM Plain Accuracy: " << svm_acc * 100 / 82332 << "%" << "\n";
	std::cout << "Logistic Regression Plain Accuracy: " << logistic_acc * 100 / 82332 << "%" << "\n";
	//std::cout << "maxValue=" << maxValue << endl;

	std::cout << endl;
	std::cout << "************************************************************************" << endl;
	std::cout << "The SVM and LR secMLClass (secure ML Classification)-CLIENT CENTRIC CASE" << endl;
	std::cout << "************************************************************************" << endl;
	
	std::cout << "Calling the client centric function for both SVM and LR" << endl;
	//system("pause");
	SVM_LR_ClientCentric(polyModulus, svm_weightsADFA, logistic_weightsADFA, dataADFA, gal_keys,
		decryptor, true_labelADFA, batch_encoder, encryptor, evaluator, 82332, 1);
	
	std::cout << endl;
	std::cout << "************************************************************************" << endl;
	std::cout << "The SVM and LR secMLClass (secure ML Classification)-SERVER CENTRIC CASE" << endl;
	std::cout << "************************************************************************" << endl;
	
	std::cout << "Calling the client centric function for both SVM and LR" << endl;
	//system("pause");
	SVM_LR_ServerCentric(polyModulus, svm_weightsADFA, logistic_weightsADFA, dataADFA, gal_keys,
		decryptor, true_labelADFA, batch_encoder, encryptor, evaluator, 82332, 1);

	
	std::cout << "Everything is working fine!" << endl;
	system("pause");
    return 0;
}