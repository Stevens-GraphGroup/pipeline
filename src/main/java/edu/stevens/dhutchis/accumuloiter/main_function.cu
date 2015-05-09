/* ### this is a .cu file### */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

#include "header_def.cuh"
#include "data_partition.cuh"		/* only for multi-GPU */

extern "C" bool* mymain(int sizeout, const char** s, const char *Chmm_path )
{
	const char *MODEL_PATH = Chmm_path;

	/* Good upto size == 800 (temporary fixed) */
	int ACTIVE_BLOCKS = 32;
	int ROW_IN_BLOCK = 16;

	/* Once we have multiple devices, display them and set them as 'current' for ready to run */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);

	int device;
	for (device = 0; device < deviceCount; ++device)
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, device);
		//printf("\nDevice %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
	}

	HMMER_PROFILE *hmm = NULL;						/* claim hmm         */
	hmm = hmmer_profile_Create(MODEL_PATH);			/* alloc mem for hmm */

	/* ************************************************ */
	/* 		1. Parameters, match, insert emission 		*/
	/* ************************************************ */

	/* get parameters (M, MU[], LAMBDA[]) */
	get_Parameters(hmm, MODEL_PATH);

	/* Get raw match & insert emission */
	get_Emission(hmm, MODEL_PATH);

	/* Transfer to prob from raw data (need do it before get degen)*/
	log_Odd_score(hmm);

	/* Get degenerate residues */
	mat_degen(hmm);

	/* Rescaling and rounding */
	mf_conversion(hmm);

	//Until now, all preparation works of MSV has been done!


	/* ************************************************************************ */
	/* 		2. Transition probability only for Viterbi and Forward filter 		*/
	/* ************************************************************************ */

#if 0

	/**/ get_transition(hmm, MODEL_PATH) ;
	/**/ get_entryScore(hmm);
	/**/ log_Trans(hmm) ;
	/**/ xTrans(hmm) ;
	/**/ vf_conversion(hmm) ;

	//Until now, all preparation works of Viterbi has been done!

#endif	/* Configure for Viterbi */

	/* **************************************************************************** */
	/* 		3. Protein sequence database read and pre-process               		*/
	/* **************************************************************************** */

	int number = 0;		/* the total number of seq in this database */

	/**/
	number = sizeout; //get_Seqnumber(DATABASE_PATH);
	//printf("\n In this database, we have %d seqs.\n", number);

	/**/
	const char** seq = (const char**)malloc(number * sizeof(char*));     					//dynamic memory for address value of each sequence: seq[i]->an address value
	int* seq_len = (int*)malloc(number * sizeof(int));       								//for cache length of each sequence

	//if(alloc_Eachseq(seq, seq_len, number, DATABASE_PATH) != 1) printf("error!\n");
	//if(fill_Eachseq(seq, seq_len, number, DATABASE_PATH) != 1) printf("error!\n");

	seq = s;	/* copy pointer from outside */

	for(int i =0;i< number;i++)		//sets seq_length
	{
		seq_len[i] = strlen(seq[i]);
	}

	/* transfer seq into int values */
	char** iSeq = (char**)malloc(number * sizeof(char*));
	for(int i = 0; i < number; i++)
		iSeq[i] = (char*)malloc(seq_len[i] * sizeof(char));
	if(seq_To_char(iSeq, seq, seq_len, number) != 1) printf("error in transfer to short from char !\n");


	/* ***************************************** */
	/* 		6. Parameters for kernel launch	 	 */
	/* ***************************************** */

	/* gird & block */
	dim3 gridDim = dim3(1, ACTIVE_BLOCKS, 1);		//This is better since it is able to allocated paras dynamically.
	dim3 blockDim = dim3(32, ROW_IN_BLOCK, 1);

	//revolution begin
	int SPAN = number / deviceCount;
	int REMAIN = number % deviceCount;

	GPUs *kernel = new GPUs[deviceCount];	//dynamic allocation..

	//can try to use multiple threads...
	for(int i = 0; i < deviceCount; i++)
	{
		cudaSetDevice(i);									/* CURRENT GPU */
		cudaDeviceReset();								  	/* reset current GPU */
		kernel[i].get_duty(SPAN, REMAIN, i, deviceCount);	/* distribution */
		kernel[i].compress(iSeq, seq_len);					/* residues pack */
		kernel[i].get_Real_Length(seq_len);					/* real Length of everyone */
		kernel[i].trans_to_1D();							/* 2D array to 1D array */
		kernel[i].host_to_device(hmm);
	}

	/* launch */
	for (int i = 0; i < deviceCount; i++)				 //Now we havent consider the time cost of data copy yet
	{
		cudaSetDevice(i);
		kernel[i].launch_MSV(gridDim, blockDim, hmm, i, ROW_IN_BLOCK);
	}
	cudaDeviceSynchronize();

	/* get back results and counter */
	bool* passed_msv;
	passed_msv = (bool*)malloc(number*sizeof(bool));

	int counter = 0;
	for (int i = 0; i < deviceCount; i++)
	{
		cudaSetDevice(i);
		kernel[i].device_to_host(hmm,passed_msv);
	}

	for(int i = 0; i < number ; i++)
	{
		if(passed_msv[i] == 1)
		{
			counter++;
		}
	}

	/* ***************************************** */
	/* 		  6. Host & Device mem free   	 	 */
	/* ***************************************** */

	/**/
	delete[] kernel;

	/* Free */
	free(seq_len);

	for(int i = 0; i < number; i++) {
		free(iSeq[i]);
	}
	free(iSeq);		/* PLEASE REMIND: free child first, then parent */

	/* HMM model */
	freeHMM(hmm);

	return passed_msv;
}

