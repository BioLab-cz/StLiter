/*
 * Basic.h
 *
 *  Created on: Apr 14, 2018
 *      Author: ycy
 */

#ifndef BASIC_H_
#define BASIC_H_

#include<cstdio>
#include <map>
#include<cstdlib>
#include<cstring>
#include <ctime>
#include<cmath>
#include <fstream>
#include <iostream>
#include "stdint.h"
#include "stdlib.h"
#include <pthread.h>
#include <sys/time.h>
using namespace std;

//#define _Multi_thread_
#define M 31
#define _BufferSize_ 257
struct RefFilePath
{
	char **pRefFilePath;
	uint32_t NumberOfPathes;
};
struct para_getN
{
	uint64_t kmerlen;
	char* FileName;
	uint32_t InOut_label;
	uint32_t isMiddle;
	uint32_t isposition;
};
struct CDBG_para
{
	char *p_ref;
	char method;
	uint32_t memory;
	uint32_t kmer_len;
	uint32_t extend_len;
	uint32_t inout_label;
	uint32_t isFinal_label;
	uint32_t ad_label;
	uint32_t thread_num;
	uint64_t Task_Size_Single_Thread;
	uint32_t pos_label;
	uint32_t sort_label;
	uint32_t rc_label;
	uint32_t pl_label;
};
struct bitKmerPara{
	uint32_t kmer1Len;
	uint32_t kmer8Len;
	uint32_t kmer64Len;
	uint32_t remainer1to8;
	uint32_t remainer1to64;
	uint32_t remainer8to64;
};
struct bit256KmerPara{
	uint32_t kmer1Len;
	uint32_t kmer64Len;
	uint32_t remainer1to64;
	uint64_t codefor1;
};

void getFileName(struct para_getN tmp);
uint64_t cal_hash_value_indirectly(char *seq,uint64_t current,uint32_t len);
uint64_t cal_hash_value_directly(char *seq,uint32_t len);
void Transition(char start,char end,uint8_t* temp_res);
void Transition_4bit(char c_ad,uint8_t* temp_res,uint32_t label_4bit);
uint32_t Judge(uint8_t temp_a);
uint32_t Judge_new(uint8_t temp_a);
uint32_t Judge_novel(uint8_t temp_a);
uint32_t Judge_4bit(uint8_t temp_a,uint32_t label_4bit);
void ReadSeq(char **seq1,uint64_t *seq_length,char* p_ref);
void rc(char **seq_or, uint64_t seq_length);
void ReadSeq_identical(char* p_ref);
void getRefFilePathes(char* pathFile, struct RefFilePath* p);

/***************************ForBitCharKmer*********************************/
void fill_char_with_four_char(uint8_t* current,char*p);
void ReadSeq_bit(uint8_t **seq1,uint64_t *seq_length,char* p_ref);
void Fetch_kmer_on_bit_sequence(uint8_t **seq,uint64_t start, struct bitKmerPara para,uint8_t *kmer);
void calMulDimHashValue(uint8_t *kmer,uint64_t *hv64,struct bitKmerPara para);
void Transition_bitchar(uint8_t start,uint8_t end,uint8_t* temp_res);
/***************************ForBitCharKmer*********************************/

/***************************For256BitCharKmer******************************/
void cal_hash_value_directly_256bit(char *seq,uint64_t * kmer,\
		struct bit256KmerPara para);
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para);
void cal_newExtendedKmerOut(uint64_t* original,uint64_t* newKmer,uint64_t j,uint32_t exLen,struct bit256KmerPara para);
void cal_newExtendedKmerIn(uint64_t* original,uint64_t* newKmer,uint64_t j,uint32_t exLen,struct bit256KmerPara para);
/***************************For256BitCharKmer******************************/

uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len);



#endif /* BASIC_H_ */
