/*
 * Task_allocate.h
 *
 *  Created on: Apr 14, 2018
 *      Author: ycy
 */

#ifndef EXTEND_KMER_H_
#define EXTEND_KMER_H_

#include "Basic.h"
#include "Hash.h"
#include "Final_Extend.h"
#include "litool.h"
struct TwoDimenHashV
{
	uint64_t pos;
	uint64_t hv1;
	uint64_t hv2;
};
struct Parallel_para_extend04_mode
{
	char *seq;
	uint8_t **HashAd;
	uint32_t kmer_len1;
	uint32_t kmer_len2;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV** p_hv_tmp;
	uint32_t thread_id;
	uint32_t thread_num;
};
struct para_calTask
{
	char *seq;
	uint32_t pos_label;
	uint32_t extend_length;
	uint32_t original_kmer_length;
	uint32_t InOut_label;
	uint32_t isFinal;
	uint32_t ad_label;
	uint64_t memory_size;
	uint32_t thread_num;
	uint64_t Task_Size_Single_Thread;
	uint32_t sort_label;
	uint32_t rc_label;
};
struct create_HashTable_para
{
	struct bit256Hash* p_root;
	uint64_t *h1;
	uint32_t thread_num;
	uint32_t thread_id;
	uint64_t Tasksize_current;
	struct bit256KmerPara para1;
};
struct para_bPlusTreeMap
{
	uint32_t T_i;
	char *seq;
	uint64_t seq_length;
	uint64_t start_line;
	uint64_t *end_line;
	uint64_t Tasksize;
	uint32_t extend_length;
	uint32_t original_kmer_length;
	uint32_t InOut_label;
	uint32_t ad_label;
	uint32_t pos_label;
	uint32_t thread_num;
	char*inputFileName;
	uint64_t Task_Size_Single_Thread;
	uint64_t **p_candidate_pos;
	uint64_t *p_candidate_pos_len;
	uint64_t *template_64;
	uint32_t rc_label;
};
struct Parallel_para_extend03_hashtable
{
	char *seq;
	uint64_t seq_length;
	struct bit256Hash* p_root;
	uint32_t kmer_len1;
	uint32_t kmer_len2;
	uint64_t task_start_pos;
	uint64_t task_len;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV* p_hv_tmp;
	struct bit256KmerPara para;
	uint64_t *p_candidate_pos;
	uint64_t *template_64;
	uint32_t pos_label;
	uint64_t *p_candidate_current;
};
struct reslt_para_plustree_hashtable
{
	uint64_t *h1;
	uint64_t start,end;
	uint64_t ** reslt;
	uint8_t ** reslt_c;
	uint64_t * reslt_size;
	uint64_t each_hashTable_size;
	uint32_t extend_length;
	uint32_t ad_label;
	uint8_t **HashAd;
	struct bit256KmerPara para1;
	struct bit256KmerPara para2;
};
void Cal_Task(struct para_calTask tmp);
void *create_HashTabel(void *arg);
void Cal_Single_Task_4bit_UsingHashTable_withLabel_256bitKmer(struct para_bPlusTreeMap tmp);
#endif /* EXTEND_KMER_H_ */
