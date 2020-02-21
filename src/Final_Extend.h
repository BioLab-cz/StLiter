#ifndef _Final_Extend_
#define _Final_Extend_

#include "Basic.h"
#include "BplusTreeBit.h"
#include "Hash.h"
#include "litool.h"
#include "Extend_Kmer.h"

struct para_Count
{
	uint32_t T_i;
	char *seq;
	uint64_t seq_length;
	uint64_t Tasksize;
	uint32_t original_kmer_length;
	uint64_t start_line;
	uint64_t *end_line;
	char* FileName;
	uint32_t thread_num;
	uint32_t pos_label;
	uint32_t InOut_label;
};
struct para_Final
{
	uint32_t T_i;
	char *seq;
	uint64_t seq_length;
	uint64_t Tasksize;
	uint32_t original_kmer_length;
	uint32_t extend_length;
	uint64_t start_line;
	uint64_t *end_line;
	uint32_t InOut_label;
	char*inputFileName;
	uint32_t ad_label;
	uint32_t thread_num;
	uint64_t Task_Size_Single_Thread;
	uint32_t pos_label;
	uint64_t **p_candidate_pos;
	uint64_t *template_64;
	uint64_t *p_candidate_pos_len;
	uint32_t rc_label;
};

/***************************************Final for 256bitKmer******************/
struct TwoDimenHashV_count_256bitKmer
{
	uint64_t* hv;
	uint64_t pos;
};
struct Parallel_para_count01_256bitKmer
{
	char *seq;
	uint64_t seq_length;
	uint32_t kmer_len1;
	uint64_t task_start_pos;
	uint64_t task_len;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV_count_256bitKmer* p_hv_tmp;
	struct bit256KmerPara para;
	uint64_t *p_candidate_current;
	uint64_t *template_64;
	uint32_t pos_label;
};
struct Parallel_para_count02_256bitKmer
{
	char *seq;
	struct NodeBit ** p_root;
	uint32_t kmer_len1;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV_count_256bitKmer** p_hv_tmp;
	uint32_t thread_id;
	uint32_t thread_num;
	struct bit256KmerPara para;
};
struct TwoDimenHashV_Final_256bitKmer
{
	uint64_t pos;
	struct nodeBit tmp;
};
struct Parallel_para_Final01_hash_256bitKmer
{
	char *seq;
	uint64_t seq_length;
	struct bit256Hash * p_root;
	uint32_t kmer_len1;
	uint32_t kmer_len2;
	uint64_t task_start_pos;
	uint64_t task_len;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV_Final_256bitKmer* p_hv_tmp;
	struct bit256KmerPara para1;
	struct bit256KmerPara para2;
	uint64_t *p_candidate_pos;
	uint64_t *p_candidate_current;
	uint64_t *template_64;
	uint32_t pos_label;
};
struct Parallel_para_Final02_hash_256bitKmer
{
	char *seq;
	struct NodeBit ** p_root_hashtable;
	uint32_t kmer_len1;
	uint32_t kmer_len2;
	uint32_t InOut_label;
	uint64_t *result_len;
	struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp;
	uint32_t thread_id;
	uint32_t thread_num;
	struct bit256KmerPara para;
};
struct reslt_para_plustree_final_256bitKmer
{
	struct NodeBit * p_root_start;
	struct NodeBit * p_root_end;
	uint64_t ** reslt;
	uint32_t ** reslt_c;
	uint64_t * reslt_size;
	uint32_t InOut_label;
	uint32_t ad_label;
	struct bit256KmerPara para;
};
void Insert_ad_256bitKmer(struct NodeBit** p_root,struct nodeBit insert_data,struct bit256KmerPara para);
void Counting_256bitKmer(struct para_Count tmp);
void Cal_Single_Task_4bit_UsingHashTable_withLabel_Final_256bitKmer(struct para_Final tmp);
/***************************************Final for 256bitKmer******************/
#endif
