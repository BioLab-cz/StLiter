#ifndef _BEING_SHORT_
#define _BEING_SHORT_

#include "Basic.h"
struct reslt_para_being
{
	uint64_t start;
	uint64_t task_len;
	uint32_t thread_id;
	uint32_t thread_num;
	uint64_t ** reslt_in;
	uint64_t ** reslt_out;
	uint8_t ** reslt_in_ad;//20200207
	uint8_t ** reslt_out_ad;//20200207
	uint64_t ** reslt_uni;
	uint8_t ** reslt_uni_ad;//20200207
	uint64_t * reslt_size_uni;
	uint64_t * reslt_size_in;
	uint64_t * reslt_size_out;
	uint8_t **p_hash_table;
	uint32_t ad_label;//20200207
};
struct Parallel_para_being
{
	char* p_seq;
	uint8_t* p_hashtable;
	uint64_t seq_length;
	uint64_t* p_hashValueTable;
	uint32_t kmer_len;
	uint32_t thread_id;
	uint32_t thread_num;
	uint32_t start;
	uint32_t len;
};

struct para_being
{
	char* seq;
	uint64_t kmer_len_being;
	char *outputFileName;
	uint32_t thread_num;
	uint64_t Task_Size;
	uint32_t ad_label;//20200207
	uint32_t rc_label;//20200207
};

void Find_Short_Branched_Kmer(struct para_being tmp);
void Find_Short_Branched_Kmer_bit(struct para_being tmp);
void Find_Short_Branched_Kmer_256bit(struct para_being tmp);
void Find_Short_Branched_Kmer_256bit_novel(struct para_being tmp);

#endif
