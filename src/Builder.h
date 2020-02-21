/*
 * Builder.h
 *
 *  Created on: May 15, 2018
 *      Author: bio
 */

#ifndef BUILDER_H_
#define BUILDER_H_
#include "Basic.h"

struct para_unipath_bit
{
	uint32_t kmer_len;
	char * p_tmp;
	struct NodeBit * p_Branch_root;
	uint64_t start;
	uint64_t end;
	uint64_t **reslt;
	uint64_t **reslt_c;
	uint64_t *reslt_size;
	uint64_t **reslt_bran_pos;
	uint64_t *reslt_size_bran_pos;
	uint32_t thread_id;
	uint32_t thread_num;
	struct bit256KmerPara para;

	uint32_t pos_label;
	uint64_t *p_candidate_pos_in;
	uint64_t *template_64;
	uint64_t *p_candidate_pos_out;
};
struct para_findUnipath
{
	char *seq;
	uint32_t original_kmer_length;
	uint32_t thread_num;
	uint32_t pos_label;
	uint32_t rc_label;
	uint32_t memory_size;
	uint32_t pl_label;

};
void FindUniPath(struct para_findUnipath tmp);
void FindUniPath_new(struct para_findUnipath tmps);
void FindUniPath_novel(struct para_findUnipath tmps);
void FindUniPath_fast(struct para_findUnipath tmps);
void FindUniPath_Novel(struct para_findUnipath tmps);

#endif /* BUILDER_H_ */
