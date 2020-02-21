/*
 * BplusTreeBit.h
 *
 *  Created on: Nov 12, 2018
 *      Author: bio
 */

#ifndef BPLUSTREEBIT_H_
#define BPLUSTREEBIT_H_

#include "Basic.h"

struct nodeBit{
	uint64_t* 	hashValue;
	uint64_t 	arrayID;
	uint8_t		ad;
	uint64_t*	p_posList;			//let the reference length is less than 2^32
	uint64_t	p_posList_length;
};
struct NodeBit{
	struct nodeBit 		data[M];
	struct NodeBit *	p_child[M];
	struct NodeBit *	parent;
	struct NodeBit *	brother;
	uint32_t 			Node_Size;
	uint32_t 			leaf_label;
	int32_t 			numAsChild;
};

void Node_initial_bit(struct NodeBit * p);
void Insert_Value_bit(struct NodeBit** p_root,\
		struct nodeBit insert_data,struct bit256KmerPara para);
void Insert_Node_bit(struct NodeBit** p_root,struct NodeBit* p_inserted_node,\
		struct nodeBit insert_data,struct NodeBit* p_middle,\
		struct bit256KmerPara para);
void Divide_Node_bit(struct NodeBit** p_root,struct NodeBit* p_divide_node,\
		struct nodeBit insert_data,\
		uint32_t pos,struct NodeBit* insert_child,uint32_t ad_label,\
		uint32_t poslist_label,struct bit256KmerPara para);
void update_root_bit(struct NodeBit** p_root,struct NodeBit* p2,\
		struct bit256KmerPara para);
void createBplusTree_bit(struct NodeBit** p_root,struct nodeBit* p_data,\
		uint32_t p_data_size,struct bit256KmerPara para);
void destory_tree_bit(struct NodeBit* p_root,struct bit256KmerPara para);
struct NodeBit* Find_Minimal_Node_bit(struct NodeBit* p_root);
int64_t MappingHashValueToID_bit(struct NodeBit* p_root,\
		uint64_t* hashValue,struct bit256KmerPara para);
void count_hashValue_bit(struct NodeBit* p_root,\
		uint64_t* hashValue,struct bit256KmerPara para);
int64_t addAD_bit(struct NodeBit* p_root,\
		uint64_t* hashValue,uint8_t c, struct bit256KmerPara para);
int64_t addAD1_bit(struct NodeBit** p_root,\
		uint64_t* hashValue,uint8_t c, struct bit256KmerPara para);
int64_t addAD2_bit(struct NodeBit* p_root,\
		uint64_t* hashValue,uint8_t c,struct bit256KmerPara para);
int64_t MappingHashValueToID_bit_unipath_count(struct NodeBit* p_root,uint64_t* hashValue,struct bit256KmerPara para);//2019.10.23

#endif /* BPLUSTREEBIT_H_ */
