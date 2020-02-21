/*
 * Hash.h
 *
 *  Created on: Oct 9, 2018
 *      Author: bio
 */

#ifndef HASH_H_
#define HASH_H_

#include "Basic.h"
#include "BplusTreeBit.h"
#include "Final_Extend.h"
#define HashSize 0x100000
//#define HashSize 0x10000
#define AddSize 1
#define MinSizeForSearch 100
#define HashFSize 0x1000000

/***************************hashForBitCharKmer*********************************/
struct bit256KmerPosition{
	uint64_t 	kmer[4];
	uint64_t 	position;
};
struct bit256Hash{
	uint64_t** p_kmer;
	uint64_t** p_pos;
	uint64_t* len;
	uint64_t* maxlen;
};
uint32_t bit256HashFunction(uint64_t* a,struct bit256KmerPara para);
uint32_t bit256HashFunction_left(uint64_t* a,struct bit256KmerPara para);
struct bit256Hash* initial256BitHashTable();
void insert256BitHashTable(struct bit256Hash *ph,uint64_t* kmer,uint64_t pos,bit256KmerPara para);
uint64_t binary256BitSearch(uint64_t * ppo,uint64_t* p,uint64_t len,uint64_t* y,bit256KmerPara para);
uint64_t search256BitHashTable(struct bit256Hash *ph,uint64_t* y,bit256KmerPara para);
void freeHash256BitTable(struct bit256Hash *ph);
/***************************hashForBitCharKmer*********************************/

uint32_t bit256hashFFunction(uint64_t *a,struct bit256KmerPara para);
struct NodeBit** bit256initialHashFTable();
struct NodeBit * bit256initialSingleHashFTable();
void bit256insertHashFTable(struct NodeBit** pk,struct nodeBit c_tmp_hashtable,struct bit256KmerPara para);
struct NodeBit * bit256Find_Minimal_Node_HashFTable(struct NodeBit** pk);
void bit256linkHashFTable(struct NodeBit** pk);
void bit256freeHashFTable(struct NodeBit **ph,bit256KmerPara para);
void insert256BitHashTable_parallel(struct bit256Hash *ph,uint64_t* kmer,uint64_t pos,bit256KmerPara para,uint32_t thread_num,uint32_t thread_id);
#endif /* HASH_H_ */
