/*
 * litool.h
 *
 *  Created on: Nov 19, 2018
 *      Author: bio
 */

#ifndef LITOOL_H_
#define LITOOL_H_

#include "Basic.h"
#include "BplusTreeBit.h"
#include "Final_Extend.h"
void Check_bitBybit(char *p1,char *p2);
void transKmerBinaryToText(char *p1,char*p2,uint32_t n);
void transCountKmerBinaryToText(char *p1,char*p2,uint32_t n);
void CheckResult_num(char *p1);
void CheckResult_content(char *p1,char *p2);
void CheckResult_bitfileIDentical(char *p1,char *p2);
void SortOutKmersUsingBplusTreebinarytobinary(char*inputFileName,char*outputFileName,uint32_t len);
void SortOutCountKmersUsingBplusTreebinarytobinary(char*inputFileName,char*outputFileName,uint32_t len);
void binaryCoutInt(uint64_t a);
void Coutsubstring(char* p,uint32_t l);
void binaryCoutChar(uint8_t a);
void SortOutCountKmersSingle(char*inputFileName,char*outputFileName,uint32_t len);
void SortKmers(char*inputFileName,char*outputFileName,char*inputFileName_ad,char*outputFileName_ad,uint32_t len);
#endif /* LITOOL_H_ */
