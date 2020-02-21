/*
 * Task_allocate.cpp
 *
 *  Created on: Apr 14, 2018
 *      Author: bio
 */
#include "Extend_Kmer.h"

void *Allocate_HashTable(void *arg)
{
	struct Parallel_para_extend03_hashtable *tmp=(struct Parallel_para_extend03_hashtable *)arg;
	uint64_t task_start_pos=tmp->task_start_pos;

	char *seq=tmp->seq;
	uint64_t seq_length=tmp->seq_length;
	struct bit256Hash * p_root=tmp->p_root;
	uint32_t kmer_len1=tmp->kmer_len1;
	uint32_t kmer_len2=tmp->kmer_len2;
	uint64_t task_len=tmp->task_len;
	uint32_t InOut_label=tmp->InOut_label;
	struct bit256KmerPara para=tmp->para;

	struct TwoDimenHashV* p_hv_tmp=tmp->p_hv_tmp;
	uint64_t p_hv_tmp_length_cur=0;

	uint32_t pos_label=tmp->pos_label;
	uint64_t *p_candidate_pos=tmp->p_candidate_pos;
	uint64_t *template_64=tmp->template_64;
	uint64_t * p_candidate_current=tmp->p_candidate_current;

	uint64_t seq_k_current[4];
	uint32_t original_kmer_length=kmer_len1;
	char*temp_seq=seq+task_start_pos;
	uint32_t extend_length=kmer_len2;
	uint64_t HashAd_first;
	uint64_t HashAd_second;


	if(InOut_label==0)
	{
		cal_hash_value_directly_256bit(temp_seq,seq_k_current,para);
		if(task_start_pos!=0)
		{
			int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
			if(arrayID!=0)
			{
				if(pos_label==1||pos_label==3)
				{
					uint64_t tmp_pos_divide=task_start_pos/64;
					uint64_t tmp_pos_mod=task_start_pos%64;
					p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];
				}

				HashAd_first=arrayID-1;
				HashAd_second=cal_hash_value_directly(temp_seq+original_kmer_length,extend_length);
				p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos;
				p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
				p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
				p_hv_tmp_length_cur++;
			}
		}
		if(pos_label==3)
		{
			for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-extend_length)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				uint64_t tmp_pos_divide_current=(task_start_pos+j)/64;
				uint64_t tmp_pos_mod_current=(task_start_pos+j)%64;
				if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
				{
					int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
					if(arrayID!=0)
					{
						p_candidate_pos[tmp_pos_divide_current]=p_candidate_pos[tmp_pos_divide_current]|template_64[tmp_pos_mod_current];

						HashAd_first=arrayID-1;
						HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
						p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
						p_hv_tmp_length_cur++;
					}
				}
			}
		}
		else if(pos_label==2)
		{
			for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-extend_length)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				uint64_t tmp_pos_divide_current=(task_start_pos+j)/64;
				uint64_t tmp_pos_mod_current=(task_start_pos+j)%64;
				if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
				{
					int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
					if(arrayID!=0)
					{
						HashAd_first=arrayID-1;
						HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
						p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
						p_hv_tmp_length_cur++;
					}
				}
			}

		}
		else if(pos_label==1)
		{
			for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-extend_length)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
				if(arrayID!=0)
				{
					uint64_t tmp_pos_divide=(task_start_pos+j)/64;
					uint64_t tmp_pos_mod=(task_start_pos+j)%64;
					p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];

					HashAd_first=arrayID-1;
					HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
					p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
					p_hv_tmp_length_cur++;
				}
			}
		}
		else
		{
			for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-extend_length)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
				if(arrayID!=0)
				{
					HashAd_first=arrayID-1;
					HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
					p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
					p_hv_tmp_length_cur++;
				}
			}
		}
	}
	else
	{
		uint64_t task_start_pos_out;
		if(task_start_pos<extend_length)
		{
			task_start_pos_out=extend_length;
		}
		else
		{
			task_start_pos_out=task_start_pos;
		}
		cal_hash_value_directly_256bit(temp_seq+task_start_pos_out-task_start_pos,seq_k_current,para);
		int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
		if(arrayID!=0)
		{
			if(pos_label==1||pos_label==3)
			{
				uint64_t tmp_pos_divide=(task_start_pos_out+original_kmer_length-1)/64;
				uint64_t tmp_pos_mod=(task_start_pos_out+original_kmer_length-1)%64;
				p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];
			}

			HashAd_first=arrayID-1;
			HashAd_second=cal_hash_value_directly(temp_seq+task_start_pos_out-task_start_pos-extend_length,extend_length);//ycy20180418am1056
			p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos_out;
			p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
			p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
			p_hv_tmp_length_cur++;
		}
		if(pos_label==3)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-1)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				uint64_t tmp_pos=task_start_pos+j+original_kmer_length-1;
				uint64_t tmp_pos_divide_current=tmp_pos/64;
				uint64_t tmp_pos_mod_current=tmp_pos%64;

				if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
				{
					int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
					if(arrayID!=0)
					{
						p_candidate_pos[tmp_pos_divide_current]=p_candidate_pos[tmp_pos_divide_current]|template_64[tmp_pos_mod_current];

						HashAd_first=arrayID-1;
						HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
						p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
						p_hv_tmp_length_cur++;
					}
				}
			}
		}
		else if(pos_label==2)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-1)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				uint64_t tmp_pos=task_start_pos+j+original_kmer_length-1;
				uint64_t tmp_pos_divide_current=tmp_pos/64;
				uint64_t tmp_pos_mod_current=tmp_pos%64;

				if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
				{
					int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
					if(arrayID!=0)
					{
						HashAd_first=arrayID-1;
						HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
						p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
						p_hv_tmp_length_cur++;
					}
				}
			}
		}
		else if(pos_label==1)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-1)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
				if(arrayID!=0)
				{
					uint64_t tmp_pos=task_start_pos+j+original_kmer_length-1;
					uint64_t tmp_pos_divide_current=tmp_pos/64;
					uint64_t tmp_pos_mod_current=tmp_pos%64;
					p_candidate_pos[tmp_pos_divide_current]=p_candidate_pos[tmp_pos_divide_current]|template_64[tmp_pos_mod_current];

					HashAd_first=arrayID-1;
					HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
					p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
					p_hv_tmp_length_cur++;
				}
			}
		}
		else
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
			{
				if(j+task_start_pos>seq_length-original_kmer_length-1)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);//ycy20180418am1056

				int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
				if(arrayID!=0)
				{
					HashAd_first=arrayID-1;
					HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					p_hv_tmp[p_hv_tmp_length_cur].hv1=HashAd_first;
					p_hv_tmp[p_hv_tmp_length_cur].hv2=HashAd_second;
					p_hv_tmp_length_cur++;
				}
			}
		}
	}
	*(tmp->result_len)=p_hv_tmp_length_cur;
	return NULL;
}
void *Assign_BplusTree_mode(void *arg)
{
	struct Parallel_para_extend04_mode *tmp=(struct Parallel_para_extend04_mode *)arg;
	char *seq=tmp->seq;
	uint8_t **HashAd=tmp->HashAd;
	uint32_t kmer_len1=tmp->kmer_len1;
	uint32_t InOut_label=tmp->InOut_label;
	uint64_t *result_len=tmp->result_len;
	struct TwoDimenHashV** p_hv_tmp=tmp->p_hv_tmp;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;

	uint32_t original_kmer_length=kmer_len1;

	if(InOut_label==0)
	{
		for(uint32_t i=0;i<thread_num;i++)
		{
			for(uint64_t k=0;k<result_len[i];k++)
			{
				if(p_hv_tmp[i][k].hv1%thread_num==thread_id)
				{
					uint64_t HashAd_second_4bit=p_hv_tmp[i][k].hv2>>1;
					uint64_t HashAd_second_4bit_label=p_hv_tmp[i][k].hv2&1;
					Transition_4bit(seq[p_hv_tmp[i][k].pos-1],HashAd[HashAd_second_4bit]+\
							p_hv_tmp[i][k].hv1,HashAd_second_4bit_label);
				}
			}
		}
	}
	else
	{
		for(uint32_t i=0;i<thread_num;i++)
		{
			for(uint64_t k=0;k<result_len[i];k++)
			{
				if(p_hv_tmp[i][k].hv1%thread_num==thread_id)
				{
					uint64_t HashAd_second_4bit=p_hv_tmp[i][k].hv2>>1;
					uint64_t HashAd_second_4bit_label=p_hv_tmp[i][k].hv2&1;
					Transition_4bit(seq[p_hv_tmp[i][k].pos+original_kmer_length],HashAd[HashAd_second_4bit]+\
							p_hv_tmp[i][k].hv1,HashAd_second_4bit_label);
				}
			}
		}
	}
	return NULL;
}
void *judge_thread_in_HashTable(void *arg)
{
	struct reslt_para_plustree_hashtable* tmp=(struct reslt_para_plustree_hashtable*)arg;

	uint64_t start=tmp->start;
	uint64_t end=tmp->end;
	uint64_t reslt_Max_size;
	uint32_t ad_label=tmp->ad_label;
	struct bit256KmerPara para1;
	struct bit256KmerPara para2;
	para1=tmp->para1;
	para2=tmp->para2;

	uint64_t *reslt;
	uint8_t *reslt_c;

	uint8_t **HashAd=tmp->HashAd;
	uint32_t extend_length=tmp->extend_length;
	uint64_t each_hashTable_size=tmp->each_hashTable_size;
	uint64_t *h1=tmp->h1;
	uint32_t reslt_size=0;
	reslt=(uint64_t *)malloc(sizeof(uint64_t)*para2.kmer64Len*each_hashTable_size*2);
	reslt_c=(uint8_t *)malloc(sizeof(uint8_t)*each_hashTable_size*2);
	reslt_Max_size=each_hashTable_size*2;
	for(uint64_t k=start;k<=end;k++)
	{
		uint64_t* h1_temp=h1+k*para1.kmer64Len;
		for(uint64_t j=0;j<each_hashTable_size;j++)
		{
			if(Judge_4bit(HashAd[j][k],0)==1)
			{
				if(reslt_size<reslt_Max_size)
				{
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerIn(h1_temp,p_reslt_current,j<<1,extend_length, para1);
					if(ad_label==1)
					{
						uint8_t c_tmp=(HashAd[j][k]>>4)&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
				else
				{
					reslt=(uint64_t *)realloc(reslt,sizeof(uint64_t)*para2.kmer64Len*(reslt_Max_size+pow(2,10)));
					if(ad_label==1)
					{
						reslt_c=(uint8_t*)realloc(reslt_c,sizeof(uint8_t)*(reslt_Max_size+pow(2,10)));
					}
					reslt_Max_size=reslt_Max_size+pow(2,10);

					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerIn(h1_temp,p_reslt_current,j<<1,extend_length, para1);
					//reslt[reslt_size]=(h1_temp<<(2*extend_length))+(j<<1);
					if(ad_label==1)
					{
						uint8_t c_tmp=(HashAd[j][k]>>4)&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
			}
			if(Judge_4bit(HashAd[j][k],1)==1)
			{
				if(reslt_size<reslt_Max_size)
				{
					//reslt[reslt_size]=(h1_temp<<(2*extend_length))+(j<<1)+1;
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerIn(h1_temp,p_reslt_current,(j<<1)+1,extend_length, para1);
					if(ad_label==1)
					{
						uint8_t c_tmp=HashAd[j][k]&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
				else
				{
					reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*para2.kmer64Len*(reslt_Max_size+pow(2,10)));
					if(ad_label==1)
					{
						reslt_c=(uint8_t*)realloc(reslt_c,sizeof(uint8_t)*(reslt_Max_size+pow(2,10)));
					}
					reslt_Max_size=reslt_Max_size+pow(2,10);

					//reslt[reslt_size]=(h1_temp<<(2*extend_length))+(j<<1)+1;
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerIn(h1_temp,p_reslt_current,(j<<1)+1,extend_length, para1);
					if(ad_label==1)
					{
						uint8_t c_tmp=HashAd[j][k]&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
			}
		}
	}

	*(tmp->reslt_size)=reslt_size;
	*(tmp->reslt)=reslt;
	*(tmp->reslt_c)=reslt_c;
	return NULL;
}
void *judge_thread_out_HashTable(void *arg)
{
	struct reslt_para_plustree_hashtable* tmp=(struct reslt_para_plustree_hashtable*)arg;

	uint64_t start=tmp->start;
	uint64_t end=tmp->end;
	uint32_t ad_label=tmp->ad_label;
	uint64_t reslt_Max_size;
	struct bit256KmerPara para1;
	para1=tmp->para1;
	struct bit256KmerPara para2;
	para2=tmp->para2;

//	cout << para1.remainer1to64<<endl;
//	cout << para1.kmer1Len<<endl;
//	cout << para1.kmer64Len<<endl;
//	cout << para2.remainer1to64<<endl;
//	cout << para2.kmer1Len<<endl;
//	cout << para2.kmer64Len<<endl;
	uint64_t *reslt;
	uint8_t *reslt_c;

	uint8_t **HashAd=tmp->HashAd;
	uint32_t extend_length=tmp->extend_length;
	uint64_t each_hashTable_size=tmp->each_hashTable_size;
	uint64_t *h1=tmp->h1;
	uint32_t reslt_size=0;
	reslt=(uint64_t *)malloc(sizeof(uint64_t)*each_hashTable_size*2*para2.kmer64Len);
	if(ad_label==1)
	{
		reslt_c=(uint8_t *)malloc(sizeof(uint8_t)*each_hashTable_size*2);
	}
	reslt_Max_size=each_hashTable_size*2;

	for(uint64_t k=start;k<=end;k++)
	{
		uint64_t *h1_temp;
		h1_temp=h1+k*para1.kmer64Len;
		for(uint64_t j=0;j<each_hashTable_size;j++)
		{
			if(Judge_4bit(HashAd[j][k],0)==1)
			{
				if(reslt_size<reslt_Max_size)
				{
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerOut(h1_temp,p_reslt_current,j<<1,extend_length,para1);
					//reslt[reslt_size]=h1_temp+((j<<1)<<(2*extend_length));
					if(ad_label==1)
					{
						uint8_t c_tmp=(HashAd[j][k]>>4)&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
				else
				{
					reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*para2.kmer64Len*(reslt_Max_size+pow(2,10)));
					if(ad_label==1)
					{
						reslt_c=(uint8_t*)realloc(reslt_c,sizeof(uint8_t)*(reslt_Max_size+pow(2,10)));
					}
					reslt_Max_size=reslt_Max_size+pow(2,10);

					//reslt[reslt_size]=h1_temp+((j<<1)<<(2*extend_length));
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerOut(h1_temp,p_reslt_current,j<<1,extend_length,para1);
					if(ad_label==1)
					{
						uint8_t c_tmp=(HashAd[j][k]>>4)&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
			}
			if(Judge_4bit(HashAd[j][k],1)==1)
			{
				if(reslt_size<reslt_Max_size)
				{
					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerOut(h1_temp,p_reslt_current,(j<<1)+1,extend_length,para1);
					//reslt[reslt_size]=h1_temp+(((j<<1)+1)<<(2*extend_length));
					if(ad_label==1)
					{
						uint8_t c_tmp=HashAd[j][k]&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
				else
				{
					reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*para2.kmer64Len*(reslt_Max_size+pow(2,10)));
					if(ad_label==1)
					{
						reslt_c=(uint8_t*)realloc(reslt_c,sizeof(uint8_t)*(reslt_Max_size+pow(2,10)));
					}
					reslt_Max_size=reslt_Max_size+pow(2,10);

					uint64_t *p_reslt_current=reslt+reslt_size*para2.kmer64Len;
					cal_newExtendedKmerOut(h1_temp,p_reslt_current,(j<<1)+1,extend_length,para1);
					//reslt[reslt_size]=h1_temp+(((j<<1)+1)<<(2*extend_length));
					if(ad_label==1)
					{
						uint8_t c_tmp=HashAd[j][k]&15;
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
			}
		}
	}

	*(tmp->reslt_size)=reslt_size;
	*(tmp->reslt)=reslt;
	if(ad_label==1)
	{
		*(tmp->reslt_c)=reslt_c;
	}
	return NULL;
}
void *create_HashTabel(void *arg)
{
	struct create_HashTable_para* tmp=(struct create_HashTable_para*)arg;

	struct bit256Hash* p_root=tmp->p_root;
	uint64_t *h1=tmp->h1;
	uint32_t thread_num=tmp->thread_num;
	uint32_t thread_id=tmp->thread_id;
	uint64_t Tasksize_current=tmp->Tasksize_current;
	struct bit256KmerPara para1=tmp->para1;

	for(uint64_t i=0;i<Tasksize_current;i++)
	{
		insert256BitHashTable_parallel(p_root,h1+i*para1.kmer64Len,i+1,para1,thread_num,thread_id);
	}
	return NULL;
}
void Cal_Single_Task_4bit_UsingHashTable_withLabel_256bitKmer(struct para_bPlusTreeMap tmp)
{
	uint32_t T_i=tmp.T_i;
	char *seq;
	char* dataset=tmp.seq;
	uint64_t seq_length;
	uint64_t start_line=tmp.start_line;
	uint64_t *end_line=tmp.end_line;
	uint64_t Tasksize=tmp.Tasksize;
	uint32_t extend_length=tmp.extend_length;
	uint32_t original_kmer_length=tmp.original_kmer_length;
	uint32_t InOut_label=tmp.InOut_label;
	uint32_t ad_label=tmp.ad_label;
	uint32_t thread_num=tmp.thread_num;
	uint64_t Task_Size_Single_Thread=tmp.Task_Size_Single_Thread;
	uint32_t each_hashTable_size=pow(4,extend_length)/2;
	uint32_t rc_label=tmp.rc_label;
	Task_Size_Single_Thread=pow(2,Task_Size_Single_Thread);
	char*inputFileName=tmp.inputFileName;
	cout << extend_length << "	" << original_kmer_length << endl;

	cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	cerr << "start the " << T_i << "th group" << endl;

	struct bit256KmerPara para1;
	para1.kmer1Len=original_kmer_length*2;
	cout <<para1.kmer1Len <<endl;
	para1.remainer1to64=para1.kmer1Len%64;
	para1.kmer64Len=para1.kmer1Len/64+(para1.remainer1to64?1:0);
	cout <<para1.kmer64Len <<endl;
	if(para1.remainer1to64==0)
	{
		para1.codefor1=0;
	}
	else
	{
		para1.codefor1=0;
		for(uint32_t i=0;i<para1.remainer1to64-1;i++)
		{
			para1.codefor1=para1.codefor1|1;
			para1.codefor1=para1.codefor1<<1;
		}
		para1.codefor1=para1.codefor1|1;
	}

	cout << extend_length << "	" << original_kmer_length << endl;

	struct bit256KmerPara para2;
	para2.kmer1Len=(original_kmer_length+extend_length)*2;
	cout <<para2.kmer1Len <<endl;
	para2.remainer1to64=para2.kmer1Len%64;
	para2.kmer64Len=para2.kmer1Len/64+(para2.remainer1to64?1:0);
	cout <<para2.kmer64Len <<endl;
	if(para2.remainer1to64==0)
	{
		para2.codefor1=0;
	}
	else
	{
		para2.codefor1=0;
		for(uint32_t i=0;i<para2.remainer1to64-1;i++)
		{
			para2.codefor1=para2.codefor1|1;
			para2.codefor1=para2.codefor1<<1;
		}
		para2.codefor1=para2.codefor1|1;
	}


	cout << extend_length << "	" << original_kmer_length << endl;
	/**********************************read kmers***********************************/
	FILE * f_in;
	f_in=fopen(inputFileName,"rb");
	uint64_t *h1;
	h1=(uint64_t*)malloc(sizeof(uint64_t)*Tasksize*para1.kmer64Len);
	fseek(f_in,(start_line)*sizeof(uint64_t)*para1.kmer64Len,0);
	uint64_t Tasksize_current=fread(h1,sizeof(uint64_t),para1.kmer64Len*Tasksize,f_in);
	Tasksize_current=Tasksize_current/para1.kmer64Len;
	*end_line=start_line+Tasksize_current;
	cout << "read task over!" << endl;
	cout <<"the task size is: "<< Tasksize_current<<endl;
	if(Tasksize_current<=0)
	{
		cout << "no more task!" << endl;
		return;
	}
	fclose(f_in);
	/**********************************read kmers***********************************/

	/**********************************alloc hashad space***************************/
	uint8_t **HashAd;
	HashAd=(uint8_t**)malloc(sizeof(uint8_t*)*each_hashTable_size);
	for(uint64_t i=0;i<each_hashTable_size;i++)
	{
		HashAd[i]=(uint8_t*)malloc(sizeof(uint8_t)*Tasksize_current);
		memset(HashAd[i],0,sizeof(uint8_t)*Tasksize_current);
	}
	/**********************************alloc hashad space***************************/

	//hash the first hash value to the first dimensional coordinate of HashAd
	struct bit256Hash* p_root=NULL;
	p_root=initial256BitHashTable();
	if(p_root==NULL)
	{
		cerr<< "error!" << endl;
	}
	struct timeval tvs_hs,tve_hs;
	gettimeofday(&tvs_hs,NULL);
//	for(uint64_t i=0;i<Tasksize_current;i++)
//	{
//		insert256BitHashTable(p_root,h1+i*para1.kmer64Len,i+1,para1);
//	}
	if(thread_num==1)
	{
		for(uint64_t i=0;i<Tasksize_current;i++)
		{
			insert256BitHashTable(p_root,h1+i*para1.kmer64Len,i+1,para1);
		}
	}
	else
	{
		pthread_t* t_hash;
		t_hash=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
		struct create_HashTable_para* p;
		p=(struct create_HashTable_para*)malloc(sizeof(struct create_HashTable_para)*thread_num);
		for(uint32_t i=0;i<thread_num;i++)
		{
			p[i].p_root=p_root;
			p[i].h1=h1;
			p[i].thread_num=thread_num;
			p[i].thread_id=i;
			p[i].Tasksize_current=Tasksize_current;
			p[i].para1=para1;
		}
		for(uint32_t j=0;j<thread_num;j++)
		{
			if(pthread_create(t_hash+j, NULL, create_HashTabel, (void*)(p+j))!=0)
			{
				cout << "error!" << endl;
			}
		}
		for(uint32_t i=0;i<thread_num;i++)
		{
			pthread_join(t_hash[i], NULL);
		}
		free(t_hash);
		free(p);
	}
	gettimeofday(&tve_hs,NULL);
	double span_hs = tve_hs.tv_sec-tvs_hs.tv_sec + (tve_hs.tv_usec-tvs_hs.tv_usec)/1000000.0;
	cout << "load hash time is: "<<span_hs<<endl;
	cerr << "hash kmer over!" << endl;

	//test
	struct timeval tvs_hs_c,tve_hs_c;
	gettimeofday(&tvs_hs_c,NULL);
	for(uint64_t i=0;i<Tasksize_current;i++)
	{
		uint64_t tmp;
		tmp=search256BitHashTable(p_root,h1+i*para1.kmer64Len,para1);
		if(tmp!=i+1)
		{
			cout << "error!" << endl;
		}
	}
	gettimeofday(&tve_hs_c,NULL);
	double span_hs_c = tve_hs_c.tv_sec-tvs_hs_c.tv_sec + (tve_hs_c.tv_usec-tvs_hs_c.tv_usec)/1000000.0;
	cout << "check hash time is: "<<span_hs_c<<endl;
	cerr << "check hash kmer over!" << endl;
	//hash the first hash value to the first dimensional coordinate of HashAd

	//fill the two-dimensional hash table!
	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		cout << ref_i+1 << ":"<< p_ref_path.NumberOfPathes <<endl;
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		char *temp_seq;
		temp_seq=seq;

		if(tmp.pos_label==1||tmp.pos_label==3)
		{
			if(T_i==0)
			{
				tmp.p_candidate_pos[ref_i]=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
				memset(tmp.p_candidate_pos[ref_i],0,sizeof(uint64_t)*(seq_length/64+1));
				tmp.p_candidate_pos_len[ref_i]=seq_length/64+1;
			}
			else
			{
				tmp.p_candidate_pos[ref_i]=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
				memset(tmp.p_candidate_pos[ref_i],0,sizeof(uint64_t)*(seq_length/64+1));
				tmp.p_candidate_pos_len[ref_i]=seq_length/64+1;

				char* posfilename;
				posfilename=(char*)malloc(sizeof(char)*14);
				struct para_getN tmp_inName;
				tmp_inName.kmerlen=original_kmer_length+extend_length;
				tmp_inName.FileName=posfilename;
				tmp_inName.InOut_label=InOut_label;
				tmp_inName.isMiddle=0;
				tmp_inName.isposition=ref_i+1;
				getFileName(tmp_inName);

				FILE* fp_pos;
				fp_pos=fopen(posfilename,"rb");
				fread(tmp.p_candidate_pos[ref_i],sizeof(uint64_t),tmp.p_candidate_pos_len[ref_i],fp_pos);
				fclose(fp_pos);
				free(posfilename);
			}
		}

		uint64_t *p_candidate_current=NULL;
		if(tmp.pos_label==2||tmp.pos_label==3)
		{
			char* posfilename;
			posfilename=(char*)malloc(sizeof(char)*14);
			struct para_getN tmp_inName;
			tmp_inName.kmerlen=original_kmer_length;
			tmp_inName.FileName=posfilename;
			tmp_inName.InOut_label=InOut_label;
			tmp_inName.isMiddle=0;
			tmp_inName.isposition=ref_i+1;
			getFileName(tmp_inName);

			FILE* fp_pos;
			fp_pos=fopen(posfilename,"rb");
			p_candidate_current=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
			uint64_t x=fread(p_candidate_current,sizeof(uint64_t),seq_length/64+1,fp_pos);
			if(x!=seq_length/64+1)
			{
				cout << "error: read p_candidate_current failed!" << endl;
			}
			fclose(fp_pos);
		}

		struct timeval tvs,tve;
		gettimeofday(&tvs,NULL);
		uint64_t HashAd_first;
		uint64_t HashAd_second;
		uint64_t seq_k_current[4];
		if(InOut_label==0)
		{
			if(thread_num>1)
			{
				//2)
				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
				uint64_t unit_long=Task_Size_Single_Thread;
				struct TwoDimenHashV** p_hv_tmp;
				p_hv_tmp=(struct TwoDimenHashV**)malloc(sizeof(struct TwoDimenHashV*)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					p_hv_tmp[i]=(struct TwoDimenHashV*)malloc(sizeof(struct TwoDimenHashV)*unit_long);
				}
				uint64_t* result_len;
				result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
				struct Parallel_para_extend03_hashtable* para_alloc;
				para_alloc=(struct Parallel_para_extend03_hashtable*)malloc(sizeof(struct Parallel_para_extend03_hashtable)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_alloc[i].seq=seq;
					para_alloc[i].seq_length=seq_length;
					para_alloc[i].p_root=p_root;
					para_alloc[i].kmer_len1=original_kmer_length;
					para_alloc[i].kmer_len2=extend_length;
					para_alloc[i].task_len=unit_long;
					para_alloc[i].InOut_label=InOut_label;
					para_alloc[i].result_len=result_len+i;
					para_alloc[i].p_hv_tmp=p_hv_tmp[i];
					para_alloc[i].para=para1;

					para_alloc[i].p_candidate_pos=tmp.p_candidate_pos[ref_i];
					para_alloc[i].template_64=tmp.template_64;
					para_alloc[i].pos_label=tmp.pos_label;
					para_alloc[i].p_candidate_current=p_candidate_current;

				}

				struct Parallel_para_extend04_mode* para_assi;
				para_assi=(struct Parallel_para_extend04_mode*)malloc(sizeof(struct Parallel_para_extend04_mode)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_assi[i].seq=seq;
					para_assi[i].HashAd=HashAd;
					para_assi[i].kmer_len1=original_kmer_length;
					para_assi[i].kmer_len2=extend_length;
					para_assi[i].InOut_label=InOut_label;
					para_assi[i].result_len=result_len;
					para_assi[i].p_hv_tmp=p_hv_tmp;
					para_assi[i].thread_id=i;
					para_assi[i].thread_num=thread_num;
				}

				uint64_t start=0;
				for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length;loop_num++)
				{
					start=loop_num*(thread_num*unit_long);
					for(uint32_t j=0;j<thread_num;j++)
					{
						*(para_alloc[j].result_len)=0;//19-01-15
						para_alloc[j].task_start_pos=start+j*unit_long;
						if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length+extend_length))
						{
							if(start+(j+1)*unit_long>seq_length-(original_kmer_length+extend_length))
							{
								para_alloc[j].task_len=seq_length-(original_kmer_length+extend_length)-para_alloc[j].task_start_pos+1;
							}
							else
							{
								para_alloc[j].task_len=unit_long;
							}
						}
					}
					for(uint32_t j=0;j<thread_num;j++)
					{
						if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length+extend_length))
						{
							if(pthread_create(t+j, NULL, Allocate_HashTable, (void*)(para_alloc+j))!=0)
							{
								cout << "error!" << endl;
							}
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length+extend_length))
						{
							pthread_join(t[i], NULL);
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(pthread_create(t+i, NULL, Assign_BplusTree_mode, (void*)(para_assi+i))!=0)
						{
							cout << "error!" << endl;
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						pthread_join(t[i], NULL);
					}
				}
				free(t);
				for(uint32_t i=0;i<thread_num;i++)
				{
					free(p_hv_tmp[i]);
				}
				free(p_hv_tmp);
				free(result_len);
				free(para_alloc);
				free(para_assi);

			}
			else
			{
				cal_hash_value_directly_256bit(temp_seq,seq_k_current,para1);
				if(tmp.pos_label==3)
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{

						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						uint64_t tmp_pos_divide_current=j/64;
						uint64_t tmp_pos_mod_current=j%64;
						if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
						{
							int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
							if(arrayID!=0)
							{
								tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]=tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]|tmp.template_64[tmp_pos_mod_current];

								HashAd_first=arrayID-1;
								HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
								uint64_t HashAd_second_4bit=HashAd_second>>1;
								uint64_t HashAd_second_4bit_label=HashAd_second&1;
								Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
										HashAd_first,HashAd_second_4bit_label);
							}
						}
					}
				}
				else if(tmp.pos_label==2)
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{

						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						uint64_t tmp_pos_divide_current=j/64;
						uint64_t tmp_pos_mod_current=j%64;
						if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
						{
							int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
							if(arrayID!=0)
							{
								HashAd_first=arrayID-1;
								HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
								uint64_t HashAd_second_4bit=HashAd_second>>1;
								uint64_t HashAd_second_4bit_label=HashAd_second&1;
								Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
										HashAd_first,HashAd_second_4bit_label);
							}
						}
					}
				}
				else if(tmp.pos_label==1)
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
						if(arrayID!=0)
						{
							uint64_t tmp_pos_divide=j/64;
							uint64_t tmp_pos_mod=j%64;
							tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];
							HashAd_first=arrayID-1;
							HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
							uint64_t HashAd_second_4bit=HashAd_second>>1;
							uint64_t HashAd_second_4bit_label=HashAd_second&1;
							Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
									HashAd_first,HashAd_second_4bit_label);
						}
					}
				}
				else
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
						if(arrayID!=0)
						{
							HashAd_first=arrayID-1;
							HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
							uint64_t HashAd_second_4bit=HashAd_second>>1;
							uint64_t HashAd_second_4bit_label=HashAd_second&1;
							Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
									HashAd_first,HashAd_second_4bit_label);
						}
					}
				}
			}
		}
		else
		{
			if(thread_num>1)
			{
				//2)
				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
				uint64_t unit_long=Task_Size_Single_Thread;
				struct TwoDimenHashV** p_hv_tmp;
				p_hv_tmp=(struct TwoDimenHashV**)malloc(sizeof(struct TwoDimenHashV*)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					p_hv_tmp[i]=(struct TwoDimenHashV*)malloc(sizeof(struct TwoDimenHashV)*unit_long);
				}
				uint64_t* result_len;
				result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
				struct Parallel_para_extend03_hashtable* para_alloc;
				para_alloc=(struct Parallel_para_extend03_hashtable*)malloc(sizeof(struct Parallel_para_extend03_hashtable)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_alloc[i].seq=seq;
					para_alloc[i].seq_length=seq_length;
					para_alloc[i].p_root=p_root;
					para_alloc[i].kmer_len1=original_kmer_length;
					para_alloc[i].kmer_len2=extend_length;
					para_alloc[i].task_len=unit_long;
					para_alloc[i].InOut_label=InOut_label;
					para_alloc[i].result_len=result_len+i;
					para_alloc[i].p_hv_tmp=p_hv_tmp[i];
					para_alloc[i].para=para1;

					para_alloc[i].p_candidate_pos=tmp.p_candidate_pos[ref_i];
					para_alloc[i].template_64=tmp.template_64;
					para_alloc[i].p_candidate_current=p_candidate_current;
					para_alloc[i].pos_label=tmp.pos_label;
				}
				struct Parallel_para_extend04_mode* para_assi;
				para_assi=(struct Parallel_para_extend04_mode*)malloc(sizeof(struct Parallel_para_extend04_mode)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_assi[i].seq=seq;
					para_assi[i].HashAd=HashAd;
					para_assi[i].kmer_len1=original_kmer_length;
					para_assi[i].kmer_len2=extend_length;
					para_assi[i].InOut_label=InOut_label;
					para_assi[i].result_len=result_len;
					para_assi[i].p_hv_tmp=p_hv_tmp;
					para_assi[i].thread_id=i;
					para_assi[i].thread_num=thread_num;
				}

				uint64_t start=0;
				for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length;loop_num++)
				{
					start=loop_num*thread_num*unit_long;
					for(uint32_t j=0;j<thread_num;j++)
					{
						*(para_alloc[j].result_len)=0;//19-01-15
						para_alloc[j].task_start_pos=start+j*unit_long;
						if(para_alloc[j].task_start_pos<=seq_length-original_kmer_length-1)
						{
							if(start+(j+1)*unit_long>seq_length-original_kmer_length-1)
							{
								para_alloc[j].task_len=seq_length-original_kmer_length-1-para_alloc[j].task_start_pos+1;
							}
							else
							{
								para_alloc[j].task_len=unit_long;
							}
						}
					}
					for(uint32_t j=0;j<thread_num;j++)
					{
						if(para_alloc[j].task_start_pos<=seq_length-original_kmer_length-1)
						{
							if(pthread_create(t+j, NULL, Allocate_HashTable, (void*)(para_alloc+j))!=0)
							{
								cout << "error!" << endl;
							}
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(para_alloc[i].task_start_pos<=seq_length-original_kmer_length-1)
						{
							pthread_join(t[i], NULL);
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(pthread_create(t+i, NULL, Assign_BplusTree_mode, (void*)(para_assi+i))!=0)
						{
							cout << "error!" << endl;
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						pthread_join(t[i], NULL);
					}
				}
				free(t);
				for(uint32_t i=0;i<thread_num;i++)
				{
					free(p_hv_tmp[i]);
				}
				free(p_hv_tmp);
				free(result_len);
				free(para_alloc);
				free(para_assi);
			}
			else
			{
				cal_hash_value_directly_256bit(temp_seq+extend_length,seq_k_current,para1);
				int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
				if(arrayID!=0)
				{
					if(tmp.pos_label==1||tmp.pos_label==3)
					{
						uint64_t tmp_pos=extend_length+original_kmer_length-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];
					}

					HashAd_first=arrayID-1;
					HashAd_second=cal_hash_value_directly(temp_seq,extend_length);//ycy20180418am1056
					uint64_t HashAd_second_4bit=HashAd_second>>1;
					uint64_t HashAd_second_4bit_label=HashAd_second&1;
					Transition_4bit(seq[extend_length+original_kmer_length],HashAd[HashAd_second_4bit]+\
							HashAd_first,HashAd_second_4bit_label);
				}
				if(tmp.pos_label==3)
				{
					for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056
						uint64_t tmp_pos=j+original_kmer_length-1;
						uint64_t tmp_pos_divide_current=tmp_pos/64;
						uint64_t tmp_pos_mod_current=tmp_pos%64;
						uint64_t tmp_template=tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current];
						if(tmp.template_64[tmp_pos_mod_current]==tmp_template)
						{
							arrayID = search256BitHashTable(p_root,seq_k_current,para1);
							if(arrayID!=0)
							{
								tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]=tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]|tmp.template_64[tmp_pos_mod_current];

								HashAd_first=arrayID-1;
								HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
								uint64_t HashAd_second_4bit=HashAd_second>>1;
								uint64_t HashAd_second_4bit_label=HashAd_second&1;
								Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
										HashAd_first,HashAd_second_4bit_label);
							}
						}
					}
				}
				else if(tmp.pos_label==2)
				{
					for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						uint64_t tmp_pos=j+original_kmer_length-1;
						uint64_t tmp_pos_divide_current=tmp_pos/64;
						uint64_t tmp_pos_mod_current=tmp_pos%64;
						if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
						{
							arrayID = search256BitHashTable(p_root,seq_k_current,para1);
							if(arrayID!=0)
							{
								HashAd_first=arrayID-1;
								HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
								uint64_t HashAd_second_4bit=HashAd_second>>1;
								uint64_t HashAd_second_4bit_label=HashAd_second&1;
								Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
										HashAd_first,HashAd_second_4bit_label);
							}
						}
					}
				}
				else if(tmp.pos_label==1)
				{
					for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						arrayID = search256BitHashTable(p_root,seq_k_current,para1);
						if(arrayID!=0)
						{
							uint64_t tmp_pos=j+original_kmer_length-1;
							uint64_t tmp_pos_divide=tmp_pos/64;
							uint64_t tmp_pos_mod=tmp_pos%64;
							tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];

							HashAd_first=arrayID-1;
							HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
							uint64_t HashAd_second_4bit=HashAd_second>>1;
							uint64_t HashAd_second_4bit_label=HashAd_second&1;
							Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
									HashAd_first,HashAd_second_4bit_label);
						}
					}
				}
				else
				{
					for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

						arrayID = search256BitHashTable(p_root,seq_k_current,para1);
						if(arrayID!=0)
						{
							HashAd_first=arrayID-1;
							HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
							uint64_t HashAd_second_4bit=HashAd_second>>1;
							uint64_t HashAd_second_4bit_label=HashAd_second&1;
							Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
									HashAd_first,HashAd_second_4bit_label);
						}
					}
				}
			}
		}
		if(tmp.pos_label==2||tmp.pos_label==3)
		{
			free(p_candidate_current);
		}
		if(tmp.pos_label==1||tmp.pos_label==3)
		{
			char* posfilename;
			posfilename=(char*)malloc(sizeof(char)*14);
			struct para_getN tmp_inName;
			tmp_inName.kmerlen=original_kmer_length+extend_length;
			tmp_inName.FileName=posfilename;
			tmp_inName.InOut_label=InOut_label;
			tmp_inName.isMiddle=0;
			tmp_inName.isposition=ref_i+1;
			getFileName(tmp_inName);

			FILE* fp_pos;
			fp_pos=fopen(posfilename,"wb");
			fwrite(tmp.p_candidate_pos[ref_i],sizeof(uint64_t),tmp.p_candidate_pos_len[ref_i],fp_pos);
			fclose(fp_pos);
			free(posfilename);
			free(tmp.p_candidate_pos[ref_i]);
			tmp.p_candidate_pos[ref_i]=NULL;
		}
		if(rc_label==1)//20200207
        {
            cout << ref_i+1 << "RC:"<< p_ref_path.NumberOfPathes <<endl;
            rc(&seq,seq_length);
            char *temp_seq;
            temp_seq=seq;

            if(tmp.pos_label==1||tmp.pos_label==3)
            {
                if(T_i==0)
                {
                    tmp.p_candidate_pos[ref_i]=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
                    memset(tmp.p_candidate_pos[ref_i],0,sizeof(uint64_t)*(seq_length/64+1));
                    tmp.p_candidate_pos_len[ref_i]=seq_length/64+1;
                }
                else
                {
                    tmp.p_candidate_pos[ref_i]=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
                    memset(tmp.p_candidate_pos[ref_i],0,sizeof(uint64_t)*(seq_length/64+1));
                    tmp.p_candidate_pos_len[ref_i]=seq_length/64+1;

                    char* posfilename;
                    posfilename=(char*)malloc(sizeof(char)*14);
                    struct para_getN tmp_inName;
                    tmp_inName.kmerlen=original_kmer_length+extend_length;
                    tmp_inName.FileName=posfilename;
                    tmp_inName.InOut_label=InOut_label;
                    tmp_inName.isMiddle=0;
                    tmp_inName.isposition=ref_i+1;
                    getFileName(tmp_inName);
                    posfilename[2]='r';
                    posfilename[3]='c';

                    FILE* fp_pos;
                    fp_pos=fopen(posfilename,"rb");
                    fread(tmp.p_candidate_pos[ref_i],sizeof(uint64_t),tmp.p_candidate_pos_len[ref_i],fp_pos);
                    fclose(fp_pos);
                    free(posfilename);
                }
            }

            uint64_t *p_candidate_current=NULL;
            if(tmp.pos_label==2||tmp.pos_label==3)
            {
                char* posfilename;
                posfilename=(char*)malloc(sizeof(char)*14);
                struct para_getN tmp_inName;
                tmp_inName.kmerlen=original_kmer_length;
                tmp_inName.FileName=posfilename;
                tmp_inName.InOut_label=InOut_label;
                tmp_inName.isMiddle=0;
                tmp_inName.isposition=ref_i+1;
                getFileName(tmp_inName);
                posfilename[2]='r';
                posfilename[3]='c';


                FILE* fp_pos;
                fp_pos=fopen(posfilename,"rb");
                p_candidate_current=(uint64_t*)malloc(sizeof(uint64_t)*(seq_length/64+1));
                uint64_t x=fread(p_candidate_current,sizeof(uint64_t),seq_length/64+1,fp_pos);
                if(x!=seq_length/64+1)
                {
                    cout << "error: read p_candidate_current failed!" << endl;
                }
                fclose(fp_pos);
            }

            struct timeval tvs,tve;
            gettimeofday(&tvs,NULL);
            uint64_t HashAd_first;
            uint64_t HashAd_second;
            uint64_t seq_k_current[4];
            if(InOut_label==0)
            {
                if(thread_num>1)
                {
                    //2)
                    pthread_t* t;
                    t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                    uint64_t unit_long=Task_Size_Single_Thread;
                    struct TwoDimenHashV** p_hv_tmp;
                    p_hv_tmp=(struct TwoDimenHashV**)malloc(sizeof(struct TwoDimenHashV*)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        p_hv_tmp[i]=(struct TwoDimenHashV*)malloc(sizeof(struct TwoDimenHashV)*unit_long);
                    }
                    uint64_t* result_len;
                    result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
                    struct Parallel_para_extend03_hashtable* para_alloc;
                    para_alloc=(struct Parallel_para_extend03_hashtable*)malloc(sizeof(struct Parallel_para_extend03_hashtable)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_alloc[i].seq=seq;
                        para_alloc[i].seq_length=seq_length;
                        para_alloc[i].p_root=p_root;
                        para_alloc[i].kmer_len1=original_kmer_length;
                        para_alloc[i].kmer_len2=extend_length;
                        para_alloc[i].task_len=unit_long;
                        para_alloc[i].InOut_label=InOut_label;
                        para_alloc[i].result_len=result_len+i;
                        para_alloc[i].p_hv_tmp=p_hv_tmp[i];
                        para_alloc[i].para=para1;

                        para_alloc[i].p_candidate_pos=tmp.p_candidate_pos[ref_i];
                        para_alloc[i].template_64=tmp.template_64;
                        para_alloc[i].pos_label=tmp.pos_label;
                        para_alloc[i].p_candidate_current=p_candidate_current;

                    }

                    struct Parallel_para_extend04_mode* para_assi;
                    para_assi=(struct Parallel_para_extend04_mode*)malloc(sizeof(struct Parallel_para_extend04_mode)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_assi[i].seq=seq;
                        para_assi[i].HashAd=HashAd;
                        para_assi[i].kmer_len1=original_kmer_length;
                        para_assi[i].kmer_len2=extend_length;
                        para_assi[i].InOut_label=InOut_label;
                        para_assi[i].result_len=result_len;
                        para_assi[i].p_hv_tmp=p_hv_tmp;
                        para_assi[i].thread_id=i;
                        para_assi[i].thread_num=thread_num;
                    }

                    uint64_t start=0;
                    for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length;loop_num++)
                    {
                        start=loop_num*(thread_num*unit_long);
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            *(para_alloc[j].result_len)=0;//19-01-15
                            para_alloc[j].task_start_pos=start+j*unit_long;
                            if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length+extend_length))
                            {
                                if(start+(j+1)*unit_long>seq_length-(original_kmer_length+extend_length))
                                {
                                    para_alloc[j].task_len=seq_length-(original_kmer_length+extend_length)-para_alloc[j].task_start_pos+1;
                                }
                                else
                                {
                                    para_alloc[j].task_len=unit_long;
                                }
                            }
                        }
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length+extend_length))
                            {
                                if(pthread_create(t+j, NULL, Allocate_HashTable, (void*)(para_alloc+j))!=0)
                                {
                                    cout << "error!" << endl;
                                }
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length+extend_length))
                            {
                                pthread_join(t[i], NULL);
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(pthread_create(t+i, NULL, Assign_BplusTree_mode, (void*)(para_assi+i))!=0)
                            {
                                cout << "error!" << endl;
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            pthread_join(t[i], NULL);
                        }
                    }
                    free(t);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        free(p_hv_tmp[i]);
                    }
                    free(p_hv_tmp);
                    free(result_len);
                    free(para_alloc);
                    free(para_assi);

                }
                else
                {
                    cal_hash_value_directly_256bit(temp_seq,seq_k_current,para1);
                    if(tmp.pos_label==3)
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {

                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            uint64_t tmp_pos_divide_current=j/64;
                            uint64_t tmp_pos_mod_current=j%64;
                            if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
                            {
                                int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                                if(arrayID!=0)
                                {
                                    tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]=tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]|tmp.template_64[tmp_pos_mod_current];

                                    HashAd_first=arrayID-1;
                                    HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
                                    uint64_t HashAd_second_4bit=HashAd_second>>1;
                                    uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                    Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
                                            HashAd_first,HashAd_second_4bit_label);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==2)
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {

                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            uint64_t tmp_pos_divide_current=j/64;
                            uint64_t tmp_pos_mod_current=j%64;
                            if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
                            {
                                int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                                if(arrayID!=0)
                                {
                                    HashAd_first=arrayID-1;
                                    HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
                                    uint64_t HashAd_second_4bit=HashAd_second>>1;
                                    uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                    Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
                                            HashAd_first,HashAd_second_4bit_label);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==1)
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                            if(arrayID!=0)
                            {
                                uint64_t tmp_pos_divide=j/64;
                                uint64_t tmp_pos_mod=j%64;
                                tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];
                                HashAd_first=arrayID-1;
                                HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
                                uint64_t HashAd_second_4bit=HashAd_second>>1;
                                uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
                                        HashAd_first,HashAd_second_4bit_label);
                            }
                        }
                    }
                    else
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                            if(arrayID!=0)
                            {
                                HashAd_first=arrayID-1;
                                HashAd_second=cal_hash_value_directly(temp_seq+j+original_kmer_length,extend_length);
                                uint64_t HashAd_second_4bit=HashAd_second>>1;
                                uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                Transition_4bit(seq[j-1],HashAd[HashAd_second_4bit]+\
                                        HashAd_first,HashAd_second_4bit_label);
                            }
                        }
                    }
                }
            }
            else
            {
                if(thread_num>1)
                {
                    //2)
                    pthread_t* t;
                    t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                    uint64_t unit_long=Task_Size_Single_Thread;
                    struct TwoDimenHashV** p_hv_tmp;
                    p_hv_tmp=(struct TwoDimenHashV**)malloc(sizeof(struct TwoDimenHashV*)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        p_hv_tmp[i]=(struct TwoDimenHashV*)malloc(sizeof(struct TwoDimenHashV)*unit_long);
                    }
                    uint64_t* result_len;
                    result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
                    struct Parallel_para_extend03_hashtable* para_alloc;
                    para_alloc=(struct Parallel_para_extend03_hashtable*)malloc(sizeof(struct Parallel_para_extend03_hashtable)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_alloc[i].seq=seq;
                        para_alloc[i].seq_length=seq_length;
                        para_alloc[i].p_root=p_root;
                        para_alloc[i].kmer_len1=original_kmer_length;
                        para_alloc[i].kmer_len2=extend_length;
                        para_alloc[i].task_len=unit_long;
                        para_alloc[i].InOut_label=InOut_label;
                        para_alloc[i].result_len=result_len+i;
                        para_alloc[i].p_hv_tmp=p_hv_tmp[i];
                        para_alloc[i].para=para1;

                        para_alloc[i].p_candidate_pos=tmp.p_candidate_pos[ref_i];
                        para_alloc[i].template_64=tmp.template_64;
                        para_alloc[i].p_candidate_current=p_candidate_current;
                        para_alloc[i].pos_label=tmp.pos_label;
                    }
                    struct Parallel_para_extend04_mode* para_assi;
                    para_assi=(struct Parallel_para_extend04_mode*)malloc(sizeof(struct Parallel_para_extend04_mode)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_assi[i].seq=seq;
                        para_assi[i].HashAd=HashAd;
                        para_assi[i].kmer_len1=original_kmer_length;
                        para_assi[i].kmer_len2=extend_length;
                        para_assi[i].InOut_label=InOut_label;
                        para_assi[i].result_len=result_len;
                        para_assi[i].p_hv_tmp=p_hv_tmp;
                        para_assi[i].thread_id=i;
                        para_assi[i].thread_num=thread_num;
                    }

                    uint64_t start=0;
                    for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length;loop_num++)
                    {
                        start=loop_num*thread_num*unit_long;
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            *(para_alloc[j].result_len)=0;//19-01-15
                            para_alloc[j].task_start_pos=start+j*unit_long;
                            if(para_alloc[j].task_start_pos<=seq_length-original_kmer_length-1)
                            {
                                if(start+(j+1)*unit_long>seq_length-original_kmer_length-1)
                                {
                                    para_alloc[j].task_len=seq_length-original_kmer_length-1-para_alloc[j].task_start_pos+1;
                                }
                                else
                                {
                                    para_alloc[j].task_len=unit_long;
                                }
                            }
                        }
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            if(para_alloc[j].task_start_pos<=seq_length-original_kmer_length-1)
                            {
                                if(pthread_create(t+j, NULL, Allocate_HashTable, (void*)(para_alloc+j))!=0)
                                {
                                    cout << "error!" << endl;
                                }
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(para_alloc[i].task_start_pos<=seq_length-original_kmer_length-1)
                            {
                                pthread_join(t[i], NULL);
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(pthread_create(t+i, NULL, Assign_BplusTree_mode, (void*)(para_assi+i))!=0)
                            {
                                cout << "error!" << endl;
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            pthread_join(t[i], NULL);
                        }
                    }
                    free(t);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        free(p_hv_tmp[i]);
                    }
                    free(p_hv_tmp);
                    free(result_len);
                    free(para_alloc);
                    free(para_assi);
                }
                else
                {
                    cal_hash_value_directly_256bit(temp_seq+extend_length,seq_k_current,para1);
                    int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                    if(arrayID!=0)
                    {
                        if(tmp.pos_label==1||tmp.pos_label==3)
                        {
                            uint64_t tmp_pos=extend_length+original_kmer_length-1;
                            uint64_t tmp_pos_divide=tmp_pos/64;
                            uint64_t tmp_pos_mod=tmp_pos%64;
                            tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];
                        }

                        HashAd_first=arrayID-1;
                        HashAd_second=cal_hash_value_directly(temp_seq,extend_length);//ycy20180418am1056
                        uint64_t HashAd_second_4bit=HashAd_second>>1;
                        uint64_t HashAd_second_4bit_label=HashAd_second&1;
                        Transition_4bit(seq[extend_length+original_kmer_length],HashAd[HashAd_second_4bit]+\
                                HashAd_first,HashAd_second_4bit_label);
                    }
                    if(tmp.pos_label==3)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056
                            uint64_t tmp_pos=j+original_kmer_length-1;
                            uint64_t tmp_pos_divide_current=tmp_pos/64;
                            uint64_t tmp_pos_mod_current=tmp_pos%64;
                            uint64_t tmp_template=tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current];
                            if(tmp.template_64[tmp_pos_mod_current]==tmp_template)
                            {
                                arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                                if(arrayID!=0)
                                {
                                    tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]=tmp.p_candidate_pos[ref_i][tmp_pos_divide_current]|tmp.template_64[tmp_pos_mod_current];

                                    HashAd_first=arrayID-1;
                                    HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
                                    uint64_t HashAd_second_4bit=HashAd_second>>1;
                                    uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                    Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
                                            HashAd_first,HashAd_second_4bit_label);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==2)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            uint64_t tmp_pos=j+original_kmer_length-1;
                            uint64_t tmp_pos_divide_current=tmp_pos/64;
                            uint64_t tmp_pos_mod_current=tmp_pos%64;
                            if(tmp.template_64[tmp_pos_mod_current]==(tmp.template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
                            {
                                arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                                if(arrayID!=0)
                                {
                                    HashAd_first=arrayID-1;
                                    HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
                                    uint64_t HashAd_second_4bit=HashAd_second>>1;
                                    uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                    Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
                                            HashAd_first,HashAd_second_4bit_label);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==1)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                            if(arrayID!=0)
                            {
                                uint64_t tmp_pos=j+original_kmer_length-1;
                                uint64_t tmp_pos_divide=tmp_pos/64;
                                uint64_t tmp_pos_mod=tmp_pos%64;
                                tmp.p_candidate_pos[ref_i][tmp_pos_divide]=tmp.p_candidate_pos[ref_i][tmp_pos_divide]|tmp.template_64[tmp_pos_mod];

                                HashAd_first=arrayID-1;
                                HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
                                uint64_t HashAd_second_4bit=HashAd_second>>1;
                                uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
                                        HashAd_first,HashAd_second_4bit_label);
                            }
                        }
                    }
                    else
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-(original_kmer_length);j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);//ycy20180418am1056

                            arrayID = search256BitHashTable(p_root,seq_k_current,para1);
                            if(arrayID!=0)
                            {
                                HashAd_first=arrayID-1;
                                HashAd_second=cal_hash_value_directly(temp_seq+j-extend_length,extend_length);//ycy20180418am1056
                                uint64_t HashAd_second_4bit=HashAd_second>>1;
                                uint64_t HashAd_second_4bit_label=HashAd_second&1;
                                Transition_4bit(seq[j+original_kmer_length],HashAd[HashAd_second_4bit]+\
                                        HashAd_first,HashAd_second_4bit_label);
                            }
                        }
                    }
                }
            }
            if(tmp.pos_label==2||tmp.pos_label==3)
            {
                free(p_candidate_current);
            }
            if(tmp.pos_label==1||tmp.pos_label==3)
            {
                char* posfilename;
                posfilename=(char*)malloc(sizeof(char)*14);
                struct para_getN tmp_inName;
                tmp_inName.kmerlen=original_kmer_length+extend_length;
                tmp_inName.FileName=posfilename;
                tmp_inName.InOut_label=InOut_label;
                tmp_inName.isMiddle=0;
                tmp_inName.isposition=ref_i+1;
                getFileName(tmp_inName);
                posfilename[2]='r';
                posfilename[3]='c';


                FILE* fp_pos;
                fp_pos=fopen(posfilename,"wb");
                fwrite(tmp.p_candidate_pos[ref_i],sizeof(uint64_t),tmp.p_candidate_pos_len[ref_i],fp_pos);
                fclose(fp_pos);
                free(posfilename);
                free(tmp.p_candidate_pos[ref_i]);
                tmp.p_candidate_pos[ref_i]=NULL;
            }
        }
		free(seq);
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout << "fill hash table time is: "<<span<<endl;
		cout <<"write hashtable is over!"<< endl;
	}

	//fill the two-dimensional hash table!

	//check each char space and output the results
	struct timeval tvs_check,tve_check;
	gettimeofday(&tvs_check,NULL);
	int statics_branch_kmer_number=0;
	FILE* out_hash_out;
	FILE* out_hash_out_ad;

	if(InOut_label==0)
	{
		char* outputFileName;
		outputFileName=(char*)malloc(sizeof(char)*11);
		struct para_getN tmp;
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName;
		tmp.InOut_label=0;
		tmp.isMiddle=0;
		tmp.isposition=0;
		getFileName(tmp);
		out_hash_out=fopen(outputFileName,"ab+");

		char* outputFileName_ad;
		struct para_getN tmp_ad;
		if(ad_label==1)
		{
			outputFileName_ad=(char*)malloc(sizeof(char)*13);
			tmp_ad.kmerlen=original_kmer_length+extend_length;
			tmp_ad.FileName=outputFileName_ad;
			tmp_ad.InOut_label=0;
			tmp_ad.isMiddle=2;
			tmp_ad.isposition=0;
			getFileName(tmp_ad);
			out_hash_out_ad=fopen(outputFileName_ad,"ab+");
		}

		if(thread_num>1)
		{
			uint64_t **reslt;
			uint8_t **reslt_c;
			reslt=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
			reslt_c=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);

			uint64_t * reslt_size;
			reslt_size=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

			struct reslt_para_plustree_hashtable* tmp;
			tmp=(struct reslt_para_plustree_hashtable*)malloc(sizeof(struct reslt_para_plustree_hashtable)*thread_num);

			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

			uint64_t unit_size=Tasksize_current/thread_num+1;

			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=j*unit_size;
				if((j+1)*unit_size>Tasksize_current-1)
				{
					tmp[j].end=Tasksize_current-1;
				}
				else
				{
					tmp[j].end=(j+1)*unit_size-1;
				}
				tmp[j].h1=h1;
				tmp[j].HashAd=HashAd;
				tmp[j].reslt=&(reslt[j]);
				tmp[j].reslt_c=&(reslt_c[j]);
				tmp[j].reslt_size=reslt_size+j;
				tmp[j].each_hashTable_size=each_hashTable_size;
				tmp[j].extend_length=extend_length;
				tmp[j].ad_label=ad_label;
				tmp[j].para1=para1;
				tmp[j].para2=para2;
				if(pthread_create(t+j, NULL, judge_thread_in_HashTable, (void*)(tmp+j))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				pthread_join(t[j], NULL);
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				for(uint64_t k=0;k<reslt_size[j];k++)
				{
					fwrite(reslt[j]+k*para2.kmer64Len,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
					if(ad_label==1)
					{
						fwrite(&(reslt_c[j][k]),sizeof(uint8_t),1,out_hash_out_ad);
					}
					statics_branch_kmer_number++;
				}
			}

			free(t);
			free(reslt_size);
			for(uint32_t i=0;i<thread_num;i++)
			{
				free(reslt[i]);
				if(ad_label==1)
				{
					free(reslt_c[i]);
				}
			}
			free(reslt);
			free(reslt_c);
			free(tmp);
		}
		else
		{
			for(uint64_t k=0;k<Tasksize_current;k++)
			{
				uint64_t* h1_temp=h1+k*para1.kmer64Len;
				for(uint64_t j=0;j<each_hashTable_size;j++)
				{
					if(Judge_4bit(HashAd[j][k],0)==1)
					{
						uint64_t newKmer[4];
						cal_newExtendedKmerIn(h1_temp,newKmer,j<<1,extend_length, para1);
						fwrite(newKmer,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
						if(ad_label==1)
						{
							uint8_t c_tmp=(HashAd[j][k]>>4)&15;
							fwrite(&c_tmp,sizeof(uint8_t),1,out_hash_out_ad);
						}
						statics_branch_kmer_number++;
					}
					if(Judge_4bit(HashAd[j][k],1)==1)
					{
						uint64_t newKmer[4];
						cal_newExtendedKmerIn(h1_temp,newKmer,(j<<1)+1,extend_length, para1);
						fwrite(newKmer,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
						if(ad_label==1)
						{
							uint8_t c_tmp=HashAd[j][k]&15;
							fwrite(&c_tmp,sizeof(uint8_t),1,out_hash_out_ad);
						}
						statics_branch_kmer_number++;
					}
				}
			}
		}
		fclose(out_hash_out);
		free(outputFileName);
		if(ad_label==1)
		{
			fclose(out_hash_out_ad);
			free(outputFileName_ad);
		}
	}
	else
	{
		char* outputFileName;
		outputFileName=(char*)malloc(sizeof(char)*11);
		struct para_getN tmp;
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName;
		tmp.InOut_label=1;
		tmp.isMiddle=0;
		tmp.isposition=0;
		getFileName(tmp);
		out_hash_out=fopen(outputFileName,"ab+");

		char* outputFileName_ad;
		struct para_getN tmp_ad;
		if(ad_label==1)
		{
			outputFileName_ad=(char*)malloc(sizeof(char)*13);
			tmp_ad.kmerlen=original_kmer_length+extend_length;
			tmp_ad.FileName=outputFileName_ad;
			tmp_ad.InOut_label=1;
			tmp_ad.isMiddle=2;
			tmp_ad.isposition=0;
			getFileName(tmp_ad);
			out_hash_out_ad=fopen(outputFileName_ad,"ab+");
		}

		if(thread_num>1)
		{
			uint64_t **reslt;
			uint8_t **reslt_c;
			reslt=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
			reslt_c=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);

			uint64_t * reslt_size;
			reslt_size=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

			struct reslt_para_plustree_hashtable* tmp;
			tmp=(struct reslt_para_plustree_hashtable*)malloc(sizeof(struct reslt_para_plustree_hashtable)*thread_num);

			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

			uint64_t unit_size=Tasksize_current/thread_num+1;

			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=j*unit_size;
				if((j+1)*unit_size>Tasksize_current-1)
				{
					tmp[j].end=Tasksize_current-1;
				}
				else
				{
					tmp[j].end=(j+1)*unit_size-1;
				}
				tmp[j].h1=h1;
				tmp[j].HashAd=HashAd;
				tmp[j].reslt=&(reslt[j]);
				tmp[j].reslt_c=&(reslt_c[j]);
				tmp[j].reslt_size=reslt_size+j;
				tmp[j].each_hashTable_size=each_hashTable_size;
				tmp[j].extend_length=extend_length;
				tmp[j].ad_label=ad_label;
				tmp[j].para1=para1;
				tmp[j].para2=para2;

				if(pthread_create(t+j, NULL, judge_thread_out_HashTable, (void*)(tmp+j))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				pthread_join(t[j], NULL);
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				for(uint64_t k=0;k<reslt_size[j];k++)
				{
					fwrite(reslt[j]+k*para2.kmer64Len,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
					if(ad_label==1)
					{
						fwrite(&(reslt_c[j][k]),sizeof(uint8_t),1,out_hash_out_ad);
					}
					statics_branch_kmer_number++;
				}
			}

			free(t);
			free(tmp);
			free(reslt_size);
			for(uint32_t i=0;i<thread_num;i++)
			{
				free(reslt[i]);
				if(ad_label==1)
				{
					free(reslt_c[i]);
				}
			}
			free(reslt);
			free(reslt_c);
		}
		else
		{
			for(uint64_t k=0;k<Tasksize_current;k++)
			{
				uint64_t* h1_temp=h1+k*para1.kmer64Len;
				for(uint64_t j=0;j<each_hashTable_size;j++)
				{
					if(Judge_4bit(HashAd[j][k],0)==1)
					{
						uint64_t newKmer[4];
						cal_newExtendedKmerOut(h1_temp,newKmer,j<<1,extend_length,para1);
						fwrite(newKmer,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
						if(ad_label==1)
						{
							uint8_t c_tmp=(HashAd[j][k]>>4)&15;
							fwrite(&c_tmp,sizeof(uint8_t),1,out_hash_out_ad);
						}
						statics_branch_kmer_number++;
					}
					if(Judge_4bit(HashAd[j][k],1)==1)
					{
						uint64_t newKmer[4];
						cal_newExtendedKmerOut(h1_temp,newKmer,(j<<1)+1,extend_length,para1);
						fwrite(newKmer,sizeof(uint64_t),para2.kmer64Len,out_hash_out);
						if(ad_label==1)
						{
							uint8_t c_tmp=HashAd[j][k]&15;
							fwrite(&c_tmp,sizeof(uint8_t),1,out_hash_out_ad);
						}
						statics_branch_kmer_number++;
					}
				}
			}
		}
		fclose(out_hash_out);
		free(outputFileName);
		if(ad_label==1)
		{
			fclose(out_hash_out_ad);
			free(outputFileName_ad);
		}
	}
	gettimeofday(&tve_check,NULL);
	double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
	cout << "check hash table time is: "<<span_check<<endl;
	cout <<"check each char space and output the results is over!"<< endl;
	cout <<"the # of branched kmer is "<<statics_branch_kmer_number<<endl;
	//check each char space and output the results

	free(h1);
	freeHash256BitTable(p_root);
	for(uint64_t i=0;i<each_hashTable_size;i++)
	{
		free(HashAd[i]);
	}
	free(HashAd);
}
void Cal_Task(struct para_calTask tmp)
{
	char *seq=tmp.seq;
	uint32_t pos_label=tmp.pos_label;
	uint32_t extend_length=tmp.extend_length;
	uint32_t original_kmer_length=tmp.original_kmer_length;
	uint32_t InOut_label=tmp.InOut_label;
	uint32_t isFinal=tmp.isFinal;
	uint32_t ad_label=tmp.ad_label;
	uint64_t memory_size=tmp.memory_size;
	uint32_t thread_num=tmp.thread_num;
	uint64_t Task_Size_Single_Thread=tmp.Task_Size_Single_Thread;
	uint32_t sort_label=tmp.sort_label;
	uint32_t rc_label=tmp.rc_label;

	struct bit256KmerPara para1;
	para1.kmer1Len=original_kmer_length*2;
	para1.remainer1to64=para1.kmer1Len%64;
	para1.kmer64Len=para1.kmer1Len/64+(para1.remainer1to64?1:0);

	uint64_t eachHashTabel=pow(2,2*extend_length-1);
	uint64_t h1_m=sizeof(uint64_t)*(para1.kmer64Len*2+1);
	uint64_t mid_m=3*sizeof(uint64_t)*thread_num*pow(2,Task_Size_Single_Thread);
	uint64_t Hash_head_m=sizeof(uint8_t*)*eachHashTabel;

	uint64_t eachNodeSize=h1_m+eachHashTabel;
	uint64_t Tasksize_PlusTree=(uint64_t)(memory_size*pow(2,30)-mid_m-Hash_head_m)/eachNodeSize;

	cout << sizeof(struct NodeBit) << endl;
	cout << Tasksize_PlusTree << endl;

	char* incAndoutc;
	incAndoutc=(char*)malloc(sizeof(char)*11);
	struct para_getN tmp_inName;
	tmp_inName.kmerlen=original_kmer_length;
	tmp_inName.FileName=incAndoutc;
	tmp_inName.InOut_label=InOut_label;
	tmp_inName.isMiddle=0;
	tmp_inName.isposition=0;
	getFileName(tmp_inName);

	FILE* int_input;
	cout << "the processed kmerfile is :"<< incAndoutc << endl;
	int_input=fopen(incAndoutc,"rb");
	if(int_input==NULL)
	{
		cout << "error: no kmer file find!" << endl;
		return;
	}

	uint64_t total_hash_number;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*para1.kmer64Len);
	fclose(int_input);
	cout << "total number of kmers in processed file is: " <<total_hash_number << endl;

	uint64_t max_F_size=0;

	if (memory_size>=15)
	{
		max_F_size=100000000;
	}
	else if (memory_size>=10)
	{
		max_F_size=50000000;
	}
	else if(memory_size>=5)
	{
		max_F_size=10000000;
	}
	else
	{
		max_F_size=5000000;
	}

	uint64_t Tasksize_PlusTree_F=min(total_hash_number,max_F_size);

	uint64_t **p_candidate_pos;
	uint64_t *p_candidate_pos_len;
	p_candidate_pos=NULL;
	p_candidate_pos_len=NULL;
	struct RefFilePath p_ref_path;
	uint64_t template_64[64];
	template_64[0]=1;
	for(uint32_t i=0;i<63;i++)
	{
		template_64[i+1]=template_64[i]<<1;
	}

	getRefFilePathes(seq, &p_ref_path);
	p_candidate_pos=(uint64_t**)malloc(sizeof(uint64_t*)*p_ref_path.NumberOfPathes);
	p_candidate_pos_len=(uint64_t*)malloc(sizeof(uint64_t)*p_ref_path.NumberOfPathes);

	for(uint32_t i=0;i<p_ref_path.NumberOfPathes;i++)
	{
		p_candidate_pos[i]=NULL;
		p_candidate_pos_len[i]=0;
	}

	uint64_t start_line=0;
	uint64_t end_line=0;
	uint32_t gn=0;

	struct timeval tvs_check_total,tve_check_total;
	gettimeofday(&tvs_check_total,NULL);

	while(end_line<total_hash_number)
	{
		struct timeval tvs_check,tve_check;
		gettimeofday(&tvs_check,NULL);

		start_line=end_line;
		if(isFinal==0)
		{
			struct para_bPlusTreeMap para_call;
			para_call.T_i=gn;
			para_call.seq=seq;
			para_call.start_line=start_line;
			para_call.end_line=&end_line;
			para_call.Tasksize=Tasksize_PlusTree;
			para_call.extend_length=extend_length;
			para_call.original_kmer_length=original_kmer_length;
			para_call.InOut_label=InOut_label;
			para_call.ad_label=ad_label;
			para_call.inputFileName=incAndoutc;
			para_call.thread_num=thread_num;
			para_call.Task_Size_Single_Thread=Task_Size_Single_Thread;
			para_call.pos_label=pos_label;
			para_call.p_candidate_pos=p_candidate_pos;
			para_call.template_64=template_64;
			para_call.p_candidate_pos_len=p_candidate_pos_len;
			para_call.rc_label=rc_label;
			Cal_Single_Task_4bit_UsingHashTable_withLabel_256bitKmer(para_call);

		}
		if(isFinal==1)
		{
			struct para_Count tmp_count;
			tmp_count.T_i=gn;
			tmp_count.seq=seq;
			tmp_count.Tasksize=Tasksize_PlusTree;
			tmp_count.original_kmer_length=original_kmer_length;
			tmp_count.start_line=start_line;
			tmp_count.end_line=&end_line;
			tmp_count.FileName=incAndoutc;
			tmp_count.thread_num=thread_num;
			tmp_count.pos_label=pos_label;
			tmp_count.InOut_label=InOut_label;
			Counting_256bitKmer(tmp_count);
		}
		if(isFinal==2)
		{

			struct para_Final tmp_final;
			tmp_final.T_i=gn;
			tmp_final.seq=seq;
			tmp_final.Tasksize=Tasksize_PlusTree_F;
			tmp_final.original_kmer_length=original_kmer_length;
			tmp_final.extend_length=extend_length;
			tmp_final.start_line=start_line;
			tmp_final.end_line=&end_line;
			tmp_final.InOut_label=InOut_label;
			tmp_final.inputFileName=incAndoutc;
			tmp_final.ad_label=ad_label;
			tmp_final.thread_num=thread_num;
			tmp_final.Task_Size_Single_Thread=Task_Size_Single_Thread;
			tmp_final.p_candidate_pos=p_candidate_pos;
			tmp_final.template_64=template_64;
			tmp_final.pos_label=pos_label;
			tmp_final.p_candidate_pos_len=p_candidate_pos_len;
			tmp_final.rc_label=rc_label;
			Cal_Single_Task_4bit_UsingHashTable_withLabel_Final_256bitKmer(tmp_final);
		}
		gn++;

		gettimeofday(&tve_check,NULL);
		double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
		cout << start_line << "," <<end_line <<","<<total_hash_number<<endl;
		cout<<(((double)end_line)/((double)total_hash_number))<< endl;
		cout <<"the" << gn-1 <<"th group time is: " << span_check<< endl;
	}
	if(isFinal==0&&InOut_label==1&&sort_label==1)
	{
		char* outputFileName;
		outputFileName=(char*)malloc(sizeof(char)*11);
		struct para_getN tmp;
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName;
		tmp.InOut_label=1;
		tmp.isMiddle=0;
		tmp.isposition=0;
		getFileName(tmp);

		char* outputFileName_ad;
		outputFileName_ad=(char*)malloc(sizeof(char)*13);
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName_ad;
		tmp.InOut_label=InOut_label;
		tmp.isMiddle=2;
		tmp.isposition=0;
		getFileName(tmp);

		char tmp_name[4]="tmp";
		char tmp_name_ad[7]="tmp_ad";

		uint32_t len=(original_kmer_length+extend_length)/32+(((original_kmer_length+extend_length)%32)?1:0);
		if(ad_label==0)
		{
			SortKmers(outputFileName,tmp_name,NULL,NULL,len);
			remove(outputFileName);
			rename(tmp_name,outputFileName);
		}
		else
		{
			SortKmers(outputFileName,tmp_name,outputFileName_ad,tmp_name_ad,len);
			remove(outputFileName);
			remove(outputFileName_ad);
			rename(tmp_name,outputFileName);
			rename(tmp_name_ad,outputFileName_ad);
		}

		free(outputFileName);
		free(outputFileName_ad);
	}
	if(isFinal==2&&sort_label==1)
	{
		char* outputFileName;
		outputFileName=(char*)malloc(sizeof(char)*11);
		struct para_getN tmp;
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName;
		tmp.InOut_label=InOut_label;
		tmp.isMiddle=0;
		tmp.isposition=0;
		getFileName(tmp);

		char* outputFileName_ad;
		outputFileName_ad=(char*)malloc(sizeof(char)*13);
		tmp.kmerlen=original_kmer_length+extend_length;
		tmp.FileName=outputFileName_ad;
		tmp.InOut_label=InOut_label;
		tmp.isMiddle=2;
		tmp.isposition=0;
		getFileName(tmp);

		char tmp_name[4]="tmp";
		char tmp_name_ad[7]="tmp_ad";
		uint32_t len=(original_kmer_length+extend_length)/32+(((original_kmer_length+extend_length)%32)?1:0);
		if(ad_label==0)
		{
			SortKmers(outputFileName,tmp_name,NULL,NULL,len);
			remove(outputFileName);
			rename(tmp_name,outputFileName);
		}
		else
		{
			SortKmers(outputFileName,tmp_name,outputFileName_ad,tmp_name_ad,len);
			remove(outputFileName);
			remove(outputFileName_ad);
			rename(tmp_name,outputFileName);
			rename(tmp_name_ad,outputFileName_ad);
		}

		free(outputFileName);
		free(outputFileName_ad);
	}
	if(isFinal==1&&sort_label==1)
	{
		char name[]="tmp";
		char name_tmp[]="tmps";
		uint32_t len_total;
		len_total=(original_kmer_length)/32+(((original_kmer_length)%32)?1:0);
		SortOutCountKmersUsingBplusTreebinarytobinary(name,name_tmp,len_total);
		remove("tmp");
	}
	free(p_candidate_pos);
	gettimeofday(&tve_check_total,NULL);
	double span_check_total = tve_check_total.tv_sec-tvs_check_total.tv_sec + (tve_check_total.tv_usec-tvs_check_total.tv_usec)/1000000.0;
	cout <<"the total time for current extending is: " << span_check_total<< endl;
}
