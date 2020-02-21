#include "Final_Extend.h"

/***************************************Final for 256bitKmer******************/
void *Allocate_BplusTree_count_256bitKmer(void *arg)
{
	struct Parallel_para_count01_256bitKmer *tmp=(struct Parallel_para_count01_256bitKmer *)arg;
	uint64_t task_start_pos=tmp->task_start_pos;
	char *seq=tmp->seq;
	uint64_t task_len=tmp->task_len;

	uint32_t kmer_len1=tmp->kmer_len1;
	uint32_t InOut_label=tmp->InOut_label;
	uint64_t *p_candidate_current=tmp->p_candidate_current;
	uint64_t *template_64=tmp->template_64;
	uint32_t pos_label=tmp->pos_label;

	struct bit256KmerPara para=tmp->para;

	struct TwoDimenHashV_count_256bitKmer* p_hv_tmp=tmp->p_hv_tmp;
	uint64_t p_hv_tmp_length_cur=0;

	uint64_t seq_k_current[4];
	char*temp_seq=seq+task_start_pos;


	cal_hash_value_directly_256bit(temp_seq,seq_k_current,para);
	uint64_t tmp_pos_divide_current;
	uint64_t tmp_pos_mod_current;
	if(pos_label==2||pos_label==3)
	{
		if(InOut_label==0)
		{
			tmp_pos_divide_current=task_start_pos/64;
			tmp_pos_mod_current=task_start_pos%64;
		}
		else
		{
			tmp_pos_divide_current=(task_start_pos+kmer_len1-1)/64;
			tmp_pos_mod_current=(task_start_pos+kmer_len1-1)%64;
		}
		if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
		{
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_hv_tmp[p_hv_tmp_length_cur].hv[j]=seq_k_current[j];
			}
			p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos;
			p_hv_tmp_length_cur++;
		}
		for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
		{
			cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);
			if(InOut_label==0)
			{
				tmp_pos_divide_current=(task_start_pos+j)/64;
				tmp_pos_mod_current=(task_start_pos+j)%64;
			}
			else
			{
				tmp_pos_divide_current=(task_start_pos+j+kmer_len1-1)/64;
				tmp_pos_mod_current=(task_start_pos+j+kmer_len1-1)%64;
			}
			if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
			{
				for(uint32_t k=0;k<para.kmer64Len;k++)
				{
					p_hv_tmp[p_hv_tmp_length_cur].hv[k]=seq_k_current[k];
				}
				p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
				p_hv_tmp_length_cur++;
			}
		}
	}
	else
	{
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			p_hv_tmp[p_hv_tmp_length_cur].hv[j]=seq_k_current[j];
		}
		p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos;
		p_hv_tmp_length_cur++;
		for(uint64_t j=1;j<task_len;j++)//check each pos without the first in the case of "branch_in"
		{
			cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para);
			for(uint32_t k=0;k<para.kmer64Len;k++)
			{
				p_hv_tmp[p_hv_tmp_length_cur].hv[k]=seq_k_current[k];
			}
			p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
			p_hv_tmp_length_cur++;
		}
	}

	*(tmp->result_len)=p_hv_tmp_length_cur;
}
void *Assign_BplusTree_count_256bitKmer(void *arg)
{
	struct Parallel_para_count02_256bitKmer *tmp=(struct Parallel_para_count02_256bitKmer *)arg;

	struct NodeBit ** p_root=tmp->p_root;
	uint64_t *result_len=tmp->result_len;
	struct TwoDimenHashV_count_256bitKmer** p_hv_tmp=tmp->p_hv_tmp;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	struct bit256KmerPara para=tmp->para;

	for(uint32_t i=0;i<thread_num;i++)
	{
		for(uint64_t k=0;k<result_len[i];k++)
		{
			if(p_hv_tmp[i][k].hv[0]%thread_num==thread_id)
			{
				count_hashValue_bit(p_root[bit256hashFFunction(p_hv_tmp[i][k].hv,para)],p_hv_tmp[i][k].hv,para);
			}
		}
	}

	return NULL;
}
void *Allocate_HashTable_Final_256bitKmer(void *arg)
{
	struct Parallel_para_Final01_hash_256bitKmer *tmp=(struct Parallel_para_Final01_hash_256bitKmer *)arg;
	uint64_t task_start_pos=tmp->task_start_pos;
	char *seq=tmp->seq;
	struct bit256Hash * p_root=tmp->p_root;
	uint32_t kmer_len1=tmp->kmer_len1;
	uint32_t kmer_len2=tmp->kmer_len2;
	uint64_t task_len=tmp->task_len;
	uint32_t InOut_label=tmp->InOut_label;
	struct bit256KmerPara para1=tmp->para1;
	struct bit256KmerPara para2=tmp->para2;
	uint64_t seq_length=tmp->seq_length;

	uint32_t pos_label=tmp->pos_label;
	uint64_t *p_candidate_pos=tmp->p_candidate_pos;
	uint64_t *template_64=tmp->template_64;
	uint64_t *p_candidate_current=tmp->p_candidate_current;

	struct TwoDimenHashV_Final_256bitKmer* p_hv_tmp=tmp->p_hv_tmp;
	uint64_t p_hv_tmp_length_cur=0;

	uint64_t seq_k_current[4];
	uint64_t seq_k_currentF[4];
	uint32_t original_kmer_length=kmer_len1;
	uint32_t extend_length=kmer_len2;
	char*temp_seq=seq+task_start_pos;

	if(InOut_label==0)
	{
		cal_hash_value_directly_256bit(temp_seq,seq_k_current,para1);
		if(task_start_pos!=0)
		{
			if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
			{
				if(pos_label==1||pos_label==3)
				{
					uint64_t tmp_pos_divide=task_start_pos/64;
					uint64_t tmp_pos_mod=task_start_pos%64;
					p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];
				}

				cal_hash_value_directly_256bit(temp_seq,seq_k_currentF,para2);
				for(uint32_t j=0;j<para2.kmer64Len;j++)
				{
					p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[j]=seq_k_currentF[j];
				}
				p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos;
				uint8_t ad_temp=0;
				Transition(seq[task_start_pos-1],'0',&ad_temp);
				p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
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
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				uint64_t tmp_pos_divide=(task_start_pos+j)/64;
				uint64_t tmp_pos_mod=(task_start_pos+j)%64;
				uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
				if(template_64[tmp_pos_mod]==tmp_template)
				{
					if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
					{
						p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];

						cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
						for(uint32_t k=0;k<para2.kmer64Len;k++)
						{
							p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
						}
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						uint8_t ad_temp=0;
						Transition(temp_seq[j-1],'0',&ad_temp);
						p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
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

				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				uint64_t tmp_pos_divide=(task_start_pos+j)/64;
				uint64_t tmp_pos_mod=(task_start_pos+j)%64;
				uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
				if(template_64[tmp_pos_mod]==tmp_template)
				{
					if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
					{
						cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
						for(uint32_t k=0;k<para2.kmer64Len;k++)
						{
							p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
						}
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						uint8_t ad_temp=0;
						Transition(temp_seq[j-1],'0',&ad_temp);
						p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
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
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
				{
					uint64_t tmp_pos_divide=(task_start_pos+j)/64;
					uint64_t tmp_pos_mod=(task_start_pos+j)%64;
					p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];

					cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
					for(uint32_t k=0;k<para2.kmer64Len;k++)
					{
						p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
					}
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					uint8_t ad_temp=0;
					Transition(temp_seq[j-1],'0',&ad_temp);
					p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
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
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);

				if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
				{
					cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
					for(uint32_t k=0;k<para2.kmer64Len;k++)
					{
						p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
					}
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					uint8_t ad_temp=0;
					Transition(temp_seq[j-1],'0',&ad_temp);
					p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
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
		cal_hash_value_directly_256bit(temp_seq+task_start_pos_out-task_start_pos,seq_k_current,para1);
		if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
		{
			if(tmp->pos_label==1||tmp->pos_label==3)
			{
				uint64_t tmp_pos=task_start_pos_out+original_kmer_length-1;
				uint64_t tmp_pos_divide=tmp_pos/64;
				uint64_t tmp_pos_mod=tmp_pos%64;
				p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];
			}

			cal_hash_value_directly_256bit(temp_seq+task_start_pos_out-task_start_pos-extend_length,seq_k_currentF,para2);
			for(uint32_t j=0;j<para2.kmer64Len;j++)
			{
				p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[j]=seq_k_currentF[j];
			}
			p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos_out;
			uint8_t ad_temp=0;
			Transition('0',seq[task_start_pos+original_kmer_length],&ad_temp);
			p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
			p_hv_tmp_length_cur++;
		}
		if(pos_label==3)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)
			{
				if(j+task_start_pos>seq_length-original_kmer_length-1)
				{
					break;
				}
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);

				uint64_t tmp_pos=task_start_pos+original_kmer_length+j-1;
				uint64_t tmp_pos_divide=tmp_pos/64;
				uint64_t tmp_pos_mod=tmp_pos%64;
				uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
				if(template_64[tmp_pos_mod]==tmp_template)
				{
					if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
					{
						p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];

						cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
						for(uint32_t k=0;k<para2.kmer64Len;k++)
						{
							p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
						}
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						uint8_t ad_temp=0;
						Transition('0',temp_seq[original_kmer_length+j],&ad_temp);
						p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
						p_hv_tmp_length_cur++;
					}
				}
			}
		}
		else if(pos_label==2)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)
			{
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				uint64_t tmp_pos=task_start_pos+original_kmer_length+j-1;
				uint64_t tmp_pos_divide=tmp_pos/64;
				uint64_t tmp_pos_mod=tmp_pos%64;
				uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
				if(template_64[tmp_pos_mod]==tmp_template)
				{
					if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
					{
						cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
						for(uint32_t k=0;k<para2.kmer64Len;k++)
						{
							p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
						}
						p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
						uint8_t ad_temp=0;
						Transition('0',temp_seq[original_kmer_length+j],&ad_temp);
						p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
						p_hv_tmp_length_cur++;
					}
				}
			}
		}
		else if(pos_label==1)
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)
			{
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
				{
					uint64_t tmp_pos=task_start_pos+original_kmer_length+j-1;
					uint64_t tmp_pos_divide=tmp_pos/64;
					uint64_t tmp_pos_mod=tmp_pos%64;
					p_candidate_pos[tmp_pos_divide]=p_candidate_pos[tmp_pos_divide]|template_64[tmp_pos_mod];

					cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
					for(uint32_t k=0;k<para2.kmer64Len;k++)
					{
						p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
					}
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					uint8_t ad_temp=0;
					Transition('0',temp_seq[original_kmer_length+j],&ad_temp);
					p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
					p_hv_tmp_length_cur++;
				}
			}
		}
		else
		{
			for(uint64_t j=task_start_pos_out-task_start_pos+1;j<task_len;j++)
			{
				cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
				if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
				{
					cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
					for(uint32_t k=0;k<para2.kmer64Len;k++)
					{
						p_hv_tmp[p_hv_tmp_length_cur].tmp.hashValue[k]=seq_k_currentF[k];
					}
					p_hv_tmp[p_hv_tmp_length_cur].pos=task_start_pos+j;
					uint8_t ad_temp=0;
					Transition('0',temp_seq[original_kmer_length+j],&ad_temp);
					p_hv_tmp[p_hv_tmp_length_cur].tmp.ad=ad_temp;
					p_hv_tmp_length_cur++;
				}
			}
		}
	}
	*(tmp->result_len)=p_hv_tmp_length_cur;
	return NULL;
}
void *Assign_HashTable_Final_256bitKmer(void *arg)
{
	struct Parallel_para_Final02_hash_256bitKmer *tmp=(struct Parallel_para_Final02_hash_256bitKmer *)arg;
	struct NodeBit ** p_root_hashtable=tmp->p_root_hashtable;
	uint64_t *result_len=tmp->result_len;
	struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp=tmp->p_hv_tmp;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	struct bit256KmerPara para=tmp->para;

	for(uint32_t i=0;i<thread_num;i++)
	{
		for(uint64_t k=0;k<result_len[i];k++)
		{
			if(bit256hashFFunction(p_hv_tmp[i][k].tmp.hashValue,para)%thread_num==thread_id)
			{
				bit256insertHashFTable(p_root_hashtable,p_hv_tmp[i][k].tmp,para);
			}
		}
	}
	return NULL;
}
void *judge_thread_plustree_final_256bitKmer(void *arg)
{
	struct reslt_para_plustree_final_256bitKmer* tmp=(struct reslt_para_plustree_final_256bitKmer*)arg;

	struct NodeBit * start=tmp->p_root_start;
	struct NodeBit * end=tmp->p_root_end;
	uint32_t InOut_label=tmp->InOut_label;
	uint32_t ad_label=tmp->ad_label;
	struct bit256KmerPara para=tmp->para;
	uint64_t reslt_Max_size;
	uint64_t *reslt;
	uint32_t *reslt_c;

	uint64_t reslt_size=0;
	reslt=(uint64_t *)malloc(sizeof(uint64_t)*pow(2,10)*para.kmer64Len);
	if(ad_label==1)
	{
		reslt_c=(uint32_t *)malloc(sizeof(uint32_t)*pow(2,10));
	}
	reslt_Max_size=pow(2,10);
	struct NodeBit *while_tmp=start;
	while(while_tmp!=end)
	{
		for(uint32_t i=0;i<while_tmp->Node_Size;i++)
		{
			if(Judge_4bit(while_tmp->data[i].ad,InOut_label))
			{
				if(reslt_size<reslt_Max_size)
				{
					for(uint32_t j=0;j<para.kmer64Len;j++)
					{
						reslt[reslt_size*para.kmer64Len+j]=while_tmp->data[i].hashValue[j];
					}
					if(ad_label==1)
					{
						uint8_t c_tmp;
						if(InOut_label==0)
						{
							c_tmp=(while_tmp->data[i].ad>>4)&15;
						}
						else
						{
							c_tmp=(while_tmp->data[i].ad)&15;
						}
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
				else
				{
					reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*para.kmer64Len*(reslt_Max_size+pow(2,10)));
					if(ad_label==1)
					{
						reslt_c=(uint32_t*)realloc(reslt_c,sizeof(uint32_t)*(reslt_Max_size+pow(2,10)));
					}
					reslt_Max_size=reslt_Max_size+pow(2,10);

					for(uint32_t j=0;j<para.kmer64Len;j++)
					{
						reslt[reslt_size*para.kmer64Len+j]=while_tmp->data[i].hashValue[j];
					}
					if(ad_label==1)
					{
						uint8_t c_tmp;
						if(InOut_label==0)
						{
							c_tmp=(while_tmp->data[i].ad>>4)&15;
						}
						else
						{
							c_tmp=(while_tmp->data[i].ad)&15;
						}
						reslt_c[reslt_size]=c_tmp;
					}
					reslt_size++;
				}
			}
		}
		while_tmp=while_tmp->brother;
	}

	*(tmp->reslt_size)=reslt_size;
	*(tmp->reslt)=reslt;
	*(tmp->reslt_c)=reslt_c;
	return NULL;
}

void Insert_ad_256bitKmer(struct NodeBit** p_root,struct nodeBit insert_data,struct bit256KmerPara para)
{
	struct NodeBit*p_root_tmp= * p_root;
	if(p_root_tmp==NULL)
	{
		cout << "heh" << endl;
	}
		//insert lead node
	while(1)
	{
		if(p_root_tmp->leaf_label==1)
		{
			break;
		}
		uint32_t child_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				uint32_t gi;
				if(i==0)
				{
					gi=0;
				}
				else
				{
					gi=i-1;
				}
				p_root_tmp=p_root_tmp->p_child[gi];
				child_label=1;
				break;
			}
		}
		if(child_label==0)
		{
			p_root_tmp=p_root_tmp->p_child[p_root_tmp->Node_Size-1];
		}
	}

	uint32_t pos;
	if(p_root_tmp->Node_Size==0)
	{
		pos=0;
	}
	else
	{
		int32_t pos_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==2)
			{
				pos=i;
				pos_label=-1;
				break;
			}
			else if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				pos=i;
				pos_label=1;
				break;
			}
		}
		if(pos_label==-1)
		{
			p_root_tmp->data[pos].ad=p_root_tmp->data[pos].ad|insert_data.ad;
			return;
		}
		else if(pos_label==0)
		{
			pos=p_root_tmp->Node_Size;
		}
	}

	if(p_root_tmp->Node_Size==M)
	{
		Divide_Node_bit(p_root,p_root_tmp,insert_data,pos,NULL,1,0,para);
	}
	else
	{
		for(uint32_t i=pos;i<p_root_tmp->Node_Size;i++)
		{
			uint32_t l=p_root_tmp->Node_Size-1-(i-pos);
			p_root_tmp->data[l+1].arrayID=p_root_tmp->data[l].arrayID;
			if(p_root_tmp->data[l+1].hashValue==NULL)
			{
				p_root_tmp->data[l+1].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_root_tmp->data[l+1].hashValue[j]=p_root_tmp->data[l].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_root_tmp->data[l+1].hashValue[j]=p_root_tmp->data[l].hashValue[j];
				}
			}
			p_root_tmp->data[l+1].ad=p_root_tmp->data[l].ad;
		}
		p_root_tmp->data[pos].arrayID=insert_data.arrayID;
		if(p_root_tmp->data[pos].hashValue==NULL)
		{
			p_root_tmp->data[pos].hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_root_tmp->data[pos].hashValue[j]=insert_data.hashValue[j];
			}
		}
		else
		{
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_root_tmp->data[pos].hashValue[j]=insert_data.hashValue[j];
			}
		}
		p_root_tmp->data[pos].ad=insert_data.ad;

		p_root_tmp->Node_Size++;

		if(pos==0&&p_root_tmp->parent!=NULL)
		{
			NodeBit* p_parent_update=p_root_tmp->parent;
			p_parent_update->data[p_root_tmp->numAsChild].arrayID=\
					p_root_tmp->data[0].arrayID;
			if(p_parent_update->data[p_root_tmp->numAsChild].hashValue==NULL)
			{
				p_parent_update->data[p_root_tmp->numAsChild].hashValue=\
						(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_parent_update->data[p_root_tmp->numAsChild].hashValue[j]=\
							p_root_tmp->data[0].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_parent_update->data[p_root_tmp->numAsChild].hashValue[j]=\
							p_root_tmp->data[0].hashValue[j];
				}
			}
			NodeBit* p1=p_root_tmp;
			NodeBit* p2=p_parent_update;
			while(p1->numAsChild==0&&p2->parent!=NULL)
			{
				p2->parent->data[p2->numAsChild].arrayID=\
						p2->data[0].arrayID;
				if(p2->parent->data[p2->numAsChild].hashValue==NULL)
				{
					p2->parent->data[p2->numAsChild].hashValue=\
							(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
					for(uint32_t j=0;j<para.kmer64Len;j++)
					{
						p2->parent->data[p2->numAsChild].hashValue[j]=\
								p2->data[0].hashValue[j];
					}
				}
				else
				{
					for(uint32_t j=0;j<para.kmer64Len;j++)
					{
						p2->parent->data[p2->numAsChild].hashValue[j]=\
								p2->data[0].hashValue[j];
					}
				}
				p1=p2;
				p2=p2->parent;
			}
		}

	}
}
void Counting_256bitKmer(struct para_Count tmp)
{
	uint32_t T_i=tmp.T_i;
	char *seq;
	char *dataset=tmp.seq;
	uint64_t seq_length;
	uint64_t Tasksize=tmp.Tasksize;
	uint32_t original_kmer_length=tmp.original_kmer_length;
	uint64_t start_line=tmp.start_line;
	uint64_t *end_line=tmp.end_line;
	char* FileName=tmp.FileName;
	uint32_t thread_num=tmp.thread_num;
	uint32_t pos_label=tmp.pos_label;
	uint32_t InOut_label=tmp.InOut_label;

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

	char *temp_seq;

	/**********************************read kmers***********************************/
	FILE * f_in;
	f_in=fopen(FileName,"rb");
	uint64_t tasksize_tmp=0;
	fseek(f_in,0,2);
	tasksize_tmp=ftell(f_in)/(sizeof(uint64_t)*para1.kmer64Len);
	if(tasksize_tmp-start_line<Tasksize)
	{
		Tasksize=tasksize_tmp-start_line;
	}
	fseek(f_in,0,0);
	uint64_t *h1;
	h1=(uint64_t*)malloc(sizeof(uint64_t)*Tasksize*para1.kmer64Len);
	fseek(f_in,(start_line)*sizeof(uint64_t)*para1.kmer64Len,0);
	uint32_t Tasksize_current=fread(h1,sizeof(uint64_t),para1.kmer64Len*Tasksize,f_in);
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

	//initial a Bplus tree for saving the hashvalues
	struct NodeBit ** p_root;
	p_root=bit256initialHashFTable();
	for(uint64_t i=0;i<Tasksize_current;i++)
	{
		struct nodeBit tmp_plus;
		tmp_plus.hashValue=h1+i*para1.kmer64Len;
		tmp_plus.arrayID=0;
		uint32_t x=bit256hashFFunction(tmp_plus.hashValue,para1);
		if(p_root[x]==NULL)
		{
			p_root[x]=bit256initialSingleHashFTable();
		}
		Insert_Value_bit(p_root+x,tmp_plus,para1);
	}
	cout << "Constructing BplusTree is over!" <<endl;
	//initial a Bplus tree for saving the hashvalues

	//count the hash value!
	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		temp_seq=seq;

		uint64_t *p_candidate_current=NULL;
		if(pos_label==2||pos_label==3)
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

		uint64_t template_64[64];
		template_64[0]=1;
		for(uint32_t i=0;i<63;i++)
		{
			template_64[i+1]=template_64[i]<<1;
		}

		if(thread_num>1)
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			uint64_t unit_long=pow(2,25);
			struct TwoDimenHashV_count_256bitKmer** p_hv_tmp;
			p_hv_tmp=(struct TwoDimenHashV_count_256bitKmer**)malloc(sizeof(struct TwoDimenHashV_count_256bitKmer*)*thread_num);
			for(uint32_t i=0;i<thread_num;i++)
			{
				p_hv_tmp[i]=(struct TwoDimenHashV_count_256bitKmer*)malloc(sizeof(struct TwoDimenHashV_count_256bitKmer)*unit_long);
				for(uint32_t j=0;j<unit_long;j++)
				{
					p_hv_tmp[i][j].hv=(uint64_t*)malloc(sizeof(uint64_t)*para1.kmer64Len);
				}
			}
			uint64_t* result_len;
			result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
			struct Parallel_para_count01_256bitKmer* para_alloc;
			para_alloc=(struct Parallel_para_count01_256bitKmer*)malloc(sizeof(struct Parallel_para_count01_256bitKmer)*thread_num);
			for(uint32_t i=0;i<thread_num;i++)
			{
				para_alloc[i].seq=seq;
				para_alloc[i].seq_length=seq_length;
				para_alloc[i].kmer_len1=original_kmer_length;
				para_alloc[i].task_len=unit_long;
				para_alloc[i].result_len=result_len+i;
				para_alloc[i].p_hv_tmp=p_hv_tmp[i];
				para_alloc[i].para=para1;
				para_alloc[i].p_candidate_current=p_candidate_current;
				para_alloc[i].template_64=template_64;
				para_alloc[i].pos_label=pos_label;
				para_alloc[i].InOut_label=InOut_label;
			}

			struct Parallel_para_count02_256bitKmer* para_assi;
			para_assi=(struct Parallel_para_count02_256bitKmer*)malloc(sizeof(struct Parallel_para_count02_256bitKmer)*thread_num);
			for(uint32_t i=0;i<thread_num;i++)
			{
				para_assi[i].seq=seq;
				para_assi[i].kmer_len1=original_kmer_length;
				para_assi[i].result_len=result_len;
				para_assi[i].p_hv_tmp=p_hv_tmp;
				para_assi[i].thread_id=i;
				para_assi[i].thread_num=thread_num;
				para_assi[i].p_root=p_root;
				para_assi[i].para=para1;
			}

			uint64_t start=0;
			for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length;loop_num++)
			{
				start=loop_num*(thread_num*unit_long);
				for(uint32_t j=0;j<thread_num;j++)
				{
					para_alloc[j].task_start_pos=start+j*unit_long;
					if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length))
					{
						if(start+(j+1)*unit_long>seq_length-(original_kmer_length))
						{
							para_alloc[j].task_len=seq_length-(original_kmer_length)-para_alloc[j].task_start_pos+1;
						}
						else
						{
							para_alloc[j].task_len=unit_long;
						}
					}
				}
				for(uint32_t j=0;j<thread_num;j++)
				{
					result_len[j]=0;
					if(para_alloc[j].task_start_pos<=seq_length-(original_kmer_length))
					{
						if(pthread_create(t+j, NULL, Allocate_BplusTree_count_256bitKmer, (void*)(para_alloc+j))!=0)
						{
							cout << "error!" << endl;
						}
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length))
					{
						pthread_join(t[i], NULL);
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
//					if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length))
//					{
						if(pthread_create(t+i, NULL, Assign_BplusTree_count_256bitKmer, (void*)(para_assi+i))!=0)
						{
							cout << "error!" << endl;
						}
//					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
//					if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length))
//					{
						pthread_join(t[i], NULL);
//					}
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
			uint64_t seq_k_current[4];
			uint64_t tmp_pos_divide_current;
			uint64_t tmp_pos_mod_current;
			cal_hash_value_directly_256bit(temp_seq,seq_k_current,para1);
			if(pos_label==2||pos_label==3)
			{
				if(InOut_label==0)
				{
					tmp_pos_divide_current=0;
					tmp_pos_mod_current=0;
				}
				else
				{
					tmp_pos_divide_current=(original_kmer_length-1)/64;
					tmp_pos_mod_current=(original_kmer_length-1)%64;
				}
				if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
				{
					count_hashValue_bit(p_root[bit256hashFFunction(seq_k_current,para1)],seq_k_current,para1);
				}
				for(uint64_t j=1;j<seq_length-(original_kmer_length)+1;j++)//check each pos without the first in the case of "branch_in"
				{
					cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
					if(InOut_label==0)
					{
						tmp_pos_divide_current=j/64;
						tmp_pos_mod_current=j%64;
					}
					else
					{
						tmp_pos_divide_current=(j+original_kmer_length-1)/64;
						tmp_pos_mod_current=(j+original_kmer_length-1)%64;
					}
					if(template_64[tmp_pos_mod_current]==(template_64[tmp_pos_mod_current]&p_candidate_current[tmp_pos_divide_current]))
					{
						count_hashValue_bit(p_root[bit256hashFFunction(seq_k_current,para1)],seq_k_current,para1);
					}
				}
			}
			else
			{
				count_hashValue_bit(p_root[bit256hashFFunction(seq_k_current,para1)],seq_k_current,para1);
				for(uint64_t j=1;j<seq_length-(original_kmer_length)+1;j++)//check each pos without the first in the case of "branch_in"
				{
					cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
					count_hashValue_bit(p_root[bit256hashFFunction(seq_k_current,para1)],seq_k_current,para1);
				}
			}
		}
		free(seq);
	}
	cout <<"write hashValue count is over!"<< endl;
	//count the hash value!

	//output the results
	FILE *out_hash;
	char name_tmp[]="tmp";
	out_hash=fopen(name_tmp,"ab+");
	struct NodeBit* p_minimal;
	p_minimal=Find_Minimal_Node_bit(p_root[0]);
	uint32_t imin;
	for(imin=0;imin<HashFSize-1;imin++)
	{
		if(p_root[imin]!=NULL)
		{
			p_minimal=Find_Minimal_Node_bit(p_root[imin]);
			break;
		}
	}
	for(uint32_t i=imin;i<HashFSize-1;)
	{
		struct NodeBit *p_leaf_hashtable_tmp=Find_Minimal_Node_bit(p_root[i]);
		while(p_leaf_hashtable_tmp->brother!=NULL)
		{
			p_leaf_hashtable_tmp=p_leaf_hashtable_tmp->brother;
		}
		uint32_t isave=i;
		for(uint32_t j=i+1;j<HashFSize;j++)
		{
			if(p_root[j]!=NULL)
			{
				p_leaf_hashtable_tmp->brother=Find_Minimal_Node_bit(p_root[j]);
				i=j;
				break;
			}
		}
		if(i==isave)
		{
			break;
		}
	}
	uint32_t y=0;
	while(p_minimal!=NULL)
	{
		y=y+p_minimal->Node_Size;
		for(uint32_t i=0;i<p_minimal->Node_Size;i++)
		{
			fwrite(p_minimal->data[i].hashValue,sizeof(uint64_t),para1.kmer64Len,out_hash);
			fwrite(&(p_minimal->data[i].arrayID),sizeof(uint64_t),1,out_hash);
			if(p_minimal->data[i].arrayID<2)
			{
				cout << "error: from kmer counter, count is small than 2!" << endl;
				return ;
			}
		}
		p_minimal=p_minimal->brother;
	}
	fclose(out_hash);
	cout <<"output the results is over!"<< endl;
	free(h1);
	bit256freeHashFTable(p_root,para1);

	//output the results
//	SortOutCountKmersSingle("tmp_tmp","tmp",para1.kmer64Len);

}
void Cal_Single_Task_4bit_UsingHashTable_withLabel_Final_256bitKmer(struct para_Final tmp)
{
	uint32_t T_i=tmp.T_i;
	char *seq;
	char * dataset=tmp.seq;
	uint64_t seq_length;
	uint64_t Tasksize=tmp.Tasksize;
	uint32_t original_kmer_length=tmp.original_kmer_length;
	uint32_t extend_length=tmp.extend_length;
	uint64_t start_line=tmp.start_line;
	uint64_t *end_line=tmp.end_line;
	uint32_t InOut_label=tmp.InOut_label;
	char*inputFileName=tmp.inputFileName;
	uint32_t ad_label=tmp.ad_label;
	uint32_t thread_num=tmp.thread_num;
	uint64_t **p_candidate_pos=tmp.p_candidate_pos;
	uint64_t * p_candidate_pos_len=tmp.p_candidate_pos_len;
	uint64_t * template_64=tmp.template_64;
	uint64_t Task_Size_Single_Thread=tmp.Task_Size_Single_Thread;
	uint32_t rc_label=tmp.rc_label;
	Task_Size_Single_Thread=pow(2,Task_Size_Single_Thread);
	cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	cerr << "start the " << T_i << "th group" << endl;

	/**********************cal para for original and final kmers*****************/
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
	/**********************cal para for original and final kmers*****************/

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

	//hash the first hash value to the first dimensional coordinate of HashAd
	struct bit256Hash* p_root=NULL;
	p_root=initial256BitHashTable();
	if(p_root==NULL)
	{
		cerr<< "error!" << endl;
	}
	struct timeval tvs_hs,tve_hs;
	gettimeofday(&tvs_hs,NULL);
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

	free(h1);

	/***************************define the hashtable tree***************************/
	struct NodeBit ** p_root_hashtable;
	p_root_hashtable=bit256initialHashFTable();
	/***************************define the hashtable tree***************************/

	struct nodeBit c_tmp_hashtable;
	c_tmp_hashtable.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para2.kmer64Len);

	/**********************fill the two-dimensional hash table!*********************/
	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	char *temp_seq;
	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
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
			uint64_t x;
			x=fread(p_candidate_current,sizeof(uint64_t),seq_length/64+((seq_length%64)?1:0),fp_pos);
			if(x!=seq_length/64+((seq_length%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos);
		}

		struct timeval tvs,tve;
		gettimeofday(&tvs,NULL);
		uint64_t seq_k_current[4];
		uint64_t seq_k_currentF[4];
		if(InOut_label==0)
		{
			if(thread_num>1)
			{
				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
				uint64_t unit_long=Task_Size_Single_Thread;
				struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp;
				p_hv_tmp=(struct TwoDimenHashV_Final_256bitKmer**)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer*)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					p_hv_tmp[i]=(struct TwoDimenHashV_Final_256bitKmer*)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer)*unit_long);
					for(uint32_t j=0;j<unit_long;j++)
					{
						p_hv_tmp[i][j].tmp.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para2.kmer64Len);
					}
				}
				uint64_t* result_len;
				result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
				struct Parallel_para_Final01_hash_256bitKmer* para_alloc;
				para_alloc=(struct Parallel_para_Final01_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final01_hash_256bitKmer)*thread_num);
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
					para_alloc[i].para1=para1;
					para_alloc[i].para2=para2;

					para_alloc[i].p_candidate_pos=p_candidate_pos[ref_i];
					para_alloc[i].template_64=template_64;
					para_alloc[i].pos_label=tmp.pos_label;
					para_alloc[i].p_candidate_current=p_candidate_current;
				}

				struct Parallel_para_Final02_hash_256bitKmer* para_assi;
				para_assi=(struct Parallel_para_Final02_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final02_hash_256bitKmer)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_assi[i].seq=seq;
					para_assi[i].p_root_hashtable=p_root_hashtable;
					para_assi[i].kmer_len1=original_kmer_length;
					para_assi[i].kmer_len2=extend_length;
					para_assi[i].InOut_label=InOut_label;
					para_assi[i].result_len=result_len;
					para_assi[i].p_hv_tmp=p_hv_tmp;
					para_assi[i].thread_id=i;
					para_assi[i].thread_num=thread_num;
					para_assi[i].para=para2;
				}

				uint64_t start=0;
				for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<=seq_length-(original_kmer_length+extend_length);loop_num++)
				{
					start=loop_num*(thread_num*unit_long);
					for(uint32_t j=0;j<thread_num;j++)
					{
						*(para_alloc[j].result_len)=0;
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
							if(pthread_create(t+j, NULL, Allocate_HashTable_Final_256bitKmer, (void*)(para_alloc+j))!=0)
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
						if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length+extend_length))
						{
							if(pthread_create(t+i, NULL, Assign_HashTable_Final_256bitKmer, (void*)(para_assi+i))!=0)
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
				}
				free(t);
				for(uint32_t i=0;i<thread_num;i++)
				{
					for(uint32_t j=0;j<unit_long;j++)
					{
						free(p_hv_tmp[i][j].tmp.hashValue);
					}
				}
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
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						uint64_t tmp_pos_divide=j/64;
						uint64_t tmp_pos_mod=j%64;
						uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
						if(tmp.template_64[tmp_pos_mod]==tmp_template)
						{
							if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
							{
								p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

								cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
								for(uint32_t k=0;k<para2.kmer64Len;k++)
								{
									c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
								}
								uint8_t ad_temp=0;
								Transition(temp_seq[j-1],'0',&ad_temp);
								c_tmp_hashtable.ad=ad_temp;
								bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
							}
						}
					}
				}
				else if(tmp.pos_label==2)
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						uint64_t tmp_pos_divide=j/64;
						uint64_t tmp_pos_mod=j%64;
						uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
						if(tmp.template_64[tmp_pos_mod]==tmp_template)
						{
							if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
							{
								cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
								for(uint32_t k=0;k<para2.kmer64Len;k++)
								{
									c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
								}
								uint8_t ad_temp=0;
								Transition(temp_seq[j-1],'0',&ad_temp);
								c_tmp_hashtable.ad=ad_temp;
								bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
							}
						}
					}
				}
				else if(tmp.pos_label==1)
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
						{
							uint64_t tmp_pos_divide=j/64;
							uint64_t tmp_pos_mod=j%64;
							p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

							cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
							for(uint32_t k=0;k<para2.kmer64Len;k++)
							{
								c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
							}
							uint8_t ad_temp=0;
							Transition(temp_seq[j-1],'0',&ad_temp);
							c_tmp_hashtable.ad=ad_temp;
							bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
						}
					}
				}
				else
				{
					for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
						{
							cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
							for(uint32_t k=0;k<para2.kmer64Len;k++)
							{
								c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
							}
							uint8_t ad_temp=0;
							Transition(temp_seq[j-1],'0',&ad_temp);
							c_tmp_hashtable.ad=ad_temp;
							bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
						}
					}
				}
			}
		}
		else
		{
			if(thread_num>1)
			{
				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
				uint64_t unit_long=Task_Size_Single_Thread;
				struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp;
				p_hv_tmp=(struct TwoDimenHashV_Final_256bitKmer**)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer*)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					p_hv_tmp[i]=(struct TwoDimenHashV_Final_256bitKmer*)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer)*unit_long);
					for(uint32_t j=0;j<unit_long;j++)
					{
						p_hv_tmp[i][j].tmp.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para2.kmer64Len);
					}
				}
				uint64_t* result_len;
				result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
				struct Parallel_para_Final01_hash_256bitKmer* para_alloc;
				para_alloc=(struct Parallel_para_Final01_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final01_hash_256bitKmer)*thread_num);
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
					para_alloc[i].para1=para1;
					para_alloc[i].para2=para2;

					para_alloc[i].p_candidate_pos=p_candidate_pos[ref_i];
					para_alloc[i].template_64=template_64;
					para_alloc[i].pos_label=tmp.pos_label;
					para_alloc[i].p_candidate_current=p_candidate_current;
				}
				struct Parallel_para_Final02_hash_256bitKmer* para_assi;
				para_assi=(struct Parallel_para_Final02_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final02_hash_256bitKmer)*thread_num);
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_assi[i].seq=seq;
					para_assi[i].p_root_hashtable=p_root_hashtable;
					para_assi[i].kmer_len1=original_kmer_length;
					para_assi[i].kmer_len2=extend_length;
					para_assi[i].InOut_label=InOut_label;
					para_assi[i].result_len=result_len;
					para_assi[i].p_hv_tmp=p_hv_tmp;
					para_assi[i].thread_id=i;
					para_assi[i].thread_num=thread_num;
					para_assi[i].para=para2;
				}

				uint64_t start=0;
				for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length-original_kmer_length;loop_num++)
				{
					start=loop_num*thread_num*unit_long;
					for(uint32_t j=0;j<thread_num;j++)
					{
						*(para_alloc[j].result_len)=0;
						para_alloc[j].task_start_pos=start+j*unit_long;
						if(para_alloc[j].task_start_pos<seq_length-(original_kmer_length))
						{
							if(start+(j+1)*unit_long>=seq_length-(original_kmer_length))
							{
								//para_alloc[j].task_len=seq_length-(original_kmer_length+extend_length)-para_alloc[j].task_start_pos+1;
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
						if(para_alloc[j].task_start_pos<seq_length-(original_kmer_length))
						{
							if(pthread_create(t+j, NULL, Allocate_HashTable_Final_256bitKmer, (void*)(para_alloc+j))!=0)
							{
								cout << "error!" << endl;
							}
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
						{
							pthread_join(t[i], NULL);
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
						{
							if(pthread_create(t+i, NULL, Assign_HashTable_Final_256bitKmer, (void*)(para_assi+i))!=0)
							{
								cout << "error!" << endl;
							}
						}
					}
					for(uint32_t i=0;i<thread_num;i++)
					{
						if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
						{
							pthread_join(t[i], NULL);
						}
					}
				}
				free(t);
				for(uint32_t i=0;i<thread_num;i++)
				{
					for(uint32_t j=0;j<unit_long;j++)
					{
						free(p_hv_tmp[i][j].tmp.hashValue);
					}
				}
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
				if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
				{
					if(tmp.pos_label==1||tmp.pos_label==3)
					{
						uint64_t tmp_pos=extend_length+original_kmer_length-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];
					}

					cal_hash_value_directly_256bit(temp_seq,seq_k_currentF,para2);
					for(uint32_t k=0;k<para2.kmer64Len;k++)
					{
						c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
					}
					c_tmp_hashtable.arrayID=0;
					uint8_t ad_temp=0;
					Transition('0',temp_seq[extend_length+original_kmer_length],&ad_temp);
					c_tmp_hashtable.ad=ad_temp;
					bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
				}
				if(tmp.pos_label==3)
				{
					for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						uint64_t tmp_pos=j+original_kmer_length-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
						if(tmp.template_64[tmp_pos_mod]==tmp_template)
						{
							if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
							{
								p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

								cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
								for(uint32_t k=0;k<para2.kmer64Len;k++)
								{
									c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
								}
								c_tmp_hashtable.arrayID=0;
								uint8_t ad_temp=0;
								Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
								c_tmp_hashtable.ad=ad_temp;
								bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
							}
						}
					}
				}
				else if(tmp.pos_label==2)
				{
					for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						uint64_t tmp_pos=j+original_kmer_length-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
						if(tmp.template_64[tmp_pos_mod]==tmp_template)
						{
							if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
							{
								cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
								for(uint32_t k=0;k<para2.kmer64Len;k++)
								{
									c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
								}
								c_tmp_hashtable.arrayID=0;
								uint8_t ad_temp=0;
								Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
								c_tmp_hashtable.ad=ad_temp;
								bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
							}
						}
					}
				}
				else if(tmp.pos_label==1)
				{
					for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
						{
							uint64_t tmp_pos=j+original_kmer_length-1;
							uint64_t tmp_pos_divide=tmp_pos/64;
							uint64_t tmp_pos_mod=tmp_pos%64;
							p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

							cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
							for(uint32_t k=0;k<para2.kmer64Len;k++)
							{
								c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
							}
							c_tmp_hashtable.arrayID=0;
							uint8_t ad_temp=0;
							Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
							c_tmp_hashtable.ad=ad_temp;
							bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
						}
					}
				}
				else
				{
					for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
					{
						cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
						if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
						{
							cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
							for(uint32_t k=0;k<para2.kmer64Len;k++)
							{
								c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
							}
							c_tmp_hashtable.arrayID=0;
							uint8_t ad_temp=0;
							Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
							c_tmp_hashtable.ad=ad_temp;
							bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
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
		if(rc_label==1)
        {
            rc(&seq,seq_length);
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
                uint64_t x;
                x=fread(p_candidate_current,sizeof(uint64_t),seq_length/64+((seq_length%64)?1:0),fp_pos);
                if(x!=seq_length/64+((seq_length%64)?1:0))
                {
                    cout << "error: read current candidate failed!" << endl;
                }
                fclose(fp_pos);
            }

            struct timeval tvs,tve;
            gettimeofday(&tvs,NULL);
            uint64_t seq_k_current[4];
            uint64_t seq_k_currentF[4];
            if(InOut_label==0)
            {
                if(thread_num>1)
                {
                    pthread_t* t;
                    t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                    uint64_t unit_long=Task_Size_Single_Thread;
                    struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp;
                    p_hv_tmp=(struct TwoDimenHashV_Final_256bitKmer**)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer*)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        p_hv_tmp[i]=(struct TwoDimenHashV_Final_256bitKmer*)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer)*unit_long);
                        for(uint32_t j=0;j<unit_long;j++)
                        {
                            p_hv_tmp[i][j].tmp.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para2.kmer64Len);
                        }
                    }
                    uint64_t* result_len;
                    result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
                    struct Parallel_para_Final01_hash_256bitKmer* para_alloc;
                    para_alloc=(struct Parallel_para_Final01_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final01_hash_256bitKmer)*thread_num);
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
                        para_alloc[i].para1=para1;
                        para_alloc[i].para2=para2;

                        para_alloc[i].p_candidate_pos=p_candidate_pos[ref_i];
                        para_alloc[i].template_64=template_64;
                        para_alloc[i].pos_label=tmp.pos_label;
                        para_alloc[i].p_candidate_current=p_candidate_current;
                    }

                    struct Parallel_para_Final02_hash_256bitKmer* para_assi;
                    para_assi=(struct Parallel_para_Final02_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final02_hash_256bitKmer)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_assi[i].seq=seq;
                        para_assi[i].p_root_hashtable=p_root_hashtable;
                        para_assi[i].kmer_len1=original_kmer_length;
                        para_assi[i].kmer_len2=extend_length;
                        para_assi[i].InOut_label=InOut_label;
                        para_assi[i].result_len=result_len;
                        para_assi[i].p_hv_tmp=p_hv_tmp;
                        para_assi[i].thread_id=i;
                        para_assi[i].thread_num=thread_num;
                        para_assi[i].para=para2;
                    }

                    uint64_t start=0;
                    for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<=seq_length-(original_kmer_length+extend_length);loop_num++)
                    {
                        start=loop_num*(thread_num*unit_long);
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            *(para_alloc[j].result_len)=0;
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
                                if(pthread_create(t+j, NULL, Allocate_HashTable_Final_256bitKmer, (void*)(para_alloc+j))!=0)
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
                            if(para_alloc[i].task_start_pos<=seq_length-(original_kmer_length+extend_length))
                            {
                                if(pthread_create(t+i, NULL, Assign_HashTable_Final_256bitKmer, (void*)(para_assi+i))!=0)
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
                    }
                    free(t);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        for(uint32_t j=0;j<unit_long;j++)
                        {
                            free(p_hv_tmp[i][j].tmp.hashValue);
                        }
                    }
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
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            uint64_t tmp_pos_divide=j/64;
                            uint64_t tmp_pos_mod=j%64;
                            uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
                            if(tmp.template_64[tmp_pos_mod]==tmp_template)
                            {
                                if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                                {
                                    p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

                                    cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
                                    for(uint32_t k=0;k<para2.kmer64Len;k++)
                                    {
                                        c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                    }
                                    uint8_t ad_temp=0;
                                    Transition(temp_seq[j-1],'0',&ad_temp);
                                    c_tmp_hashtable.ad=ad_temp;
                                    bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==2)
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            uint64_t tmp_pos_divide=j/64;
                            uint64_t tmp_pos_mod=j%64;
                            uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
                            if(tmp.template_64[tmp_pos_mod]==tmp_template)
                            {
                                if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                                {
                                    cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
                                    for(uint32_t k=0;k<para2.kmer64Len;k++)
                                    {
                                        c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                    }
                                    uint8_t ad_temp=0;
                                    Transition(temp_seq[j-1],'0',&ad_temp);
                                    c_tmp_hashtable.ad=ad_temp;
                                    bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==1)
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                            {
                                uint64_t tmp_pos_divide=j/64;
                                uint64_t tmp_pos_mod=j%64;
                                p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

                                cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
                                for(uint32_t k=0;k<para2.kmer64Len;k++)
                                {
                                    c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                }
                                uint8_t ad_temp=0;
                                Transition(temp_seq[j-1],'0',&ad_temp);
                                c_tmp_hashtable.ad=ad_temp;
                                bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                            }
                        }
                    }
                    else
                    {
                        for(uint64_t j=1;j<seq_length-(original_kmer_length+extend_length)+1;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                            {
                                cal_hash_value_directly_256bit(temp_seq+j,seq_k_currentF,para2);
                                for(uint32_t k=0;k<para2.kmer64Len;k++)
                                {
                                    c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                }
                                uint8_t ad_temp=0;
                                Transition(temp_seq[j-1],'0',&ad_temp);
                                c_tmp_hashtable.ad=ad_temp;
                                bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                            }
                        }
                    }
                }
            }
            else
            {
                if(thread_num>1)
                {
                    pthread_t* t;
                    t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                    uint64_t unit_long=Task_Size_Single_Thread;
                    struct TwoDimenHashV_Final_256bitKmer** p_hv_tmp;
                    p_hv_tmp=(struct TwoDimenHashV_Final_256bitKmer**)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer*)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        p_hv_tmp[i]=(struct TwoDimenHashV_Final_256bitKmer*)malloc(sizeof(struct TwoDimenHashV_Final_256bitKmer)*unit_long);
                        for(uint32_t j=0;j<unit_long;j++)
                        {
                            p_hv_tmp[i][j].tmp.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para2.kmer64Len);
                        }
                    }
                    uint64_t* result_len;
                    result_len=(uint64_t*)malloc(sizeof(uint64_t)*thread_num);
                    struct Parallel_para_Final01_hash_256bitKmer* para_alloc;
                    para_alloc=(struct Parallel_para_Final01_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final01_hash_256bitKmer)*thread_num);
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
                        para_alloc[i].para1=para1;
                        para_alloc[i].para2=para2;

                        para_alloc[i].p_candidate_pos=p_candidate_pos[ref_i];
                        para_alloc[i].template_64=template_64;
                        para_alloc[i].pos_label=tmp.pos_label;
                        para_alloc[i].p_candidate_current=p_candidate_current;
                    }
                    struct Parallel_para_Final02_hash_256bitKmer* para_assi;
                    para_assi=(struct Parallel_para_Final02_hash_256bitKmer*)malloc(sizeof(struct Parallel_para_Final02_hash_256bitKmer)*thread_num);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_assi[i].seq=seq;
                        para_assi[i].p_root_hashtable=p_root_hashtable;
                        para_assi[i].kmer_len1=original_kmer_length;
                        para_assi[i].kmer_len2=extend_length;
                        para_assi[i].InOut_label=InOut_label;
                        para_assi[i].result_len=result_len;
                        para_assi[i].p_hv_tmp=p_hv_tmp;
                        para_assi[i].thread_id=i;
                        para_assi[i].thread_num=thread_num;
                        para_assi[i].para=para2;
                    }

                    uint64_t start=0;
                    for(uint32_t loop_num=0;loop_num*(thread_num*unit_long)<seq_length-original_kmer_length;loop_num++)
                    {
                        start=loop_num*thread_num*unit_long;
                        for(uint32_t j=0;j<thread_num;j++)
                        {
                            *(para_alloc[j].result_len)=0;
                            para_alloc[j].task_start_pos=start+j*unit_long;
                            if(para_alloc[j].task_start_pos<seq_length-(original_kmer_length))
                            {
                                if(start+(j+1)*unit_long>=seq_length-(original_kmer_length))
                                {
                                    //para_alloc[j].task_len=seq_length-(original_kmer_length+extend_length)-para_alloc[j].task_start_pos+1;
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
                            if(para_alloc[j].task_start_pos<seq_length-(original_kmer_length))
                            {
                                if(pthread_create(t+j, NULL, Allocate_HashTable_Final_256bitKmer, (void*)(para_alloc+j))!=0)
                                {
                                    cout << "error!" << endl;
                                }
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
                            {
                                pthread_join(t[i], NULL);
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
                            {
                                if(pthread_create(t+i, NULL, Assign_HashTable_Final_256bitKmer, (void*)(para_assi+i))!=0)
                                {
                                    cout << "error!" << endl;
                                }
                            }
                        }
                        for(uint32_t i=0;i<thread_num;i++)
                        {
                            if(para_alloc[i].task_start_pos<seq_length-(original_kmer_length))
                            {
                                pthread_join(t[i], NULL);
                            }
                        }
                    }
                    free(t);
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        for(uint32_t j=0;j<unit_long;j++)
                        {
                            free(p_hv_tmp[i][j].tmp.hashValue);
                        }
                    }
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
                    if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                    {
                        if(tmp.pos_label==1||tmp.pos_label==3)
                        {
                            uint64_t tmp_pos=extend_length+original_kmer_length-1;
                            uint64_t tmp_pos_divide=tmp_pos/64;
                            uint64_t tmp_pos_mod=tmp_pos%64;
                            p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];
                        }

                        cal_hash_value_directly_256bit(temp_seq,seq_k_currentF,para2);
                        for(uint32_t k=0;k<para2.kmer64Len;k++)
                        {
                            c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                        }
                        c_tmp_hashtable.arrayID=0;
                        uint8_t ad_temp=0;
                        Transition('0',temp_seq[extend_length+original_kmer_length],&ad_temp);
                        c_tmp_hashtable.ad=ad_temp;
                        bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                    }
                    if(tmp.pos_label==3)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            uint64_t tmp_pos=j+original_kmer_length-1;
                            uint64_t tmp_pos_divide=tmp_pos/64;
                            uint64_t tmp_pos_mod=tmp_pos%64;
                            uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
                            if(tmp.template_64[tmp_pos_mod]==tmp_template)
                            {
                                if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                                {
                                    p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

                                    cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
                                    for(uint32_t k=0;k<para2.kmer64Len;k++)
                                    {
                                        c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                    }
                                    c_tmp_hashtable.arrayID=0;
                                    uint8_t ad_temp=0;
                                    Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
                                    c_tmp_hashtable.ad=ad_temp;
                                    bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==2)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            uint64_t tmp_pos=j+original_kmer_length-1;
                            uint64_t tmp_pos_divide=tmp_pos/64;
                            uint64_t tmp_pos_mod=tmp_pos%64;
                            uint64_t tmp_template=tmp.template_64[tmp_pos_mod]&p_candidate_current[tmp_pos_divide];
                            if(tmp.template_64[tmp_pos_mod]==tmp_template)
                            {
                                if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                                {
                                    cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
                                    for(uint32_t k=0;k<para2.kmer64Len;k++)
                                    {
                                        c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                    }
                                    c_tmp_hashtable.arrayID=0;
                                    uint8_t ad_temp=0;
                                    Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
                                    c_tmp_hashtable.ad=ad_temp;
                                    bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                                }
                            }
                        }
                    }
                    else if(tmp.pos_label==1)
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                            {
                                uint64_t tmp_pos=j+original_kmer_length-1;
                                uint64_t tmp_pos_divide=tmp_pos/64;
                                uint64_t tmp_pos_mod=tmp_pos%64;
                                p_candidate_pos[ref_i][tmp_pos_divide]=p_candidate_pos[ref_i][tmp_pos_divide]|template_64[tmp_pos_mod];

                                cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
                                for(uint32_t k=0;k<para2.kmer64Len;k++)
                                {
                                    c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                }
                                c_tmp_hashtable.arrayID=0;
                                uint8_t ad_temp=0;
                                Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
                                c_tmp_hashtable.ad=ad_temp;
                                bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
                            }
                        }
                    }
                    else
                    {
                        for(uint64_t j=extend_length+1;j<seq_length-original_kmer_length;j++)//check each pos without the first in the case of "branch_in"
                        {
                            cal_hash_value_indirectly_256bit(temp_seq+j,seq_k_current,seq_k_current,para1);
                            if(search256BitHashTable(p_root,seq_k_current,para1)!=0)
                            {
                                cal_hash_value_directly_256bit(temp_seq+j-extend_length,seq_k_currentF,para2);
                                for(uint32_t k=0;k<para2.kmer64Len;k++)
                                {
                                    c_tmp_hashtable.hashValue[k]=seq_k_currentF[k];
                                }
                                c_tmp_hashtable.arrayID=0;
                                uint8_t ad_temp=0;
                                Transition('0',temp_seq[j+original_kmer_length],&ad_temp);
                                c_tmp_hashtable.ad=ad_temp;
                                bit256insertHashFTable(p_root_hashtable,c_tmp_hashtable,para2);
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
	cout <<"write hashtable is over!"<< endl;
	/**********************fill the two-dimensional hash table!*********************/

	/*******************check each char space and output the results****************/
	int statics_branch_kmer_number=0;
	FILE *out_hash;
	struct NodeBit * p_leaf_hashtable;

	char* outputFileName;
	outputFileName=(char*)malloc(sizeof(char)*11);
	struct para_getN tmp_outname;
	tmp_outname.kmerlen=original_kmer_length+extend_length;
	tmp_outname.FileName=outputFileName;
	tmp_outname.InOut_label=InOut_label;
	tmp_outname.isMiddle=0;
	tmp_outname.isposition=0;
	getFileName(tmp_outname);
	out_hash=fopen(outputFileName,"ab+");

	FILE *out_hash_ad;
	char* outputFileName_ad;
	struct para_getN tmp_ad;
	if(ad_label==1)
	{
		outputFileName_ad=(char*)malloc(sizeof(char)*13);
		tmp_ad.kmerlen=original_kmer_length+extend_length;
		tmp_ad.FileName=outputFileName_ad;
		tmp_ad.InOut_label=InOut_label;
		tmp_ad.isMiddle=2;
		tmp_ad.isposition=0;
		getFileName(tmp_ad);
		out_hash_ad=fopen(outputFileName_ad,"ab+");
	}


	bit256linkHashFTable(p_root_hashtable);
	p_leaf_hashtable=bit256Find_Minimal_Node_HashFTable(p_root_hashtable);

	if(thread_num>1)
	{
		struct NodeBit *p_leaf_tmp=p_leaf_hashtable;
		uint64_t total_leaf_num=0;
		while(p_leaf_tmp!=NULL)
		{
			total_leaf_num++;
			p_leaf_tmp=p_leaf_tmp->brother;
		}
		uint64_t unit_leaf_num=total_leaf_num/thread_num+1;

		uint64_t **reslt;
		uint32_t **reslt_c;
		reslt=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_c=(uint32_t **)malloc(sizeof(uint32_t *)*thread_num);

		uint64_t * reslt_size;
		reslt_size=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);


		struct reslt_para_plustree_final_256bitKmer* tmp;
		tmp=(struct reslt_para_plustree_final_256bitKmer*)malloc(sizeof(struct reslt_para_plustree_final_256bitKmer)*thread_num);

		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		p_leaf_tmp=p_leaf_hashtable;
		for(uint32_t j=0;j<thread_num-1;j++)
		{
			tmp[j].p_root_start=p_leaf_tmp;
			for(uint64_t loop_i=0;loop_i<unit_leaf_num;loop_i++)
			{
				p_leaf_tmp=p_leaf_tmp->brother;
			}
			tmp[j].p_root_end=p_leaf_tmp;
		}
		tmp[thread_num-1].p_root_start=p_leaf_tmp;
		tmp[thread_num-1].p_root_end=NULL;
		for(uint32_t j=0;j<thread_num;j++)
		{
			tmp[j].reslt=&(reslt[j]);
			tmp[j].reslt_c=&(reslt_c[j]);
			tmp[j].reslt_size=reslt_size+j;
			tmp[j].InOut_label=InOut_label;
			tmp[j].ad_label=ad_label;
			tmp[j].para=para2;
			if(pthread_create(t+j, NULL, judge_thread_plustree_final_256bitKmer, (void*)(tmp+j))!=0)
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
				fwrite(reslt[j]+k*para2.kmer64Len,sizeof(uint64_t),para2.kmer64Len,out_hash);
				if(ad_label==1)
				{
					fwrite(&(reslt_c[j][k]),sizeof(uint8_t),1,out_hash_ad);
				}
				statics_branch_kmer_number++;
			}
		}

		for(uint32_t j=0;j<thread_num;j++)
		{
			free(reslt[j]);
			if(ad_label==1)
			{
				free(reslt_c[j]);
			}
		}

		free(reslt);
		free(reslt_c);
		free(reslt_size);
		free(tmp);
		free(t);
	}
	else
	{
		while(p_leaf_hashtable!=NULL)
		{
			for(uint32_t i=0;i<p_leaf_hashtable->Node_Size;i++)
			{
				if(Judge_4bit(p_leaf_hashtable->data[i].ad,InOut_label))
				{
					fwrite(p_leaf_hashtable->data[i].hashValue,sizeof(uint64_t),para2.kmer64Len,out_hash);

					if(ad_label==1)
					{
						uint8_t c_tmp;
						if(InOut_label==0)
						{
							c_tmp=(p_leaf_hashtable->data[i].ad>>4)&15;
						}
						else
						{
							c_tmp=(p_leaf_hashtable->data[i].ad)&15;
						}
						fwrite(&(c_tmp),sizeof(uint8_t),1,out_hash_ad);
					}

					statics_branch_kmer_number++;
				}
			}
			p_leaf_hashtable=p_leaf_hashtable->brother;
		}
	}
	fclose(out_hash);
	free(outputFileName);
	if(ad_label==1)
	{
		fclose(out_hash_ad);
		free(outputFileName_ad);
	}

	cout <<"check each char space and output the results is over!"<<endl;
	cout <<"the # of branched kmer is "<<statics_branch_kmer_number<<endl;
	/*******************check each char space and output the results****************/

	free(c_tmp_hashtable.hashValue);
	bit256freeHashFTable(p_root_hashtable,para2);
//	destory_tree_bit(p_root_hashtable,para2);
	freeHash256BitTable(p_root);
}
/***************************************Final for 256bitKmer******************/
