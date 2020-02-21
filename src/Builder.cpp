/*
 * Builder.cpp
 *
 *  Created on: May 15, 2018
 *      Author: bio
 */

#include "Builder.h"
#include "Hash.h"

//1)method 1. the parallel method in method 1 cost too much memory.
void *Find_unipath_para(void *arg)
{
	struct para_unipath_bit* tmp=(struct para_unipath_bit*)arg;

	uint32_t kmer_len=tmp->kmer_len;
	char * p_tmp=tmp->p_tmp;
	struct NodeBit * p_Branch_root=tmp->p_Branch_root;
	uint64_t start_t=tmp->start;
	uint64_t end_t=tmp->end;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	struct bit256KmerPara para=tmp->para;

	uint64_t* p_candidate_pos_in=tmp->p_candidate_pos_in;
	uint64_t* p_candidate_pos_out=tmp->p_candidate_pos_out;
	uint64_t* template_64=tmp->template_64;
	uint64_t pos_label=tmp->pos_label;

	uint64_t*reslt_c;
	uint64_t *reslt;
	uint64_t reslt_size=0;
	uint64_t *reslt_bran_pos;
	uint64_t reslt_size_bran_pos=0;
	uint64_t reslt_Max_size;
	uint64_t reslt_Max_size_bran_pos;
	reslt=(uint64_t *)malloc(sizeof(uint64_t)*pow(2,10));
	reslt_c=(uint64_t *)malloc(sizeof(uint64_t)*pow(2,10));
	reslt_Max_size=pow(2,10);

	reslt_bran_pos=(uint64_t *)malloc(sizeof(uint64_t)*pow(2,10));
	reslt_Max_size_bran_pos=pow(2,10);

	uint64_t start_l;
	uint64_t start,end;
	uint64_t seq_k_current[4];
	if(thread_id==0)
	{
		start=start_t;
		end=start_t;
		cal_hash_value_directly_256bit(p_tmp+start_t,seq_k_current,para);
	}
	else
	{
		cal_hash_value_directly_256bit(p_tmp+start_t,seq_k_current,para);
		if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
		{
			start=start_t+1;
			end=start_t+1;
		}
		else
		{
			for(uint64_t i=1+start_t;i<=end_t;i++)
			{
				cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);
				if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
				{
					start=i+1;
					end=i+1;
					break;
				}
			}
		}
	}
	uint64_t label_pre=1;
	if(thread_id==0)
	{
		start_l=start+1;
	}
	else
	{
		start_l=start;
	}

	for(uint64_t i=start_l;i<=end_t;i++)
	{
		cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);
		if(label_pre==1)
		{
			//save the bran-unip edge position
			if(reslt_size_bran_pos<reslt_Max_size_bran_pos)
			{
				reslt_bran_pos[reslt_size_bran_pos]=i-1;
				reslt_size_bran_pos++;
			}
			else
			{
				reslt_bran_pos=(uint64_t*)realloc(reslt_bran_pos,sizeof(uint64_t)*(reslt_Max_size_bran_pos+1));
				reslt_Max_size_bran_pos++;

				reslt_bran_pos[reslt_size_bran_pos]=i-1;
				reslt_size_bran_pos++;
			}
		}

		if(pos_label==2||pos_label==3)
		{
			uint64_t tmp_pos_divide_in=i/64;
			uint64_t tmp_pos_mod_in=i%64;
			uint64_t tmp_template=template_64[tmp_pos_mod_in]&p_candidate_pos_in[tmp_pos_divide_in];
			if(template_64[tmp_pos_mod_in]==tmp_template)
			{
				if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
				{
					//label_pre=0 is have no pre branch
					if(label_pre==0)
					{
						end=i;
						if(reslt_size<reslt_Max_size)
						{
							reslt[reslt_size]=start;
							reslt_c[reslt_size]=end;
							reslt_size++;
						}
						else
						{
							reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*(reslt_Max_size+1));
							reslt_c=(uint64_t*)realloc(reslt_c,sizeof(uint64_t)*(reslt_Max_size+1));
							reslt_Max_size=reslt_Max_size+1;

							reslt[reslt_size]=start;
							reslt_c[reslt_size]=end;
							reslt_size++;
						}
					}
					start=i+1;
					label_pre=1;
				}
				else
				{
					label_pre=0;
				}
			}
			else
			{
				uint64_t tmp_pos=i+kmer_len-1;
				uint64_t tmp_pos_divide=tmp_pos/64;
				uint64_t tmp_pos_mod=tmp_pos%64;
				uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_pos_out[tmp_pos_divide];
				if(template_64[tmp_pos_mod]==tmp_template)
				{
					if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
					{
						//label_pre=0 is have no pre branch
						if(label_pre==0)
						{
							end=i;
							if(reslt_size<reslt_Max_size)
							{
								reslt[reslt_size]=start;
								reslt_c[reslt_size]=end;
								reslt_size++;
							}
							else
							{
								reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*(reslt_Max_size+1));
								reslt_c=(uint64_t*)realloc(reslt_c,sizeof(uint64_t)*(reslt_Max_size+1));
								reslt_Max_size=reslt_Max_size+1;

								reslt[reslt_size]=start;
								reslt_c[reslt_size]=end;
								reslt_size++;
							}
						}
						start=i+1;
						label_pre=1;
					}
					else
					{
						label_pre=0;
					}
				}
				else
				{
					label_pre=0;
				}
			}
		}
		else
		{
			if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
			{
				//label_pre=0 is have no pre branch
				if(label_pre==0)
				{
					end=i;
					if(reslt_size<reslt_Max_size)
					{
						reslt[reslt_size]=start;
						reslt_c[reslt_size]=end;
						reslt_size++;
					}
					else
					{
						reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*(reslt_Max_size+1));
						reslt_c=(uint64_t*)realloc(reslt_c,sizeof(uint64_t)*(reslt_Max_size+1));
						reslt_Max_size=reslt_Max_size+1;

						reslt[reslt_size]=start;
						reslt_c[reslt_size]=end;
						reslt_size++;
					}
				}
				start=i+1;
				label_pre=1;
			}
			else
			{
				label_pre=0;
			}
		}
	}

	if(thread_id!=thread_num-1)
	{
		if(label_pre==1)
		{
			if(reslt_size_bran_pos<reslt_Max_size_bran_pos)
			{
				reslt_bran_pos[reslt_size_bran_pos]=end_t;
				reslt_size_bran_pos++;
			}
			else
			{
				reslt_bran_pos=(uint64_t*)realloc(reslt_bran_pos,sizeof(uint64_t)*(reslt_Max_size_bran_pos+1));
				reslt_Max_size_bran_pos++;

				reslt_bran_pos[reslt_size_bran_pos]=end_t;
				reslt_size_bran_pos++;
			}
		}
		uint64_t i=1;
		while(1)
		{
			cal_hash_value_indirectly_256bit(p_tmp+end_t+i,seq_k_current,seq_k_current,para);

			if(MappingHashValueToID_bit(p_Branch_root,seq_k_current,para)!=-1)
			{
				//label_pre=0 is have no pre branch
				if(label_pre==0)
				{
					end=i+end_t;
					if(reslt_size<reslt_Max_size)
					{
						reslt[reslt_size]=start;
						reslt_c[reslt_size]=end;
						reslt_size++;
					}
					else
					{
						reslt=(uint64_t*)realloc(reslt,sizeof(uint64_t)*(reslt_Max_size+1));
						reslt_c=(uint64_t*)realloc(reslt_c,sizeof(uint64_t)*(reslt_Max_size+1));
						reslt_Max_size=reslt_Max_size+1;

						reslt[reslt_size]=start;
						reslt_c[reslt_size]=end;
						reslt_size++;
					}
				}
				start=i+1;
				label_pre=1;
				break;
			}
			else
			{
				label_pre=0;
			}

			i++;
		}
	}

	*(tmp->reslt_bran_pos)=reslt_bran_pos;
	*(tmp->reslt_size_bran_pos)=reslt_size_bran_pos;
	*(tmp->reslt_size)=reslt_size;
	*(tmp->reslt)=reslt;
	*(tmp->reslt_c)=reslt_c;
	return NULL;
}
void Insert_pos_256bitKmer(struct NodeBit** p_root,struct nodeBit insert_data,struct bit256KmerPara para)
{
	struct NodeBit*p_root_tmp= * p_root;
	if(p_root_tmp==NULL)
	{
		cout << "error:the plus tree is used without initialized!" << endl;
		return;
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
//			p_root_tmp->data[pos].ad=p_root_tmp->data[pos].ad|insert_data.ad;
			if(p_root_tmp->data[pos].p_posList_length==0)
			{
				uint32_t tmp_len=p_root_tmp->data[pos].p_posList_length;
				p_root_tmp->data[pos].p_posList=(uint64_t*)malloc(sizeof(uint64_t)*(tmp_len+1));
				p_root_tmp->data[pos].p_posList_length++;

				p_root_tmp->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
			}
			else
			{
				uint32_t tmp_len=p_root_tmp->data[pos].p_posList_length;
				p_root_tmp->data[pos].p_posList=(uint64_t*)realloc(p_root_tmp->data[pos].p_posList,sizeof(uint64_t)*(tmp_len+1));
				p_root_tmp->data[pos].p_posList_length++;

				p_root_tmp->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
			}
			return;
		}
		else if(pos_label==0)
		{
			pos=p_root_tmp->Node_Size;
		}
	}

	if(p_root_tmp->Node_Size==M)
	{
		Divide_Node_bit(p_root,p_root_tmp,insert_data,pos,NULL,0,1,para);
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
//			p_root_tmp->data[l+1].ad=p_root_tmp->data[l].ad;
			p_root_tmp->data[l+1].p_posList_length=p_root_tmp->data[l].p_posList_length;
			p_root_tmp->data[l+1].p_posList=p_root_tmp->data[l].p_posList;
		}

		p_root_tmp->data[pos].p_posList=NULL;
		p_root_tmp->data[pos].p_posList_length=0;

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
//		p_root_tmp->data[pos].ad=insert_data.ad;
		if(p_root_tmp->data[pos].p_posList_length==0)
		{
			uint32_t tmp_len=p_root_tmp->data[pos].p_posList_length;
			p_root_tmp->data[pos].p_posList=(uint64_t*)malloc(sizeof(uint64_t)*(tmp_len+1));
			p_root_tmp->data[pos].p_posList_length++;

			p_root_tmp->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
		}
		else
		{
			uint32_t tmp_len=p_root_tmp->data[pos].p_posList_length;
			p_root_tmp->data[pos].p_posList=(uint64_t*)realloc(p_root_tmp->data[pos].p_posList,sizeof(uint64_t)*(tmp_len+1));
			p_root_tmp->data[pos].p_posList_length++;

			p_root_tmp->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
		}

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
void FindUniPath(struct para_findUnipath tmps)
{
	char* seq=tmps.seq;
	uint32_t kmer_len=tmps.original_kmer_length;
	uint32_t thread_num=tmps.thread_num;
	uint32_t pos_label=tmps.pos_label;
	uint32_t rc_label=tmps.rc_label;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	if(para.remainer1to64==0)
	{
		para.codefor1=0;
	}
	else
	{
		para.codefor1=0;
		for(uint32_t i=0;i<para.remainer1to64-1;i++)
		{
			para.codefor1=para.codefor1|1;
			para.codefor1=para.codefor1<<1;
		}
		para.codefor1=para.codefor1|1;
	}

	struct bit256KmerPara para_edge;
	para_edge.kmer1Len=kmer_len*2+2;
	para_edge.remainer1to64=para_edge.kmer1Len%64;
	para_edge.kmer64Len=para_edge.kmer1Len/64+(para_edge.remainer1to64?1:0);
	if(para_edge.remainer1to64==0)
	{
		para_edge.codefor1=0;
	}
	else
	{
		para_edge.codefor1=0;
		for(uint32_t i=0;i<para_edge.remainer1to64-1;i++)
		{
			para_edge.codefor1=para_edge.codefor1|1;
			para_edge.codefor1=para_edge.codefor1<<1;
		}
		para_edge.codefor1=para_edge.codefor1|1;
	}

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*20);
	struct para_getN tmp;
	tmp.kmerlen=kmer_len;
	tmp.FileName=inputFileName;
	tmp.InOut_label=0;
	tmp.isposition=0;
	tmp.isMiddle=1;
	getFileName(tmp);

	cout << "start constructing BplusTree..." <<endl;

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
	cout << total_hash_number << endl;
	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*para.kmer64Len);
	fseek(int_input,0,0);
	uint64_t x;
	x=fread(h1,sizeof(uint64_t),total_hash_number*para.kmer64Len,int_input);

	struct NodeBit * p_root;
	p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	struct nodeBit c_tmp;
	c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		for(uint64_t j=0;j<para.kmer64Len;j++)
		{
			c_tmp.hashValue[j]=h1[para.kmer64Len*i+j];
		}
		c_tmp.arrayID=i+1;
		Insert_Value_bit(&p_root,c_tmp,para);
	}
	fclose(int_input);

	for(uint64_t i=0;i<total_hash_number;i++)
	{
		uint64_t j;
		j=MappingHashValueToID_bit(p_root,h1+para.kmer64Len*i,para);
		if(j!=(i+1))
		{
			cout << "error: read kmer wrong!" << endl;
		}
	}
	cout << "Constructing BplusTree is over!" <<endl;

	char * p_tmp;
	char* dataset=seq;
	uint64_t start,end;
	start=0;
	end=0;

	struct NodeBit * p_root_unipath;
	p_root_unipath=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root_unipath);
	struct nodeBit c_tmp_unipath;
	c_tmp_unipath.p_posList_length=0;
	c_tmp_unipath.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);

//	struct NodeBit * p_root_edge;
//	p_root_edge=(struct NodeBit*)malloc(sizeof(struct NodeBit));
//	Node_initial_bit(p_root_edge);
//	struct nodeBit c_tmp_edge;
//	c_tmp_edge.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para_edge.kmer64Len);
//	c_tmp_edge.p_posList=(uint64_t *)malloc(sizeof(uint64_t)*1);
//	c_tmp_edge.p_posList_length=1;

	uint64_t label_pre=1;
	uint64_t R_len;
	uint64_t seq_k_current[4];
	cout << "start finding unipath ..." << endl;
	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	uint64_t template_64[64];
	template_64[0]=1;
	for(uint32_t i=0;i<63;i++)
	{
		template_64[i+1]=template_64[i]<<1;
	}

	uint64_t number_unipath=0;
	uint64_t total_br_num=0;
	uint64_t total_kmer_number;
	total_kmer_number=total_hash_number;
	for(uint64_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
	    cout << ref_i+1 << ":" << p_ref_path.NumberOfPathes << endl;
		ReadSeq(&p_tmp,&R_len,p_ref_path.pRefFilePath[ref_i]);

		uint64_t *p_candidate_current_in=NULL;
		uint64_t *p_candidate_current_out=NULL;

		if(tmps.pos_label==2||tmps.pos_label==3)
		{
			char* posfilename;
			posfilename=(char*)malloc(sizeof(char)*14);
			struct para_getN tmp_inName;
			tmp_inName.kmerlen=kmer_len;
			tmp_inName.FileName=posfilename;
			tmp_inName.InOut_label=0;
			tmp_inName.isMiddle=0;
			tmp_inName.isposition=ref_i+1;
			getFileName(tmp_inName);

			FILE* fp_pos_in;
			fp_pos_in=fopen(posfilename,"rb");
			if(fp_pos_in==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			uint64_t x;
			x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_in);

			tmp_inName.InOut_label=1;
			getFileName(tmp_inName);

			FILE* fp_pos_out;
			fp_pos_out=fopen(posfilename,"rb");
			if(fp_pos_out==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_out);
		}


		if(thread_num>1)
		{
			uint64_t **reslt;
			uint64_t **reslt_c;

			uint64_t **reslt_bran_pos;
			reslt=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
			reslt_c=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
			reslt_bran_pos=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

			uint64_t * reslt_size;
			reslt_size=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

			uint64_t * reslt_size_bran_pos;
			reslt_size_bran_pos=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

			struct para_unipath_bit* tmp;
			tmp=(struct para_unipath_bit*)malloc(sizeof(struct para_unipath_bit)*thread_num);

			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

			uint64_t unit_size=(R_len-kmer_len+1)/thread_num+1;

			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=j*unit_size;
				if((j+1)*unit_size>R_len-kmer_len)
				{
					tmp[j].end=R_len-kmer_len;
				}
				else
				{
					tmp[j].end=(j+1)*unit_size-1;
				}
				tmp[j].kmer_len=kmer_len;
				tmp[j].p_tmp=p_tmp;
				tmp[j].p_Branch_root=p_root;
				tmp[j].reslt=&(reslt[j]);
				tmp[j].reslt_c=&(reslt_c[j]);
				tmp[j].reslt_size=reslt_size+j;
				tmp[j].reslt_bran_pos=&(reslt_bran_pos[j]);
				tmp[j].reslt_size_bran_pos=reslt_size_bran_pos+j;
				tmp[j].thread_id=j;
				tmp[j].thread_num=thread_num;
				tmp[j].para=para;
				tmp[j].pos_label=pos_label;
				tmp[j].p_candidate_pos_in=p_candidate_current_in;
				tmp[j].p_candidate_pos_out=p_candidate_current_out;
				tmp[j].template_64=template_64;
				if(pthread_create(t+j, NULL, Find_unipath_para, (void*)(tmp+j))!=0)
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
					cal_hash_value_directly_256bit(p_tmp+reslt[j][k],c_tmp_unipath.hashValue,para);
					if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
					{
						c_tmp_unipath.arrayID=(((reslt[j][k]<<25)|((reslt_c[j][k]-reslt[j][k])<<8))|ref_i);
						Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
						total_kmer_number=total_kmer_number+(reslt_c[j][k]-reslt[j][k]);
						number_unipath++;
					}
				}
			}

			cout << "unipath finding is over!" << endl;


//			for(uint32_t j=0;j<thread_num;j++)
//			{
//				total_br_num+=reslt_size_bran_pos[j];
//				for(uint64_t k=0;k<reslt_size_bran_pos[j];k++)
//				{
//					cal_hash_value_directly_256bit(p_tmp+reslt_bran_pos[j][k],c_tmp_edge.hashValue,para_edge);
//					c_tmp_edge.p_posList[0]=(reslt_bran_pos[j][k]<<8)|(ref_i);
//					Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
//				}
//			}

			cout << "position list finding is over!" << endl;

			free(t);
			free(reslt_size);
			for(uint32_t i=0;i<thread_num;i++)
			{
				free(reslt_bran_pos[i]);
				free(reslt[i]);
				free(reslt_c[i]);
			}
			free(reslt_bran_pos);
			free(reslt_size_bran_pos);
			free(reslt);
			free(reslt_c);
			free(tmp);
		}
		else
		{
			cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
			for(uint64_t i=1;i<=R_len-kmer_len;i++)
			{
				cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

//				if(label_pre==1)
//				{
//					total_br_num++;
//					//save the bran-unip edge position
//					cal_hash_value_directly_256bit(p_tmp+i-1,c_tmp_edge.hashValue,para_edge);
//					c_tmp_edge.p_posList[0]=((i-1)<<8)|(ref_i);
//					Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
//				}

				if(pos_label==2||pos_label==3)
				{
					uint64_t tmp_pos_divide_in=i/64;
					uint64_t tmp_pos_mod_in=i%64;
					uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
					if(template_64[tmp_pos_mod_in]==tmp_template_in)
					{
						if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
						{
							//label_pre=0 is have no pre branch
							if(label_pre==0)
							{
								cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
								if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
								{
									c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|ref_i);
									Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
									//对unipath计数
//									MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
									total_kmer_number=total_kmer_number+(i-start);
									number_unipath++;
								}
							}
							start=i+1;
							label_pre=1;
						}
						else
						{
							label_pre=0;
						}
					}
					else
					{
						uint64_t tmp_pos=i+kmer_len-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
						if(template_64[tmp_pos_mod]==tmp_template)
						{
							if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
							{
								//label_pre=0 is have no pre branch
								if(label_pre==0)
								{
									cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
									if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
									{
										c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|ref_i);
										Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
										//对unipath计数
//										MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
										total_kmer_number=total_kmer_number+(i-start);
										number_unipath++;
									}
								}
								start=i+1;
								label_pre=1;
							}
							else
							{
								label_pre=0;
							}
						}
						else
						{
							label_pre=0;
						}
					}
				}
				else
				{
					if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
					{
						//label_pre=0 is have no pre branch
						if(label_pre==0)
						{
							cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
							if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
							{
								c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|ref_i);
								Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
								total_kmer_number=total_kmer_number+(i-start);
								number_unipath++;
							}
						}
						start=i+1;
						label_pre=1;
					}
					else
					{
						label_pre=0;
					}
				}
			}
		}
		if(pos_label==2||pos_label==3)
		{
			free(p_candidate_current_in);
			free(p_candidate_current_out);
		}
        if(rc_label==1)
        {
            cout << ref_i+1 << "RC:" << p_ref_path.NumberOfPathes << endl;
            start=0;
           	end=0;
           	label_pre=1;
            rc(&p_tmp,R_len);

            uint64_t *p_candidate_current_in=NULL;
            uint64_t *p_candidate_current_out=NULL;

            if(tmps.pos_label==2||tmps.pos_label==3)
            {
                char* posfilename;
                posfilename=(char*)malloc(sizeof(char)*14);
                struct para_getN tmp_inName;
                tmp_inName.kmerlen=kmer_len;
                tmp_inName.FileName=posfilename;
                tmp_inName.InOut_label=0;
                tmp_inName.isMiddle=0;
                tmp_inName.isposition=ref_i+1;
                getFileName(tmp_inName);
                posfilename[2]='r';
                posfilename[3]='c';

                FILE* fp_pos_in;
                fp_pos_in=fopen(posfilename,"rb");
                if(fp_pos_in==NULL)
                {
                    cout << "error: open candidate file failed!" << endl;
                }
                p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
                uint64_t x;
                x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
                if(x!=R_len/64+((R_len%64)?1:0))
                {
                    cout << "error: read current candidate failed!" << endl;
                }
                fclose(fp_pos_in);

                tmp_inName.InOut_label=1;
                getFileName(tmp_inName);
                posfilename[2]='r';
                posfilename[3]='c';

                FILE* fp_pos_out;
                fp_pos_out=fopen(posfilename,"rb");
                if(fp_pos_out==NULL)
                {
                    cout << "error: open candidate file failed!" << endl;
                }
                p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
                x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
                if(x!=R_len/64+((R_len%64)?1:0))
                {
                    cout << "error: read current candidate failed!" << endl;
                }
                fclose(fp_pos_out);
            }


            if(thread_num>1)
            {
                uint64_t **reslt;
                uint64_t **reslt_c;

                uint64_t **reslt_bran_pos;
                reslt=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
                reslt_c=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
                reslt_bran_pos=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

                uint64_t * reslt_size;
                reslt_size=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

                uint64_t * reslt_size_bran_pos;
                reslt_size_bran_pos=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

                struct para_unipath_bit* tmp;
                tmp=(struct para_unipath_bit*)malloc(sizeof(struct para_unipath_bit)*thread_num);

                pthread_t* t;
                t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

                uint64_t unit_size=(R_len-kmer_len+1)/thread_num+1;

                for(uint32_t j=0;j<thread_num;j++)
                {
                    tmp[j].start=j*unit_size;
                    if((j+1)*unit_size>R_len-kmer_len)
                    {
                        tmp[j].end=R_len-kmer_len;
                    }
                    else
                    {
                        tmp[j].end=(j+1)*unit_size-1;
                    }
                    tmp[j].kmer_len=kmer_len;
                    tmp[j].p_tmp=p_tmp;
                    tmp[j].p_Branch_root=p_root;
                    tmp[j].reslt=&(reslt[j]);
                    tmp[j].reslt_c=&(reslt_c[j]);
                    tmp[j].reslt_size=reslt_size+j;
                    tmp[j].reslt_bran_pos=&(reslt_bran_pos[j]);
                    tmp[j].reslt_size_bran_pos=reslt_size_bran_pos+j;
                    tmp[j].thread_id=j;
                    tmp[j].thread_num=thread_num;
                    tmp[j].para=para;
                    tmp[j].pos_label=pos_label;
                    tmp[j].p_candidate_pos_in=p_candidate_current_in;
                    tmp[j].p_candidate_pos_out=p_candidate_current_out;
                    tmp[j].template_64=template_64;
                    if(pthread_create(t+j, NULL, Find_unipath_para, (void*)(tmp+j))!=0)
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
                        cal_hash_value_directly_256bit(p_tmp+reslt[j][k],c_tmp_unipath.hashValue,para);
                        if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
                        {
                            c_tmp_unipath.arrayID=(((reslt[j][k]<<25)|((reslt_c[j][k]-reslt[j][k])<<8))|ref_i);
                            Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
                            total_kmer_number=total_kmer_number+(reslt_c[j][k]-reslt[j][k]);
                            number_unipath++;
                        }
                    }
                }

                cout << "unipath finding is over!" << endl;


    //			for(uint32_t j=0;j<thread_num;j++)
    //			{
    //				total_br_num+=reslt_size_bran_pos[j];
    //				for(uint64_t k=0;k<reslt_size_bran_pos[j];k++)
    //				{
    //					cal_hash_value_directly_256bit(p_tmp+reslt_bran_pos[j][k],c_tmp_edge.hashValue,para_edge);
    //					c_tmp_edge.p_posList[0]=(reslt_bran_pos[j][k]<<8)|(ref_i);
    //					Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
    //				}
    //			}

                cout << "position list finding is over!" << endl;

                free(t);
                free(reslt_size);
                for(uint32_t i=0;i<thread_num;i++)
                {
                    free(reslt_bran_pos[i]);
                    free(reslt[i]);
                    free(reslt_c[i]);
                }
                free(reslt_bran_pos);
                free(reslt_size_bran_pos);
                free(reslt);
                free(reslt_c);
                free(tmp);
            }
            else
            {
                cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
                for(uint64_t i=1;i<=R_len-kmer_len;i++)
                {
                    cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

    //				if(label_pre==1)
    //				{
    //					total_br_num++;
    //					//save the bran-unip edge position
    //					cal_hash_value_directly_256bit(p_tmp+i-1,c_tmp_edge.hashValue,para_edge);
    //					c_tmp_edge.p_posList[0]=((i-1)<<8)|(ref_i);
    //					Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
    //				}

                    if(pos_label==2||pos_label==3)
                    {
                        uint64_t tmp_pos_divide_in=i/64;
                        uint64_t tmp_pos_mod_in=i%64;
                        uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
                        if(template_64[tmp_pos_mod_in]==tmp_template_in)
                        {
                            if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
                            {
                                //label_pre=0 is have no pre branch
                                if(label_pre==0)
                                {
                                    cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
                                    if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
                                    {
                                        c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|(p_ref_path.NumberOfPathes+ref_i));
                                        Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
                                        //对unipath计数
    //									MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
                                        total_kmer_number=total_kmer_number+(i-start);
                                        number_unipath++;
                                    }
                                }
                                start=i+1;
                                label_pre=1;
                            }
                            else
                            {
                                label_pre=0;
                            }
                        }
                        else
                        {
                            uint64_t tmp_pos=i+kmer_len-1;
                            uint64_t tmp_pos_divide=tmp_pos/64;
                            uint64_t tmp_pos_mod=tmp_pos%64;
                            uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
                            if(template_64[tmp_pos_mod]==tmp_template)
                            {
                                if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
                                {
                                    //label_pre=0 is have no pre branch
                                    if(label_pre==0)
                                    {
                                        cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
                                        if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
                                        {
                                            c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|(p_ref_path.NumberOfPathes+ref_i));
                                            Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
                                            //对unipath计数
    //										MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
                                            total_kmer_number=total_kmer_number+(i-start);
                                            number_unipath++;
                                        }
                                    }
                                    start=i+1;
                                    label_pre=1;
                                }
                                else
                                {
                                    label_pre=0;
                                }
                            }
                            else
                            {
                                label_pre=0;
                            }
                        }
                    }
                    else
                    {
                        if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
                        {
                            //label_pre=0 is have no pre branch
                            if(label_pre==0)
                            {
                                cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
                                if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
                                {
                                    c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|(p_ref_path.NumberOfPathes+ref_i));
                                    Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
                                    total_kmer_number=total_kmer_number+(i-start);
                                    number_unipath++;
                                }
                            }
                            start=i+1;
                            label_pre=1;
                        }
                        else
                        {
                            label_pre=0;
                        }
                    }
                }
            }
            if(pos_label==2||pos_label==3)
            {
                free(p_candidate_current_in);
                free(p_candidate_current_out);
            }
        }
		free(p_tmp);
	}

	cout << "unipath and position list are ready!" << endl;

	struct NodeBit * p_leaf_hashtable;
	p_leaf_hashtable=Find_Minimal_Node_bit(p_root_unipath);

	char init[8]="unipath";
	char swap[4];
	strcpy(inputFileName,init);
	if(kmer_len<100)
	{
		inputFileName[7]='0';
		sprintf(swap,"%d",kmer_len);
		strcpy(inputFileName+8,swap);
	}
	else
	{
		sprintf(swap,"%d",kmer_len);
		strcpy(inputFileName+7,swap);
	}
	inputFileName[10]='\0';
	FILE* outUnipath;
	outUnipath=fopen(inputFileName,"wb+");

	FILE* outUnipath_c;
	outUnipath_c=fopen("unipath_counter","wb+");

	while(p_leaf_hashtable!=NULL)
	{
		for(uint32_t i=0;i<p_leaf_hashtable->Node_Size;i++)
		{
			fwrite(&(p_leaf_hashtable->data[i].arrayID),sizeof(uint64_t),1,outUnipath);
			fwrite(&(p_leaf_hashtable->data[i].p_posList_length),sizeof(uint64_t),1,outUnipath_c);
		}
		p_leaf_hashtable=p_leaf_hashtable->brother;
	}

//	struct NodeBit * p_pos_hashtable;
//	p_pos_hashtable=Find_Minimal_Node_bit(p_root_edge);

//	char init_pos[8]="poslist";
//	char swap_pos[3];
//	strcpy(inputFileName,init_pos);
//	sprintf(swap_pos,"%d",kmer_len);
//	strcpy(inputFileName+7,swap_pos);
//	inputFileName[9]='\0';
//	FILE* outPosList;
//	outPosList=fopen(inputFileName,"wb+");
//
//	uint64_t branch_emergy_times=0;
//	while(p_pos_hashtable!=NULL)
//	{
//		for(uint32_t i=0;i<p_pos_hashtable->Node_Size;i++)
//		{
//			branch_emergy_times=p_pos_hashtable->data[i].p_posList_length+branch_emergy_times;
//			fwrite(p_pos_hashtable->data[i].hashValue,sizeof(uint64_t),para_edge.kmer64Len,outPosList);
//			fwrite(&(p_pos_hashtable->data[i].p_posList_length),sizeof(uint64_t),1,outPosList);
//			fwrite(p_pos_hashtable->data[i].p_posList,sizeof(uint64_t),p_pos_hashtable->data[i].p_posList_length,outPosList);
//		}
//		p_pos_hashtable=p_pos_hashtable->brother;
//	}

//	FILE* out_total_kmer_number;
//	out_total_kmer_number=fopen("totalKmerNumber","wb");
//	fwrite(&total_kmer_number,sizeof(uint64_t),1,out_total_kmer_number);
//	fwrite(&branch_emergy_times,sizeof(uint64_t),1,out_total_kmer_number);
//	cout << total_kmer_number << endl;
//	cout << branch_emergy_times << endl;
//	fclose(out_total_kmer_number);

//	destory_tree_bit(p_root_edge,para_edge);
	destory_tree_bit(p_root_unipath,para);
	destory_tree_bit(p_root,para);
	free(inputFileName);
	fclose(outUnipath);
	fclose(outUnipath_c);
//	fclose(outPosList);
	cout << "the total times of branch kmer being in ref is:" <<total_br_num << endl;
	cout << "the number of unipath is:" << number_unipath << endl;
	cout << "finding unipath done!" << endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the finding time is: "<<span<< endl;
}

struct para_branch_node
{
    uint32_t kmer_len;
	char * p_tmp;
	struct bit256Hash * p_root;
	struct bit256KmerPara para;
	uint32_t pos_label;
	uint64_t *template_64;
	uint64_t *p_candidate_current_in;
	uint64_t *p_candidate_current_out;
	uint64_t start;
	uint64_t end;
	uint64_t *p_branch;
};
void *Find_branch_node_para(void *arg)
{
    struct para_branch_node* tmp=(struct para_branch_node*)arg;
    uint32_t kmer_len=tmp->kmer_len;
	char * p_tmp=tmp->p_tmp;
	struct bit256Hash* p_root=tmp->p_root;
	struct bit256KmerPara para=tmp->para;

	uint32_t pos_label=tmp->pos_label;
	uint64_t *template_64=tmp->template_64;
	uint64_t *p_candidate_current_in=tmp->p_candidate_current_in;
	uint64_t *p_candidate_current_out=tmp->p_candidate_current_out;

	uint64_t start=tmp->start;
	uint64_t end=tmp->end;
	uint64_t *p_branch=tmp->p_branch;

    uint64_t seq_k_current[4];

    uint64_t start1;
    if(start==0)
    {
		cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
		start1=start+1;
    }
    else
    {
		cal_hash_value_directly_256bit(p_tmp+start-1,seq_k_current,para);
		start1=start;
    }
    for(uint64_t i=start1;i<=end;i++)
    {
        cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);
        if(pos_label==2||pos_label==3)
        {
            uint64_t tmp_pos_divide_in=i/64;
            uint64_t tmp_pos_mod_in=i%64;
            uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
            if(template_64[tmp_pos_mod_in]==tmp_template_in)
            {
            	if(search256BitHashTable(p_root,seq_k_current,para)!=0)
                {
//                    uint64_t tmp_pos_divide=i/64;
//                    uint64_t tmp_pos_mod=i%64;

                    p_branch[tmp_pos_divide_in]=p_branch[tmp_pos_divide_in]|template_64[tmp_pos_mod_in];
                }
            }
            else
            {
                uint64_t tmp_pos=i+kmer_len-1;
                uint64_t tmp_pos_divide=tmp_pos/64;
                uint64_t tmp_pos_mod=tmp_pos%64;
                uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
                if(template_64[tmp_pos_mod]==tmp_template)
                {
                	if(search256BitHashTable(p_root,seq_k_current,para)!=0)
                    {
//                        uint64_t tmp_pos_divide=i/64;
//                        uint64_t tmp_pos_mod=i%64;

                        p_branch[tmp_pos_divide_in]=p_branch[tmp_pos_divide_in]|template_64[tmp_pos_mod_in];

                    }
                }
            }
        }
        else
        {
        	if(search256BitHashTable(p_root,seq_k_current,para)!=0)
            {
                uint64_t tmp_pos_divide=i/64;
                uint64_t tmp_pos_mod=i%64;

                p_branch[tmp_pos_divide]=p_branch[tmp_pos_divide]|template_64[tmp_pos_mod];

            }
        }
    }
}
struct new_para_unipath_bit
{
    uint32_t kmer_len;
	char * p_tmp;
	uint64_t R_len;
//	uint64_t ** p_unipath;
//	uint64_t * p_unipath_len;
//	uint64_t * p_unipath_Max_len;

	struct bit256Hash* p_root_unipath;
	struct NodeBit ** p_root_edge;
	uint32_t thread_id;
	uint32_t thread_num;
	struct bit256KmerPara para;
	struct bit256KmerPara para_edge;
	uint64_t *template_64;
	uint32_t ref_i;
	uint32_t rc_label;
	uint32_t ref_num;
	uint64_t *p_branch;
	uint32_t pl_label;
};
void *Find_unipath_para_new(void *arg)
{
    struct new_para_unipath_bit* tmp=(struct new_para_unipath_bit*)arg;
    uint32_t kmer_len=tmp->kmer_len;
	char * p_tmp=tmp->p_tmp;
	uint64_t R_len=tmp->R_len;
//	uint64_t * p_unipath=*tmp->p_unipath;
//	uint64_t p_unipath_len=*tmp->p_unipath_len;
//	uint64_t p_unipath_Max_len=*tmp->p_unipath_Max_len;
	struct bit256Hash* p_root_unipath=tmp->p_root_unipath;
	struct NodeBit * p_root_edge=*tmp->p_root_edge;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	struct bit256KmerPara para=tmp->para;
	struct bit256KmerPara para_edge=tmp->para_edge;
	uint64_t *template_64=tmp->template_64;
	uint32_t ref_i=tmp->ref_i;
	uint32_t rc_label=tmp->rc_label;
	uint32_t ref_num=tmp->ref_num;
	uint64_t *p_branch=tmp->p_branch;
	uint32_t pl_label=tmp->pl_label;

    uint64_t start,end;
    start=0;
    end=0;
    struct nodeBit c_tmp_unipath;
    c_tmp_unipath.p_posList_length=0;
    c_tmp_unipath.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
    uint64_t label_pre=1;
    uint64_t seq_k_current[4];

//	struct NodeBit * p_root_edge;
//	p_root_edge=(struct NodeBit*)malloc(sizeof(struct NodeBit));
//	Node_initial_bit(p_root_edge);
	struct nodeBit c_tmp_edge;
	c_tmp_edge.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para_edge.kmer64Len);
	c_tmp_edge.p_posList=(uint64_t *)malloc(sizeof(uint64_t)*1);
	c_tmp_edge.p_posList_length=1;

    cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
    for(uint64_t i=1;i<=R_len-kmer_len;i++)
    {
        cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);
        if(pl_label==1)
        {
			if(label_pre==1)
			{
//				total_br_num++;
				//save the bran-unip edge position
				cal_hash_value_directly_256bit(p_tmp+i-1,c_tmp_edge.hashValue,para_edge);
				c_tmp_edge.p_posList[0]=((i-1)<<8)|(ref_i);
				Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
			}
        }

        uint64_t tmp_pos_divide_in=i/64;
        uint64_t tmp_pos_mod_in=i%64;
        uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_branch[tmp_pos_divide_in];
        if(template_64[tmp_pos_mod_in]==tmp_template_in)
		{
			//label_pre=0 is have no pre branch
			if(label_pre==0)
			{
				cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
				uint32_t shift=bit256HashFunction(c_tmp_unipath.hashValue,para);
//					if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)

				if(shift%thread_num==thread_id)
				{
					if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
					{
						if(rc_label==1)
						{
							c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|(ref_i+ref_num));
						}
						else
						{
							c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|ref_i);
						}
//						Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
						insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
						//对unipath计数
//									MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
//                                total_kmer_number=total_kmer_number+(i-start);
//                                number_unipath++;
					}

				}
			}

			start=i+1;
			label_pre=1;
		}
		else
		{
			label_pre=0;
		}
    }

    *(tmp->p_root_edge)=p_root_edge;
}

//2)method 2. the parallel method in method 2 cost acceptable memory.but process all branching node one time. therefor it require at least memories for saving all branching k-mers on a single bplustree.
//3)method 3. the parallel method in method 2 cost acceptable memory.but process all branching node one time. method 3 can process large branching k-mer files.
struct para_calbranchvec
{
	char * file_branch_path;
	uint64_t * p_branch_vec;
	char * p_tmp;
	uint64_t R_len;
	uint32_t kmer_len;
	uint32_t rc_label;
	uint32_t pos_label;
	uint32_t thread_num;
	uint32_t task_size;
	uint32_t ref_i;
};
void calbranchvec(struct para_calbranchvec tmp)
{
	char * file_branch_path=tmp.file_branch_path;
	uint64_t * p_branch_vec=tmp.p_branch_vec;
	char * p_tmp=tmp.p_tmp;
	uint64_t R_len=tmp.R_len;
	uint32_t kmer_len=tmp.kmer_len;
	uint32_t rc_label=tmp.rc_label;
	uint32_t pos_label=tmp.pos_label;
	uint32_t thread_num=tmp.thread_num;
	uint32_t task_size=tmp.task_size;
	uint32_t ref_i=tmp.ref_i;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	if(para.remainer1to64==0)
	{
		para.codefor1=0;
	}
	else
	{
		para.codefor1=0;
		for(uint32_t i=0;i<para.remainer1to64-1;i++)
		{
			para.codefor1=para.codefor1|1;
			para.codefor1=para.codefor1<<1;
		}
		para.codefor1=para.codefor1|1;
	}

	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*20);
	struct para_getN tmp_para;
	tmp_para.kmerlen=kmer_len;
	tmp_para.FileName=inputFileName;
	tmp_para.InOut_label=0;
	tmp_para.isposition=0;
	tmp_para.isMiddle=1;
	getFileName(tmp_para);

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
	fseek(int_input,0,0);
	free(inputFileName);
	cout << total_hash_number << endl;

	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*task_size*para.kmer64Len);
    uint64_t template_64[64];
    template_64[0]=1;
    for(uint32_t i=0;i<63;i++)
    {
        template_64[i+1]=template_64[i]<<1;
    }
    uint64_t seq_k_current[4];

	uint32_t loop_num;
	loop_num=total_hash_number/task_size+((total_hash_number%task_size==0)?0:1);
	cout << total_hash_number << endl;
	cout << task_size << endl;
	cout << "loop_num:" << loop_num << endl;
	for(uint32_t i=0;i<loop_num;i++)
	{
		cout << i+1 << ":" << loop_num << endl;
		uint64_t x;
		uint32_t task_size_loop;
		if(i==loop_num-1)
		{
			task_size_loop=total_hash_number-i*task_size;
		}
		else
		{
			task_size_loop=task_size;
		}
		cout << task_size_loop << endl;
		x=fread(h1,sizeof(uint64_t),task_size_loop*para.kmer64Len,int_input);
		cout << x << endl;

//		struct NodeBit * p_root;
//		p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
//		Node_initial_bit(p_root);
//		struct nodeBit c_tmp;
//		c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
//		for(uint64_t i=0;i<task_size_loop;i++)
//		{
//			for(uint64_t j=0;j<para.kmer64Len;j++)
//			{
//				c_tmp.hashValue[j]=h1[para.kmer64Len*i+j];
//			}
//			c_tmp.arrayID=i+1;
//			Insert_Value_bit(&p_root,c_tmp,para);
//		}
//		for(uint64_t i=0;i<task_size_loop;i++)
//		{
//			uint64_t j;
//			j=MappingHashValueToID_bit(p_root,h1+para.kmer64Len*i,para);
//			if(j!=(i+1))
//			{
//				cout << "error: read kmer wrong!" << endl;
//			}
//		}

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
			for(uint64_t j=0;j<task_size_loop;j++)
			{
				insert256BitHashTable(p_root,h1+j*para.kmer64Len,j+1,para);
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
				p[i].Tasksize_current=task_size_loop;
				p[i].para1=para;
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
//		for(uint64_t j=0;j<task_size_loop;j++)
//		{
//			uint64_t tmp;
//			tmp=search256BitHashTable(p_root,h1+j*para.kmer64Len,para);
//			if(tmp!=j+1)
//			{
//				cout << j << endl;
//				cout << "error!" << endl;
//				char a;
//				cin >> a;
//			}
//		}
		gettimeofday(&tve_hs_c,NULL);
		double span_hs_c = tve_hs_c.tv_sec-tvs_hs_c.tv_sec + (tve_hs_c.tv_usec-tvs_hs_c.tv_usec)/1000000.0;
		cout << "check hash time is: "<<span_hs_c<<endl;
		cerr << "check hash kmer over!" << endl;
		if(thread_num==1)
		{
			if(pos_label==2 || pos_label==3)
			{
				uint64_t *p_candidate_current_in=NULL;
				uint64_t *p_candidate_current_out=NULL;

				char* posfilename;
				posfilename=(char*)malloc(sizeof(char)*14);
				struct para_getN tmp_inName;
				tmp_inName.kmerlen=kmer_len;
				tmp_inName.FileName=posfilename;
				tmp_inName.InOut_label=0;
				tmp_inName.isMiddle=0;
				tmp_inName.isposition=ref_i+1;
				getFileName(tmp_inName);
				if(rc_label==1)
				{
					posfilename[2]='r';
					posfilename[3]='c';
				}

				FILE* fp_pos_in;
				fp_pos_in=fopen(posfilename,"rb");
				if(fp_pos_in==NULL)
				{
					cout << "error: open candidate file failed!" << endl;
				}
				p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
				uint64_t x;
				x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
				if(x!=R_len/64+((R_len%64)?1:0))
				{
					cout << "error: read current candidate failed!" << endl;
				}
				fclose(fp_pos_in);

				tmp_inName.InOut_label=1;
				getFileName(tmp_inName);
				if(rc_label==1)
				{
					posfilename[2]='r';
					posfilename[3]='c';
				}

				FILE* fp_pos_out;
				fp_pos_out=fopen(posfilename,"rb");
				if(fp_pos_out==NULL)
				{
					cout << "error: open candidate file failed!" << endl;
				}
				p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
				x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
				if(x!=R_len/64+((R_len%64)?1:0))
				{
					cout << "error: read current candidate failed!" << endl;
				}
				fclose(fp_pos_out);

				cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
				for(uint64_t i=1;i<=R_len-kmer_len;i++)
				{
					cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

					uint64_t tmp_pos_divide_in=i/64;
					uint64_t tmp_pos_mod_in=i%64;
					uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
					if(template_64[tmp_pos_mod_in]==tmp_template_in)
					{
//						if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
						if(search256BitHashTable(p_root,seq_k_current,para)!=0)
						{
							//label_pre=0 is have no pre branch
							p_branch_vec[tmp_pos_divide_in]=p_branch_vec[tmp_pos_divide_in]|template_64[tmp_pos_mod_in];
						}
					}
					else
					{
						uint64_t tmp_pos=i+kmer_len-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
						if(template_64[tmp_pos_mod]==tmp_template)
						{
//							if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
							if(search256BitHashTable(p_root,seq_k_current,para)!=0)
							{
								//label_pre=0 is have no pre branch
								p_branch_vec[tmp_pos_divide_in]=p_branch_vec[tmp_pos_divide_in]|template_64[tmp_pos_mod_in];
							}
						}
					}
				}
				free(p_candidate_current_in);
				free(p_candidate_current_out);
			}
			else
			{
				cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
				for(uint64_t i=1;i<=R_len-kmer_len;i++)
				{
					cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

//					if(MappingHashValueToID_bit(p_root,seq_k_current,para)!=-1)
					if(search256BitHashTable(p_root,seq_k_current,para)!=0)
					{
						//label_pre=0 is have no pre branch
						uint64_t tmp_pos_divide_in=i/64;
						uint64_t tmp_pos_mod_in=i%64;
						p_branch_vec[tmp_pos_divide_in]=p_branch_vec[tmp_pos_divide_in]|template_64[tmp_pos_mod_in];
					}
				}
			}
		}
		else
		{
			uint64_t *p_candidate_current_in=NULL;
			uint64_t *p_candidate_current_out=NULL;
			if(pos_label==2 || pos_label==3)
			{
				char* posfilename;
				posfilename=(char*)malloc(sizeof(char)*14);
				struct para_getN tmp_inName;
				tmp_inName.kmerlen=kmer_len;
				tmp_inName.FileName=posfilename;
				tmp_inName.InOut_label=0;
				tmp_inName.isMiddle=0;
				tmp_inName.isposition=ref_i+1;
				getFileName(tmp_inName);
				if(rc_label==1)
				{
					posfilename[2]='r';
					posfilename[3]='c';
				}

				FILE* fp_pos_in;
				fp_pos_in=fopen(posfilename,"rb");
				if(fp_pos_in==NULL)
				{
					cout << "error: open candidate file failed!" << endl;
				}
				p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
				uint64_t x;
				x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
				if(x!=R_len/64+((R_len%64)?1:0))
				{
					cout << "error: read current candidate failed!" << endl;
				}
				fclose(fp_pos_in);

				tmp_inName.InOut_label=1;
				getFileName(tmp_inName);
				if(rc_label==1)
				{
					posfilename[2]='r';
					posfilename[3]='c';
				}

				FILE* fp_pos_out;
				fp_pos_out=fopen(posfilename,"rb");
				if(fp_pos_out==NULL)
				{
					cout << "error: open candidate file failed!" << endl;
				}
				p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
				x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
				if(x!=R_len/64+((R_len%64)?1:0))
				{
					cout << "error: read current candidate failed!" << endl;
				}
				fclose(fp_pos_out);
				free(posfilename);
			}

			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

			struct para_branch_node* para_unipath;
			para_unipath=(struct para_branch_node*)malloc(sizeof(struct para_branch_node)*thread_num);

			uint64_t unit_len=R_len/thread_num+1;
			uint64_t res=unit_len%64;
			if(res!=0)
			{
				unit_len=unit_len+(64-res);
			}

			for(uint32_t j=0;j<thread_num;j++)
			{
				para_unipath[j].kmer_len=kmer_len;
				para_unipath[j].para=para;
				para_unipath[j].pos_label=pos_label;
				para_unipath[j].p_candidate_current_in=p_candidate_current_in;
				para_unipath[j].p_candidate_current_out=p_candidate_current_out;
				para_unipath[j].p_root=p_root;
				para_unipath[j].p_tmp=p_tmp;
				para_unipath[j].template_64=template_64;
				para_unipath[j].start=j*unit_len;
				para_unipath[j].end=(j+1)*unit_len-1;
				if(para_unipath[j].end>R_len-kmer_len)
				{
					para_unipath[j].end=R_len-kmer_len;
				}
				para_unipath[j].p_branch=p_branch_vec;

				if(pthread_create(t+j, NULL, Find_branch_node_para, (void*)(para_unipath+j))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				pthread_join(t[j], NULL);
			}
			free(t);
			free(para_unipath);
			if(pos_label==2||pos_label==3)
			{
				free(p_candidate_current_in);
				free(p_candidate_current_out);
			}
		}
		freeHash256BitTable(p_root);
	}
	free(h1);
	fclose(int_input);
}
void FindUniPath_novel(struct para_findUnipath tmps)
{
	char* seq=tmps.seq;
	uint32_t kmer_len=tmps.original_kmer_length;
	uint32_t thread_num=tmps.thread_num;
	uint32_t pos_label=tmps.pos_label;
	uint32_t rc_label=tmps.rc_label;
	uint32_t memory_size=tmps.memory_size;
	uint32_t pl_label=tmps.pl_label;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	if(para.remainer1to64==0)
	{
		para.codefor1=0;
	}
	else
	{
		para.codefor1=0;
		for(uint32_t i=0;i<para.remainer1to64-1;i++)
		{
			para.codefor1=para.codefor1|1;
			para.codefor1=para.codefor1<<1;
		}
		para.codefor1=para.codefor1|1;
	}

	struct bit256KmerPara para_edge;
	para_edge.kmer1Len=kmer_len*2+2;
	para_edge.remainer1to64=para_edge.kmer1Len%64;
	para_edge.kmer64Len=para_edge.kmer1Len/64+(para_edge.remainer1to64?1:0);
	if(para_edge.remainer1to64==0)
	{
		para_edge.codefor1=0;
	}
	else
	{
		para_edge.codefor1=0;
		for(uint32_t i=0;i<para_edge.remainer1to64-1;i++)
		{
			para_edge.codefor1=para_edge.codefor1|1;
			para_edge.codefor1=para_edge.codefor1<<1;
		}
		para_edge.codefor1=para_edge.codefor1|1;
	}

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	char* branchfilename;
	branchfilename=(char*)malloc(sizeof(char)*20);
	struct para_getN tmp;
	tmp.kmerlen=kmer_len;
	tmp.FileName=branchfilename;
	tmp.InOut_label=0;
	tmp.isposition=0;
	tmp.isMiddle=1;
	getFileName(tmp);


	if(thread_num==1)
	{
		struct NodeBit * p_root_edge;
		struct nodeBit c_tmp_edge;

		if(pl_label==1)
		{
			p_root_edge=(struct NodeBit*)malloc(sizeof(struct NodeBit));
			Node_initial_bit(p_root_edge);
			c_tmp_edge.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para_edge.kmer64Len);
			c_tmp_edge.p_posList=(uint64_t *)malloc(sizeof(uint64_t)*1);
			c_tmp_edge.p_posList_length=1;
		}

		char * p_tmp;
		char* dataset=seq;
		uint64_t start,end;
		start=0;
		end=0;

		struct bit256Hash* p_root_unipath=NULL;
		p_root_unipath=initial256BitHashTable();
		if(p_root_unipath==NULL)
		{
			cerr<< "error!" << endl;
		}
//		insert256BitHashTable(p_root,h1+j*para.kmer64Len,j+1,para);
//		tmp=search256BitHashTable(p_root,h1+j*para.kmer64Len,para);

//		struct NodeBit * p_root_unipath;
//		p_root_unipath=(struct NodeBit*)malloc(sizeof(struct NodeBit));
//		Node_initial_bit(p_root_unipath);
		struct nodeBit c_tmp_unipath;
		c_tmp_unipath.p_posList_length=0;
		c_tmp_unipath.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);

		uint64_t label_pre=1;
		uint64_t R_len;
		uint64_t seq_k_current[4];
		cout << "start finding unipath ..." << endl;
		struct RefFilePath p_ref_path;
		getRefFilePathes(dataset, &p_ref_path);

		uint64_t template_64[64];
		template_64[0]=1;
		for(uint32_t i=0;i<63;i++)
		{
			template_64[i+1]=template_64[i]<<1;
		}
		uint32_t ref_num=p_ref_path.NumberOfPathes;

		for(uint64_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
		{
			cout << ref_i+1 << ":" << p_ref_path.NumberOfPathes << endl;
			ReadSeq(&p_tmp,&R_len,p_ref_path.pRefFilePath[ref_i]);

			uint64_t *p_branch;
			p_branch=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			memset(p_branch,0,sizeof(uint64_t)*(R_len/64+1));

			struct para_calbranchvec para_branch_vec;
			para_branch_vec.R_len=R_len;
			para_branch_vec.file_branch_path=branchfilename;
			para_branch_vec.kmer_len=kmer_len;
			para_branch_vec.p_branch_vec=p_branch;
			para_branch_vec.p_tmp=p_tmp;
			para_branch_vec.pos_label=pos_label;
			para_branch_vec.rc_label=0;
			para_branch_vec.ref_i=ref_i;
			para_branch_vec.thread_num=thread_num;
			para_branch_vec.task_size=80000000*memory_size;

			calbranchvec(para_branch_vec);

		    cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
		    for(uint64_t i=1;i<=R_len-kmer_len;i++)
		    {
		        cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

		        if(pl_label==1)
		        {
					if(label_pre==1)
					{
//						total_br_num++;
						//save the bran-unip edge position
						cal_hash_value_directly_256bit(p_tmp+i-1,c_tmp_edge.hashValue,para_edge);
						c_tmp_edge.p_posList[0]=((i-1)<<8)|(ref_i);
						Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
					}
		        }


		        uint64_t tmp_pos_divide_in=i/64;
		        uint64_t tmp_pos_mod_in=i%64;
		        uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_branch[tmp_pos_divide_in];
		        if(template_64[tmp_pos_mod_in]==tmp_template_in)
				{
					//label_pre=0 is have no pre branch
					if(label_pre==0)
					{
						cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
//						if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
						if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
						{
							c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|ref_i);
//							Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
							insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
							//对unipath计数
	//									MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
	//                                total_kmer_number=total_kmer_number+(i-start);
	//                                number_unipath++;
						}
					}
					start=i+1;
					label_pre=1;
				}
				else
				{
					label_pre=0;
				}
		    }


			if(rc_label==1)
			{
				cout << ref_i+1 << "RC:" << p_ref_path.NumberOfPathes << endl;
				rc(&p_tmp,R_len);
				memset(p_branch,0,sizeof(uint64_t)*(R_len/64+1));
				start=0;
				end=0;
				label_pre=1;

				struct para_calbranchvec para_branch_vec;
				para_branch_vec.R_len=R_len;
				para_branch_vec.file_branch_path=branchfilename;
				para_branch_vec.kmer_len=kmer_len;
				para_branch_vec.p_branch_vec=p_branch;
				para_branch_vec.p_tmp=p_tmp;
				para_branch_vec.pos_label=pos_label;
				para_branch_vec.rc_label=1;
				para_branch_vec.ref_i=ref_i;
				para_branch_vec.thread_num=thread_num;
				para_branch_vec.task_size=80000000*memory_size;

				calbranchvec(para_branch_vec);

				cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
				for(uint64_t i=1;i<=R_len-kmer_len;i++)
				{
					cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

					if(pl_label==1)
					{
						if(label_pre==1)
						{
//							total_br_num++;
							//save the bran-unip edge position
							cal_hash_value_directly_256bit(p_tmp+i-1,c_tmp_edge.hashValue,para_edge);
							c_tmp_edge.p_posList[0]=((i-1)<<8)|(ref_i);
							Insert_pos_256bitKmer(&p_root_edge,c_tmp_edge,para_edge);
						}
					}


					uint64_t tmp_pos_divide_in=i/64;
					uint64_t tmp_pos_mod_in=i%64;
					uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_branch[tmp_pos_divide_in];
					if(template_64[tmp_pos_mod_in]==tmp_template_in)
					{
						//label_pre=0 is have no pre branch
						if(label_pre==0)
						{
							cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
//							if(MappingHashValueToID_bit(p_root_unipath,c_tmp_unipath.hashValue,para)==-1)
							if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
							{
								c_tmp_unipath.arrayID=(((start<<25)|((i-start)<<8))|(ref_i+ref_num));
//								Insert_Value_bit(&p_root_unipath,c_tmp_unipath,para);
								insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
								//对unipath计数
		//									MappingHashValueToID_bit_unipath_count(p_root_unipath,c_tmp_unipath.hashValue,para);
		//                                total_kmer_number=total_kmer_number+(i-start);
		//                                number_unipath++;
							}
						}
						start=i+1;
						label_pre=1;
					}
					else
					{
						label_pre=0;
					}
				}

			}
			free(p_tmp);
			free(p_branch);
		}

		cout << "unipath and position list are ready!" << endl;

//		struct NodeBit * p_leaf_hashtable;
//		p_leaf_hashtable=Find_Minimal_Node_bit(p_root_unipath);

		char *inputFileName=(char*)malloc(sizeof(char)*20);
		char init[8]="unipath";
		char swap[4];
		strcpy(inputFileName,init);
		if(kmer_len<100)
		{
			inputFileName[7]='0';
			sprintf(swap,"%d",kmer_len);
			strcpy(inputFileName+8,swap);
		}
		else
		{
			sprintf(swap,"%d",kmer_len);
			strcpy(inputFileName+7,swap);
		}
		inputFileName[10]='\0';
		FILE* outUnipath;
		outUnipath=fopen(inputFileName,"wb+");

		FILE* outUnipath_c;
		outUnipath_c=fopen("unipath_counter","wb+");

		for(uint32_t i=0;i<HashSize;i++)
		{
			if(p_root_unipath->len[i]!=0)
			{
				fwrite(p_root_unipath->p_pos[i],sizeof(uint64_t),p_root_unipath->len[i],outUnipath);
			}
		}
//		while(p_leaf_hashtable!=NULL)
//		{
//			for(uint32_t i=0;i<p_leaf_hashtable->Node_Size;i++)
//			{
//				fwrite(&(p_leaf_hashtable->data[i].arrayID),sizeof(uint64_t),1,outUnipath);
//				fwrite(&(p_leaf_hashtable->data[i].p_posList_length),sizeof(uint64_t),1,outUnipath_c);
//			}
//			p_leaf_hashtable=p_leaf_hashtable->brother;
//		}

		if(pl_label==1)
		{
			struct NodeBit * p_pos_hashtable;
			p_pos_hashtable=Find_Minimal_Node_bit(p_root_edge);

			char init_pos[8]="poslist";
			char swap_pos[3];
			strcpy(inputFileName,init_pos);
			sprintf(swap_pos,"%d",kmer_len);
			strcpy(inputFileName+7,swap_pos);
			inputFileName[9]='\0';
			FILE* outPosList;
			outPosList=fopen(inputFileName,"wb+");

			uint64_t branch_emergy_times=0;
			while(p_pos_hashtable!=NULL)
			{
				for(uint32_t i=0;i<p_pos_hashtable->Node_Size;i++)
				{
					branch_emergy_times=p_pos_hashtable->data[i].p_posList_length+branch_emergy_times;
					fwrite(p_pos_hashtable->data[i].hashValue,sizeof(uint64_t),para_edge.kmer64Len,outPosList);
					fwrite(&(p_pos_hashtable->data[i].p_posList_length),sizeof(uint64_t),1,outPosList);
					fwrite(p_pos_hashtable->data[i].p_posList,sizeof(uint64_t),p_pos_hashtable->data[i].p_posList_length,outPosList);
				}
				p_pos_hashtable=p_pos_hashtable->brother;
			}

			FILE* out_total_kmer_number;
			out_total_kmer_number=fopen("totalKmerNumber","wb");
//			fwrite(&total_kmer_number,sizeof(uint64_t),1,out_total_kmer_number);
			fwrite(&branch_emergy_times,sizeof(uint64_t),1,out_total_kmer_number);
//			cout << total_kmer_number << endl;
			cout << branch_emergy_times << endl;
			fclose(out_total_kmer_number);
			destory_tree_bit(p_root_edge,para_edge);
		}

//		destory_tree_bit(p_root_unipath,para);
		freeHash256BitTable(p_root_unipath);
		free(inputFileName);
		fclose(outUnipath);
		fclose(outUnipath_c);
	//	fclose(outPosList);
//		cout << "the total times of branch kmer being in ref is:" <<total_br_num << endl;
//		cout << "the number of unipath is:" << number_unipath << endl;
		cout << "finding unipath done!" << endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"the finding time is: "<<span<< endl;
	}
	else
	{
		char * p_tmp;
		char* dataset=seq;

		struct bit256Hash* p_root_unipath=NULL;
		p_root_unipath=initial256BitHashTable();
		if(p_root_unipath==NULL)
		{
			cerr<< "error!" << endl;
		}

		struct NodeBit ** p_root_edge;
		if(pl_label==1)
		{
			p_root_edge=(struct NodeBit**)malloc(sizeof(struct NodeBit)*thread_num);
			for(uint32_t i=0;i<thread_num;i++)
			{
				p_root_edge[i]=(struct NodeBit*)malloc(sizeof(struct NodeBit));
				Node_initial_bit(p_root_edge[i]);
			}
		}

		uint64_t R_len;
		cout << "start finding unipath ..." << endl;
		struct RefFilePath p_ref_path;
		getRefFilePathes(dataset, &p_ref_path);

		uint64_t template_64[64];
		template_64[0]=1;
		for(uint32_t i=0;i<63;i++)
		{
			template_64[i+1]=template_64[i]<<1;
		}

		uint64_t *p_branch=NULL;
		for(uint64_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
		{
			cout << ref_i+1 << ":" << p_ref_path.NumberOfPathes << endl;
			ReadSeq(&p_tmp,&R_len,p_ref_path.pRefFilePath[ref_i]);

			p_branch=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			memset(p_branch,0,sizeof(uint64_t)*(R_len/64+1));

			struct para_calbranchvec para_branch_vec;
			para_branch_vec.R_len=R_len;
			para_branch_vec.file_branch_path=branchfilename;
			para_branch_vec.kmer_len=kmer_len;
			para_branch_vec.p_branch_vec=p_branch;
			para_branch_vec.p_tmp=p_tmp;
			para_branch_vec.pos_label=pos_label;
			para_branch_vec.rc_label=0;
			para_branch_vec.ref_i=ref_i;
			para_branch_vec.thread_num=thread_num;
			para_branch_vec.task_size=80000000*memory_size;

			calbranchvec(para_branch_vec);

			cout << "done!" << endl;

			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

			struct new_para_unipath_bit* para_unipath1;
			para_unipath1=(struct new_para_unipath_bit*)malloc(sizeof(struct new_para_unipath_bit)*thread_num);
			for(uint32_t j=0;j<thread_num;j++)
			{
				para_unipath1[j].kmer_len=kmer_len;
				para_unipath1[j].para=para;
				para_unipath1[j].p_tmp=p_tmp;
				para_unipath1[j].template_64=template_64;
				para_unipath1[j].p_branch=p_branch;
				para_unipath1[j].p_root_unipath=p_root_unipath;
				para_unipath1[j].p_root_edge=&(p_root_edge[j]);
				para_unipath1[j].rc_label=0;
				para_unipath1[j].ref_i=ref_i;
				para_unipath1[j].ref_num=p_ref_path.NumberOfPathes;
				para_unipath1[j].thread_id=j;
				para_unipath1[j].thread_num=thread_num;
				para_unipath1[j].R_len=R_len;
				para_unipath1[j].pl_label=pl_label;
				para_unipath1[j].para_edge=para_edge;

				if(pthread_create(t+j, NULL, Find_unipath_para_new, (void*)(para_unipath1+j))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				pthread_join(t[j], NULL);
			}

			if(t!=NULL){free(t);}
			if(para_unipath1!=NULL){free(para_unipath1);}
			if(rc_label==1)
			{
				cout << ref_i+1 << "RC:" << p_ref_path.NumberOfPathes << endl;
				rc(&p_tmp,R_len);
				memset(p_branch,0,sizeof(uint64_t)*(R_len/64+1));

				struct para_calbranchvec para_branch_vec;
				para_branch_vec.R_len=R_len;
				para_branch_vec.file_branch_path=branchfilename;
				para_branch_vec.kmer_len=kmer_len;
				para_branch_vec.p_branch_vec=p_branch;
				para_branch_vec.p_tmp=p_tmp;
				para_branch_vec.pos_label=pos_label;
				para_branch_vec.rc_label=1;
				para_branch_vec.ref_i=ref_i;
				para_branch_vec.thread_num=thread_num;
				para_branch_vec.task_size=80000000*memory_size;

				calbranchvec(para_branch_vec);

				pthread_t* t;
				t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

				struct new_para_unipath_bit* para_unipath1;
				para_unipath1=(struct new_para_unipath_bit*)malloc(sizeof(struct new_para_unipath_bit)*thread_num);
				for(uint32_t j=0;j<thread_num;j++)
				{
					para_unipath1[j].kmer_len=kmer_len;
					para_unipath1[j].para=para;
					para_unipath1[j].p_tmp=p_tmp;
					para_unipath1[j].template_64=template_64;
					para_unipath1[j].p_branch=p_branch;
					para_unipath1[j].p_root_unipath=p_root_unipath;
					para_unipath1[j].p_root_edge=&(p_root_edge[j]);
					para_unipath1[j].rc_label=1;
					para_unipath1[j].ref_i=ref_i;
					para_unipath1[j].ref_num=p_ref_path.NumberOfPathes;
					para_unipath1[j].thread_id=j;
					para_unipath1[j].thread_num=thread_num;
					para_unipath1[j].R_len=R_len;
					para_unipath1[j].pl_label=pl_label;
					para_unipath1[j].para_edge=para_edge;

					if(pthread_create(t+j, NULL, Find_unipath_para_new, (void*)(para_unipath1+j))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t j=0;j<thread_num;j++)
				{
					pthread_join(t[j], NULL);
				}

				if(t!=NULL){free(t);}
				if(para_unipath1!=NULL){free(para_unipath1);}
			}
			if(p_branch!=NULL){free(p_branch);}
			if(p_tmp!=NULL){free(p_tmp);}
		}
		cout << "unipath and position list are ready!" << endl;

		char *inputFileName=(char*)malloc(sizeof(char)*20);
		char init[8]="unipath";
		char swap[4];
		strcpy(inputFileName,init);
		if(kmer_len<100)
		{
			inputFileName[7]='0';
			sprintf(swap,"%d",kmer_len);
			strcpy(inputFileName+8,swap);
		}
		else
		{
			sprintf(swap,"%d",kmer_len);
			strcpy(inputFileName+7,swap);
		}
		inputFileName[10]='\0';
		FILE* outUnipath;
		outUnipath=fopen(inputFileName,"wb+");

		FILE* outUnipath_c;
		outUnipath_c=fopen("unipath_counter","wb+");

		for(uint32_t i=0;i<HashSize;i++)
		{
			if(p_root_unipath->len[i]!=0)
			{
				fwrite(p_root_unipath->p_pos[i],sizeof(uint64_t),p_root_unipath->len[i],outUnipath);
			}
		}

		if(pl_label==1)
		{
			for(uint32_t i=0;i<thread_num;i++)
			{
				struct NodeBit * p_pos_hashtable;
				p_pos_hashtable=Find_Minimal_Node_bit(p_root_edge[i]);

				char init_pos[8]="poslist";
				char swap_pos[3];
				strcpy(inputFileName,init_pos);
				sprintf(swap_pos,"%d",kmer_len);
				strcpy(inputFileName+7,swap_pos);
				inputFileName[9]='\0';
				FILE* outPosList;
				outPosList=fopen(inputFileName,"wb+");

				uint64_t branch_emergy_times=0;
				while(p_pos_hashtable!=NULL)
				{
					for(uint32_t i=0;i<p_pos_hashtable->Node_Size;i++)
					{
						branch_emergy_times=p_pos_hashtable->data[i].p_posList_length+branch_emergy_times;
						fwrite(p_pos_hashtable->data[i].hashValue,sizeof(uint64_t),para_edge.kmer64Len,outPosList);
						fwrite(&(p_pos_hashtable->data[i].p_posList_length),sizeof(uint64_t),1,outPosList);
						fwrite(p_pos_hashtable->data[i].p_posList,sizeof(uint64_t),p_pos_hashtable->data[i].p_posList_length,outPosList);
					}
					p_pos_hashtable=p_pos_hashtable->brother;
				}

				FILE* out_total_kmer_number;
				out_total_kmer_number=fopen("totalKmerNumber","wb");
				fwrite(&branch_emergy_times,sizeof(uint64_t),1,out_total_kmer_number);
				cout << branch_emergy_times << endl;
				fclose(out_total_kmer_number);

				destory_tree_bit(p_root_edge[i],para_edge);
			}
			free(p_root_edge);
		}
		freeHash256BitTable(p_root_unipath);
		free(inputFileName);
		fclose(outUnipath);
		fclose(outUnipath_c);
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"the finding time is: "<<span<< endl;
	}
}

//3) method
struct para_Find_Branch
{
    uint32_t kmer_len;
	char * p_tmp;
	struct bit256Hash * p_root;
	struct bit256KmerPara para;
	uint32_t pos_label;
	uint32_t thread_id;
	uint32_t rc_label;
	uint64_t *template_64;
	uint64_t *p_candidate_current_in;
	uint64_t *p_candidate_current_out;
	uint64_t *p_candidate_current;
	uint64_t start;
	uint64_t end;
	struct bit256Hash ** p_unipath;
	uint32_t candidate_label;
	uint32_t ref_i;
	uint32_t ref_num;
};
void *Find_Branch(void *arg)
{
    struct para_Find_Branch* tmp=(struct para_Find_Branch*)arg;
    uint32_t kmer_len=tmp->kmer_len;
	char * p_tmp=tmp->p_tmp;
	struct bit256Hash* p_root=tmp->p_root;
	struct bit256KmerPara para=tmp->para;

	uint32_t pos_label=tmp->pos_label;
	uint32_t rc_label=tmp->rc_label;
	uint64_t *template_64=tmp->template_64;
	uint64_t *p_candidate_current_in=tmp->p_candidate_current_in;
	uint64_t *p_candidate_current_out=tmp->p_candidate_current_out;
	uint64_t *p_candidate_current=tmp->p_candidate_current;
	uint32_t candidate_label=tmp->candidate_label;
	uint32_t ref_i=tmp->ref_i;
	uint32_t ref_num=tmp->ref_num;
	uint32_t thread_id=tmp->thread_id;

	struct bit256Hash* p_root_unipath;
	p_root_unipath=initial256BitHashTable();

	uint64_t start1=tmp->start;
	uint64_t end1=tmp->end;

    uint64_t seq_k_current[4];
    uint64_t start;
	uint64_t end;
	uint32_t label_pre;
	start=start1;
	end=start1;
	label_pre=1;
	struct nodeBit c_tmp_unipath;
	c_tmp_unipath.p_posList_length=0;
	c_tmp_unipath.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);

    uint64_t start2;
    if(start1==0)
    {
		cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
		start2=start1+1;
    }
    else
    {
		cal_hash_value_directly_256bit(p_tmp+start1-1,seq_k_current,para);
		start2=start1;
    }
    for(uint64_t i=start2;i<=end1;i++)
    {
        cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);
        if(pos_label==2||pos_label==3)
        {
        	if(candidate_label==0)
        	{
        		uint64_t tmp_pos_divide_in=i/64;
				uint64_t tmp_pos_mod_in=i%64;
				uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
				if(template_64[tmp_pos_mod_in]==tmp_template_in)
				{
					if(search256BitHashTable(p_root,seq_k_current,para)!=0)
					{
						if(label_pre==0)
						{
							cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
							if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
							{
								end=i;
								if(rc_label==1)
								{
									c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
								}
								else
								{
									c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
								}
								insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
							}
						}

						start=i+1;
						label_pre=1;
					}
					else
					{
						label_pre=0;
					}
				}
				else
				{
					uint64_t tmp_pos=i+kmer_len-1;
					uint64_t tmp_pos_divide=tmp_pos/64;
					uint64_t tmp_pos_mod=tmp_pos%64;
					uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
					if(template_64[tmp_pos_mod]==tmp_template)
					{
						if(search256BitHashTable(p_root,seq_k_current,para)!=0)
						{
							if(label_pre==0)
							{
								cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
								if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
								{
									end=i;
									if(rc_label==1)
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
									}
									else
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
									}
									insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
								}
							}

							start=i+1;
							label_pre=1;
						}
						else
						{
							label_pre=0;
						}
					}
					else
					{
						label_pre=0;
					}
				}
        	}
        	else
        	{
				uint64_t tmp_pos_divide_in=i/64;
				uint64_t tmp_pos_mod_in=i%64;
				uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current[tmp_pos_divide_in];
				if(template_64[tmp_pos_mod_in]==tmp_template_in)
				{
					if(search256BitHashTable(p_root,seq_k_current,para)!=0)
					{
						if(label_pre==0)
						{
							cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
							if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
							{
								end=i;
								if(rc_label==1)
								{
									c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
								}
								else
								{
									c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
								}
								insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
							}
						}

						start=i+1;
						label_pre=1;
					}
					else
					{
						label_pre=0;
					}
				}
				else
				{
					label_pre=0;
				}
        	}

        }
        else
        {
        	if(search256BitHashTable(p_root,seq_k_current,para)!=0)
            {
        		if(label_pre==0)
				{
					cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
					if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
					{
						end=i;
						if(rc_label==1)
						{
							c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
						}
						else
						{
							c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
						}
						insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
					}
				}

				start=i+1;
				label_pre=1;
            }
        	else
        	{
        		label_pre=0;
        	}
        }
    }

    tmp->p_unipath[thread_id]=p_root_unipath;
}
struct para_CalBranchVec
{
	char * p_tmp;
	uint64_t R_len;
	uint32_t kmer_len;
	uint32_t rc_label;
	uint32_t pos_label;
	uint32_t thread_num;
	uint32_t ref_i;
	uint32_t ref_num;
	uint32_t candidate_label;
	struct bit256Hash* p_root;
	struct bit256Hash** p_root_unipath;
};
void CalBranchVec(struct para_CalBranchVec tmp)
{
	char * p_tmp=tmp.p_tmp;
	uint64_t R_len=tmp.R_len;
	uint32_t kmer_len=tmp.kmer_len;
	uint32_t rc_label=tmp.rc_label;
	uint32_t pos_label=tmp.pos_label;
	uint32_t thread_num=tmp.thread_num;
	uint32_t ref_i=tmp.ref_i;
	uint32_t ref_num=tmp.ref_num;
	struct bit256Hash* p_root=tmp.p_root;
	uint32_t candidate_label=tmp.candidate_label;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	if(para.remainer1to64==0)
	{
		para.codefor1=0;
	}
	else
	{
		para.codefor1=0;
		for(uint32_t i=0;i<para.remainer1to64-1;i++)
		{
			para.codefor1=para.codefor1|1;
			para.codefor1=para.codefor1<<1;
		}
		para.codefor1=para.codefor1|1;
	}

	struct bit256Hash* p_root_unipath;
	struct bit256Hash** p_root_unipath_arr;

	uint64_t template_64[64];
	template_64[0]=1;
	for(uint32_t i=0;i<63;i++)
	{
		template_64[i+1]=template_64[i]<<1;
	}
	uint64_t seq_k_current[4];

	if(thread_num==1)
	{
		p_root_unipath=initial256BitHashTable();

		uint64_t start;
		uint64_t end;
		uint32_t label_pre;
		start=0;
		end=0;
		label_pre=1;
		struct nodeBit c_tmp_unipath;
		c_tmp_unipath.p_posList_length=0;
		c_tmp_unipath.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);

		if(pos_label==2 || pos_label==3)
		{
			uint64_t *p_candidate_current_in=NULL;
			uint64_t *p_candidate_current_out=NULL;

			char* posfilename;
			posfilename=(char*)malloc(sizeof(char)*14);
			struct para_getN tmp_inName;
			tmp_inName.kmerlen=kmer_len;
			tmp_inName.FileName=posfilename;
			tmp_inName.InOut_label=0;
			tmp_inName.isMiddle=0;
			tmp_inName.isposition=ref_i+1;
			getFileName(tmp_inName);
			if(rc_label==1)
			{
				posfilename[2]='r';
				posfilename[3]='c';
			}

			FILE* fp_pos_in;
			fp_pos_in=fopen(posfilename,"rb");
			if(fp_pos_in==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			uint64_t x;
			x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_in);

			tmp_inName.InOut_label=1;
			getFileName(tmp_inName);
			if(rc_label==1)
			{
				posfilename[2]='r';
				posfilename[3]='c';
			}

			FILE* fp_pos_out;
			fp_pos_out=fopen(posfilename,"rb");
			if(fp_pos_out==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_out);

			if(candidate_label==0)
			{
				cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
				for(uint64_t i=1;i<=R_len-kmer_len;i++)
				{
					cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

					uint64_t tmp_pos_divide_in=i/64;
					uint64_t tmp_pos_mod_in=i%64;
					uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current_in[tmp_pos_divide_in];
					if(template_64[tmp_pos_mod_in]==tmp_template_in)
					{
						if(search256BitHashTable(p_root,seq_k_current,para)!=0)
						{
							if(label_pre==0)
							{
								cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
								if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
								{
									end=i;
									if(rc_label==1)
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
									}
									else
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
									}
									insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
								}
							}

							start=i+1;
							label_pre=1;
						}
						else
						{
							label_pre=0;
						}
					}
					else
					{
						uint64_t tmp_pos=i+kmer_len-1;
						uint64_t tmp_pos_divide=tmp_pos/64;
						uint64_t tmp_pos_mod=tmp_pos%64;
						uint64_t tmp_template=template_64[tmp_pos_mod]&p_candidate_current_out[tmp_pos_divide];
						if(template_64[tmp_pos_mod]==tmp_template)
						{
							if(search256BitHashTable(p_root,seq_k_current,para)!=0)
							{
								if(label_pre==0)
								{
									cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
									if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
									{
										end=i;
										if(rc_label==1)
										{
											c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
										}
										else
										{
											c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
										}
										insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
									}
								}

								start=i+1;
								label_pre=1;
							}
							else
							{
								label_pre=0;
							}
						}
						else
						{
							label_pre=0;
						}
					}
				}
				free(p_candidate_current_in);
				free(p_candidate_current_out);

			}
			else
			{
				uint64_t *p_candidate_current=NULL;
				p_candidate_current=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
				uint32_t l_shift_size=kmer_len-1;
				uint32_t r_shift_size=64-l_shift_size;

				for(uint32_t i=0;i<R_len/64;i++)
				{
					p_candidate_current[i]=p_candidate_current_in[i]|((p_candidate_current_out[i]<<l_shift_size)|(p_candidate_current_out[i+1]>>r_shift_size));
				}
				p_candidate_current[R_len/64]=p_candidate_current_in[R_len/64]|(p_candidate_current_out[R_len/64]<<l_shift_size);
				free(p_candidate_current_in);
				free(p_candidate_current_out);

				cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
				for(uint64_t i=1;i<=R_len-kmer_len;i++)
				{
					cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

					uint64_t tmp_pos_divide_in=i/64;
					uint64_t tmp_pos_mod_in=i%64;
					uint64_t tmp_template_in=template_64[tmp_pos_mod_in]&p_candidate_current[tmp_pos_divide_in];
					if(template_64[tmp_pos_mod_in]==tmp_template_in)
					{
						if(search256BitHashTable(p_root,seq_k_current,para)!=0)
						{
							if(label_pre==0)
							{
								cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
								if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
								{
									end=i;
									if(rc_label==1)
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
									}
									else
									{
										c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
									}
									insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
								}
							}

							start=i+1;
							label_pre=1;
						}
						else
						{
							label_pre=0;
						}
					}
					else
					{
						label_pre=0;
					}
				}
				free(p_candidate_current);
			}
		}
		else
		{
			cal_hash_value_directly_256bit(p_tmp,seq_k_current,para);
			for(uint64_t i=1;i<=R_len-kmer_len;i++)
			{
				cal_hash_value_indirectly_256bit(p_tmp+i,seq_k_current,seq_k_current,para);

				if(search256BitHashTable(p_root,seq_k_current,para)!=0)
				{
					if(label_pre==0)
					{
						cal_hash_value_directly_256bit(p_tmp+start,c_tmp_unipath.hashValue,para);
						if(search256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,para)==0)
						{
							end=i;
							if(rc_label==1)
							{
								c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|(ref_i+ref_num));
							}
							else
							{
								c_tmp_unipath.arrayID=(((start<<25)|((end-start)<<8))|ref_i);
							}
							insert256BitHashTable(p_root_unipath,c_tmp_unipath.hashValue,c_tmp_unipath.arrayID,para);
						}
					}

					start=i+1;
					label_pre=1;
				}
				else
				{
					label_pre=0;
				}
			}
		}
	}
	else
	{
		p_root_unipath_arr=(struct bit256Hash**)malloc(sizeof(struct bit256Hash*)*thread_num);

		uint64_t *p_candidate_current_in=NULL;
		uint64_t *p_candidate_current_out=NULL;
		if(pos_label==2 || pos_label==3)
		{
			char* posfilename;
			posfilename=(char*)malloc(sizeof(char)*14);
			struct para_getN tmp_inName;
			tmp_inName.kmerlen=kmer_len;
			tmp_inName.FileName=posfilename;
			tmp_inName.InOut_label=0;
			tmp_inName.isMiddle=0;
			tmp_inName.isposition=ref_i+1;
			getFileName(tmp_inName);
			if(rc_label==1)
			{
				posfilename[2]='r';
				posfilename[3]='c';
			}

			FILE* fp_pos_in;
			fp_pos_in=fopen(posfilename,"rb");
			if(fp_pos_in==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_in=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			uint64_t x;
			x=fread(p_candidate_current_in,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_in);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_in);

			tmp_inName.InOut_label=1;
			getFileName(tmp_inName);
			if(rc_label==1)
			{
				posfilename[2]='r';
				posfilename[3]='c';
			}

			FILE* fp_pos_out;
			fp_pos_out=fopen(posfilename,"rb");
			if(fp_pos_out==NULL)
			{
				cout << "error: open candidate file failed!" << endl;
			}
			p_candidate_current_out=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			x=fread(p_candidate_current_out,sizeof(uint64_t),R_len/64+((R_len%64)?1:0),fp_pos_out);
			if(x!=R_len/64+((R_len%64)?1:0))
			{
				cout << "error: read current candidate failed!" << endl;
			}
			fclose(fp_pos_out);
			free(posfilename);
		}

		uint64_t *p_candidate_current=NULL;
		if(candidate_label==1)
		{
			p_candidate_current=(uint64_t*)malloc(sizeof(uint64_t)*(R_len/64+1));
			uint32_t l_shift_size=kmer_len-1;
			uint32_t r_shift_size=64-l_shift_size;

			for(uint32_t i=0;i<R_len/64;i++)
			{
				p_candidate_current[i]=p_candidate_current_in[i]|((p_candidate_current_out[i]<<l_shift_size)|(p_candidate_current_out[i+1]>>r_shift_size));
			}
			p_candidate_current[R_len/64]=p_candidate_current_in[R_len/64]|(p_candidate_current_out[R_len/64]<<l_shift_size);
			free(p_candidate_current_in);
			free(p_candidate_current_out);
			p_candidate_current_in=NULL;
			p_candidate_current_out=NULL;
		}
		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		struct para_Find_Branch* para_unipath;
		para_unipath=(struct para_Find_Branch*)malloc(sizeof(struct para_Find_Branch)*thread_num);

		uint64_t unit_len=R_len/thread_num+1;
		uint64_t res=unit_len%64;
		if(res!=0)
		{
			unit_len=unit_len+(64-res);
		}

		for(uint32_t j=0;j<thread_num;j++)
		{
			para_unipath[j].kmer_len=kmer_len;
			para_unipath[j].para=para;
			para_unipath[j].pos_label=pos_label;
			para_unipath[j].p_candidate_current_in=p_candidate_current_in;
			para_unipath[j].p_candidate_current_out=p_candidate_current_out;
			para_unipath[j].p_candidate_current=p_candidate_current;
			para_unipath[j].p_root=p_root;
			para_unipath[j].p_tmp=p_tmp;
			para_unipath[j].template_64=template_64;
			para_unipath[j].start=j*unit_len;
			para_unipath[j].p_unipath=p_root_unipath_arr;
			para_unipath[j].end=(j+1)*unit_len-1;
			para_unipath[j].thread_id=j;
			para_unipath[j].ref_i=ref_i;
			para_unipath[j].ref_num=ref_num;
			para_unipath[j].candidate_label=candidate_label;
			para_unipath[j].rc_label=rc_label;


			if(j!=thread_num-1)
			{
				uint32_t tmp_i;
				tmp_i=0;
				while(1)
				{
					cal_hash_value_directly_256bit(p_tmp+para_unipath[j].end+tmp_i,seq_k_current,para);
					if(search256BitHashTable(p_root,seq_k_current,para)!=0)
					{
						para_unipath[j].end=para_unipath[j].end+tmp_i;
						break;
					}
					tmp_i++;
				}
			}
			if(para_unipath[j].end>R_len-kmer_len)
			{
				para_unipath[j].end=R_len-kmer_len;
			}

			if(pthread_create(t+j, NULL, Find_Branch, (void*)(para_unipath+j))!=0)
			{
				cout << "error!" << endl;
			}
		}
		for(uint32_t j=0;j<thread_num;j++)
		{
			pthread_join(t[j], NULL);
		}
		free(t);
		free(para_unipath);
		if(pos_label==2||pos_label==3)
		{
			free(p_candidate_current_in);
			free(p_candidate_current_out);
		}
	}
	if(thread_num==1)
	{
		tmp.p_root_unipath[0]=p_root_unipath;
	}
	else
	{
		for(uint32_t i=0;i<thread_num;i++)
		{
			tmp.p_root_unipath[i]=p_root_unipath_arr[i];
		}
	}

}
void FindUniPath_Novel(struct para_findUnipath tmps)
{
	char* seq=tmps.seq;
	uint32_t kmer_len=tmps.original_kmer_length;
	uint32_t thread_num=tmps.thread_num;
	uint32_t pos_label=tmps.pos_label;
	uint32_t rc_label=tmps.rc_label;
	uint32_t pl_label=tmps.pl_label;
	uint32_t candidate_label=0;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	if(para.remainer1to64==0)
	{
		para.codefor1=0;
	}
	else
	{
		para.codefor1=0;
		for(uint32_t i=0;i<para.remainer1to64-1;i++)
		{
			para.codefor1=para.codefor1|1;
			para.codefor1=para.codefor1<<1;
		}
		para.codefor1=para.codefor1|1;
	}

	struct bit256KmerPara para_edge;
	para_edge.kmer1Len=kmer_len*2+2;
	para_edge.remainer1to64=para_edge.kmer1Len%64;
	para_edge.kmer64Len=para_edge.kmer1Len/64+(para_edge.remainer1to64?1:0);
	if(para_edge.remainer1to64==0)
	{
		para_edge.codefor1=0;
	}
	else
	{
		para_edge.codefor1=0;
		for(uint32_t i=0;i<para_edge.remainer1to64-1;i++)
		{
			para_edge.codefor1=para_edge.codefor1|1;
			para_edge.codefor1=para_edge.codefor1<<1;
		}
		para_edge.codefor1=para_edge.codefor1|1;
	}

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	char* branchfilename;
	branchfilename=(char*)malloc(sizeof(char)*20);
	struct para_getN tmp;
	tmp.kmerlen=kmer_len;
	tmp.FileName=branchfilename;
	tmp.InOut_label=0;
	tmp.isposition=0;
	tmp.isMiddle=1;
	getFileName(tmp);

	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*20);
	struct para_getN tmp_para;
	tmp_para.kmerlen=kmer_len;
	tmp_para.FileName=inputFileName;
	tmp_para.InOut_label=0;
	tmp_para.isposition=0;
	tmp_para.isMiddle=1;
	getFileName(tmp_para);

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
	fseek(int_input,0,0);
	free(inputFileName);
	cout << total_hash_number << endl;

	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*para.kmer64Len);
	uint64_t x;
	uint32_t task_size_loop;
	task_size_loop=total_hash_number;
	x=fread(h1,sizeof(uint64_t),task_size_loop*para.kmer64Len,int_input);
	if(x!=task_size_loop*para.kmer64Len)
	{
		cout << "error: read branching kmer wrong!" << endl;
	}

	struct bit256Hash* p_root=NULL;
	p_root=initial256BitHashTable();
	if(p_root==NULL)
	{
		cerr<< "error: initial BitHashTable failed!" << endl;
	}
	struct timeval tvs_hs,tve_hs;
	gettimeofday(&tvs_hs,NULL);
	if(thread_num==1)
	{
		for(uint64_t j=0;j<task_size_loop;j++)
		{
			insert256BitHashTable(p_root,h1+j*para.kmer64Len,j+1,para);
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
			p[i].Tasksize_current=task_size_loop;
			p[i].para1=para;
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
//	struct timeval tvs_hs_c,tve_hs_c;
//	gettimeofday(&tvs_hs_c,NULL);
//	for(uint64_t j=0;j<task_size_loop;j++)
//	{
//		uint64_t tmp;
//		tmp=search256BitHashTable(p_root,h1+j*para.kmer64Len,para);
//		if(tmp!=j+1)
//		{
//			cout << j << endl;
//			cout << "error!" << endl;
//			char a;
//			cin >> a;
//		}
//	}
//	fclose(int_input);
//	gettimeofday(&tve_hs_c,NULL);
//	double span_hs_c = tve_hs_c.tv_sec-tvs_hs_c.tv_sec + (tve_hs_c.tv_usec-tvs_hs_c.tv_usec)/1000000.0;
//	cout << "check hash time is: "<<span_hs_c<<endl;
//	cerr << "check hash kmer over!" << endl;
	free(h1);

	char * p_tmp;
	char* dataset=seq;
	struct bit256Hash** p_root_unipath;
	p_root_unipath=(struct bit256Hash**)malloc(sizeof(struct bit256Hash*)*thread_num);

	uint64_t R_len;
	cout << "start finding unipath ..." << endl;
	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	uint64_t template_64[64];
	template_64[0]=1;
	for(uint32_t i=0;i<63;i++)
	{
		template_64[i+1]=template_64[i]<<1;
	}

	inputFileName=(char*)malloc(sizeof(char)*20);
	char init[8]="unipath";
	char swap[4];
	strcpy(inputFileName,init);
	if(kmer_len<100)
	{
		inputFileName[7]='0';
		sprintf(swap,"%d",kmer_len);
		strcpy(inputFileName+8,swap);
	}
	else
	{
		sprintf(swap,"%d",kmer_len);
		strcpy(inputFileName+7,swap);
	}
	inputFileName[10]='\0';
	FILE* outUnipath;
	outUnipath=fopen(inputFileName,"wb+");

	FILE* outUnipath_c;
	outUnipath_c=fopen("unipath_counter","wb+");
	free(inputFileName);

	for(uint64_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		cout << ref_i+1 << ":" << p_ref_path.NumberOfPathes << endl;
		ReadSeq(&p_tmp,&R_len,p_ref_path.pRefFilePath[ref_i]);

		struct para_CalBranchVec para_branch_vec;
		para_branch_vec.R_len=R_len;
		para_branch_vec.kmer_len=kmer_len;
		para_branch_vec.p_tmp=p_tmp;
		para_branch_vec.pos_label=pos_label;
		para_branch_vec.rc_label=0;
		para_branch_vec.ref_i=ref_i;
		para_branch_vec.thread_num=thread_num;
		para_branch_vec.p_root=p_root;
		para_branch_vec.candidate_label=candidate_label;
		para_branch_vec.p_root_unipath=p_root_unipath;
		para_branch_vec.ref_num=p_ref_path.NumberOfPathes;

		CalBranchVec(para_branch_vec);

		for(uint32_t j=0;j<thread_num;j++)
		{
			for(uint32_t i=0;i<HashSize;i++)
			{
				if(p_root_unipath[j]->len[i]!=0)
				{
					fwrite(p_root_unipath[j]->p_pos[i],sizeof(uint64_t),p_root_unipath[j]->len[i],outUnipath);
				}
			}
		}
		for(uint32_t j=0;j<thread_num;j++)
		{
			freeHash256BitTable(p_root_unipath[j]);
		}

		if(rc_label==1)
		{
			cout << ref_i+1 << "RC:" << p_ref_path.NumberOfPathes << endl;
			rc(&p_tmp,R_len);

			struct para_CalBranchVec para_branch_vec;
			para_branch_vec.R_len=R_len;
			para_branch_vec.kmer_len=kmer_len;
			para_branch_vec.p_tmp=p_tmp;
			para_branch_vec.pos_label=pos_label;
			para_branch_vec.rc_label=1;
			para_branch_vec.ref_i=ref_i;
			para_branch_vec.thread_num=thread_num;
			para_branch_vec.p_root=p_root;
			para_branch_vec.candidate_label=candidate_label;
			para_branch_vec.p_root_unipath=p_root_unipath;
			para_branch_vec.ref_num=p_ref_path.NumberOfPathes;

			CalBranchVec(para_branch_vec);

			for(uint32_t j=0;j<thread_num;j++)
			{
				for(uint32_t i=0;i<HashSize;i++)
				{
					if(p_root_unipath[j]->len[i]!=0)
					{
						fwrite(p_root_unipath[j]->p_pos[i],sizeof(uint64_t),p_root_unipath[j]->len[i],outUnipath);
					}
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				freeHash256BitTable(p_root_unipath[j]);
			}
		}
		free(p_tmp);
	}
	cout << "unipath and position list are ready!" << endl;

	fclose(outUnipath);
	fclose(outUnipath_c);
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the finding time is: "<<span<< endl;

}
