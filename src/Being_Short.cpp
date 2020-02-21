#include "Being_Short.h"

void *judge_thread_being(void *arg)
{
	struct reslt_para_being* tmp=(struct reslt_para_being*)arg;

	uint64_t start=tmp->start;
	uint64_t task_len=tmp->task_len;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	uint8_t *p_hash_table=*(tmp->p_hash_table);

	uint32_t add_size;
	add_size=pow(2,20);

	uint64_t * reslt_in;
	uint64_t * reslt_out;
	reslt_in=(uint64_t *)malloc(sizeof(uint64_t)*add_size);
	reslt_out=(uint64_t *)malloc(sizeof(uint64_t)*add_size);

	uint8_t * reslt_in_ad;
	uint8_t * reslt_out_ad;

	if(tmp->ad_label==1)//20200207
    {
        reslt_in_ad=(uint8_t *)malloc(sizeof(uint8_t)*add_size);
        reslt_out_ad=(uint8_t *)malloc(sizeof(uint8_t)*add_size);
    }

	uint64_t reslt_size_in=0;
	uint64_t reslt_size_out=0;
	uint64_t reslt_Msize_in=add_size;
	uint64_t reslt_Msize_out=add_size;

	uint64_t start_thread;
	uint64_t end_thread;

	start_thread=start+thread_id*(task_len/thread_num);
	end_thread=start+(thread_id+1)*(task_len/thread_num)-1;
	if(thread_id==thread_num-1)
	{
		end_thread=start+task_len-1;
	}

	for(uint64_t i=start_thread;i<=end_thread;i++)
	{
		uint32_t tempkind;
		tempkind=Judge(p_hash_table[i]);
		if(tempkind==1||tempkind==3){
			if(reslt_size_in<reslt_Msize_in)
			{
				reslt_in[reslt_size_in]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad[reslt_size_in]=(p_hash_table[i]>>4)&15;
                }
				reslt_size_in++;
			}
			else
			{
				reslt_in=(uint64_t *)realloc(reslt_in,sizeof(uint64_t)*(reslt_Msize_in+add_size));

				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad=(uint8_t *)realloc(reslt_in_ad,sizeof(uint8_t)*(reslt_Msize_in+add_size));
                }
				reslt_Msize_in=reslt_Msize_in+add_size;

				reslt_in[reslt_size_in]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad[reslt_size_in]=(p_hash_table[i]>>4)&15;
                }
				reslt_size_in++;
			}
		}
		if(tempkind==2||tempkind==3)
		{
			if(reslt_size_out<reslt_Msize_out)
			{
				reslt_out[reslt_size_out]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad[reslt_size_out]=p_hash_table[i]&15;
                }
                reslt_size_out++;
			}
			else
			{
				reslt_out=(uint64_t *)realloc(reslt_out,sizeof(uint64_t)*(reslt_Msize_out+add_size));
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad=(uint8_t *)realloc(reslt_out_ad,sizeof(uint8_t)*(reslt_Msize_out+add_size));
                }
				reslt_Msize_out=reslt_Msize_out+add_size;

				reslt_out[reslt_size_out]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad[reslt_size_out]=p_hash_table[i]&15;
                }
				reslt_size_out++;
			}
		}
	}

	*(tmp->reslt_size_in)=reslt_size_in;
	*(tmp->reslt_size_out)=reslt_size_out;

	*(tmp->reslt_in)=reslt_in;
	*(tmp->reslt_out)=reslt_out;
	if(tmp->ad_label==1)
    {
        *(tmp->reslt_in_ad)=reslt_in_ad;
        *(tmp->reslt_out_ad)=reslt_out_ad;
    }
	return NULL;
}
void *judge_thread_being_all(void *arg)
{
	struct reslt_para_being* tmp=(struct reslt_para_being*)arg;

	uint64_t start=tmp->start;
	uint64_t task_len=tmp->task_len;
	uint32_t thread_id=tmp->thread_id;
	uint32_t thread_num=tmp->thread_num;
	uint8_t *p_hash_table=*(tmp->p_hash_table);

	uint32_t add_size;
	add_size=pow(2,20);

	uint64_t * reslt_in;
	uint64_t * reslt_out;
	uint64_t * reslt_uni;
	reslt_in=(uint64_t *)malloc(sizeof(uint64_t)*add_size);
	reslt_out=(uint64_t *)malloc(sizeof(uint64_t)*add_size);
	reslt_uni=(uint64_t *)malloc(sizeof(uint64_t)*add_size);

	uint8_t * reslt_in_ad;
	uint8_t * reslt_out_ad;
	uint8_t * reslt_uni_ad;

	if(tmp->ad_label==1)//20200207
    {
        reslt_in_ad=(uint8_t *)malloc(sizeof(uint8_t)*add_size);
        reslt_out_ad=(uint8_t *)malloc(sizeof(uint8_t)*add_size);
        reslt_uni_ad=(uint8_t *)malloc(sizeof(uint8_t)*add_size);
    }

	uint64_t reslt_size_in=0;
	uint64_t reslt_size_out=0;
	uint64_t reslt_size_uni=0;
	uint64_t reslt_Msize_in=add_size;
	uint64_t reslt_Msize_out=add_size;
	uint64_t reslt_Msize_uni=add_size;

	uint64_t start_thread;
	uint64_t end_thread;

	start_thread=start+thread_id*(task_len/thread_num);
	end_thread=start+(thread_id+1)*(task_len/thread_num)-1;
	if(thread_id==thread_num-1)
	{
		end_thread=start+task_len-1;
	}

	for(uint64_t i=start_thread;i<=end_thread;i++)
	{
		uint32_t tempkind;
		tempkind=Judge_novel(p_hash_table[i]);
		if(tempkind==0)
		{
			if(reslt_size_uni<reslt_Msize_uni)
			{
				reslt_uni[reslt_size_uni]=i;
				if(tmp->ad_label==1)//20200207
				{
					reslt_uni_ad[reslt_size_uni]=p_hash_table[i];
				}
				reslt_size_uni++;
			}
			else
			{
				reslt_uni=(uint64_t *)realloc(reslt_uni,sizeof(uint64_t)*(reslt_Msize_uni+add_size));

				if(tmp->ad_label==1)//20200207
				{
					reslt_uni_ad=(uint8_t *)realloc(reslt_uni_ad,sizeof(uint8_t)*(reslt_Msize_uni+add_size));
				}
				reslt_Msize_uni=reslt_Msize_uni+add_size;

				reslt_uni[reslt_size_uni]=i;
				if(tmp->ad_label==1)//20200207
				{
					reslt_uni_ad[reslt_size_uni]=p_hash_table[i];
				}
				reslt_size_uni++;
			}

		}
		if(tempkind==1||tempkind==3){
			if(reslt_size_in<reslt_Msize_in)
			{
				reslt_in[reslt_size_in]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad[reslt_size_in]=(p_hash_table[i]>>4)&15;
                }
				reslt_size_in++;
			}
			else
			{
				reslt_in=(uint64_t *)realloc(reslt_in,sizeof(uint64_t)*(reslt_Msize_in+add_size));

				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad=(uint8_t *)realloc(reslt_in_ad,sizeof(uint8_t)*(reslt_Msize_in+add_size));
                }
				reslt_Msize_in=reslt_Msize_in+add_size;

				reslt_in[reslt_size_in]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_in_ad[reslt_size_in]=(p_hash_table[i]>>4)&15;
                }
				reslt_size_in++;
			}
		}
		if(tempkind==2||tempkind==3)
		{
			if(reslt_size_out<reslt_Msize_out)
			{
				reslt_out[reslt_size_out]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad[reslt_size_out]=p_hash_table[i]&15;
                }
                reslt_size_out++;
			}
			else
			{
				reslt_out=(uint64_t *)realloc(reslt_out,sizeof(uint64_t)*(reslt_Msize_out+add_size));
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad=(uint8_t *)realloc(reslt_out_ad,sizeof(uint8_t)*(reslt_Msize_out+add_size));
                }
				reslt_Msize_out=reslt_Msize_out+add_size;

				reslt_out[reslt_size_out]=i;
				if(tmp->ad_label==1)//20200207
                {
                    reslt_out_ad[reslt_size_out]=p_hash_table[i]&15;
                }
				reslt_size_out++;
			}
		}
	}

	*(tmp->reslt_size_in)=reslt_size_in;
	*(tmp->reslt_size_out)=reslt_size_out;
	*(tmp->reslt_size_uni)=reslt_size_uni;

	*(tmp->reslt_in)=reslt_in;
	*(tmp->reslt_out)=reslt_out;
	*(tmp->reslt_uni)=reslt_uni;
	if(tmp->ad_label==1)
    {
        *(tmp->reslt_in_ad)=reslt_in_ad;
        *(tmp->reslt_out_ad)=reslt_out_ad;
        *(tmp->reslt_uni_ad)=reslt_uni_ad;
    }
	return NULL;
}
void *Parallel_cal_hashValue(void *arg)
{
	struct Parallel_para_being * p=(struct Parallel_para_being *)arg;
	char* p_seq=p->p_seq;
	uint64_t* p_hashValueTable=p->p_hashValueTable;
	uint32_t kmer_len=p->kmer_len;
	uint32_t thread_id=p->thread_id;
	uint32_t thread_num=p->thread_num;
	uint32_t start=p->start;
	uint32_t len=p->len;

	uint64_t start_thread;
	uint64_t end_thread;
	start_thread=start+(len/thread_num)*(thread_id-1);
	if(thread_id==thread_num)
	{
		end_thread=start+len-1;
	}
	else
	{
		end_thread=start+(len/thread_num)*(thread_id)-1;
	}
	uint64_t current;
	current=cal_hash_value_directly(p_seq+start_thread,kmer_len);
	p_hashValueTable[start_thread-start]=current;
	for(uint64_t i=start_thread+1;i<=end_thread;i++)
	{
		current=cal_hash_value_indirectly(p_seq+i,current,kmer_len);
		p_hashValueTable[i-start]=current;
	}
	return NULL;
}
void *Parallel_assign_ad(void *arg)
{
	struct Parallel_para_being * p=(struct Parallel_para_being *)arg;
	char* p_seq=p->p_seq;
	uint8_t* p_hashtable=p->p_hashtable;
	uint64_t seq_length=p->seq_length;
	uint64_t* p_hashValueTable=p->p_hashValueTable;
	uint32_t kmer_len=p->kmer_len;
	uint64_t thread_id=p->thread_id;
	uint64_t start=p->start;
	uint32_t len=p->len;

	if(thread_id<5)
	{
		char c;
		if(thread_id==1)
		{
			c='A';
		}
		if(thread_id==2)
		{
			c='C';
		}
		if(thread_id==3)
		{
			c='G';
		}
		if(thread_id==4)
		{
			c='T';
		}

		if(start==0)
		{
			if(p_seq[start]==c)
			{
				Transition('0',p_seq[start+kmer_len],(p_hashtable+p_hashValueTable[0]));
			}
			if(thread_id==1&&p_seq[start]!='A'&&p_seq[start]!='C'&&p_seq[start]!='G'&&p_seq[start]!='T')
			{
				Transition('0',p_seq[start+kmer_len],(p_hashtable+p_hashValueTable[0]));
			}
			for(uint64_t i=1;i<len;i++)
			{
				if(p_seq[start+i]==c)
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
				if(thread_id==1&&p_seq[start+i]!='A'&&p_seq[start+i]!='C'&&p_seq[start+i]!='G'&&p_seq[start+i]!='T')
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
			}
		}
		else if(start+len-1==seq_length-kmer_len)
		{
			for(uint64_t i=0;i<len-1;i++)
			{
				if(p_seq[start+i]==c)
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
				if(thread_id==1&&p_seq[start+i]!='A'&&p_seq[start+i]!='C'&&p_seq[start+i]!='G'&&p_seq[start+i]!='T')
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
			}
			if(p_seq[start+len-1]==c)
			{
				Transition(p_seq[start+len-1-1],'0',(p_hashtable+p_hashValueTable[len-1]));
			}
			if(thread_id==1&&p_seq[start+len-1]!='A'&&p_seq[start+len-1]!='C'&&p_seq[start+len-1]!='G'&&p_seq[start+len-1]!='T')
			{
				Transition(p_seq[start+len-1-1],'0',(p_hashtable+p_hashValueTable[len-1]));
			}
		}
		else
		{
			for(uint64_t i=0;i<len;i++)
			{
				if(p_seq[start+i]==c)
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
				if(thread_id==1&&p_seq[start+i]!='A'&&p_seq[start+i]!='C'&&p_seq[start+i]!='G'&&p_seq[start+i]!='T')
				{
					Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
				}
			}
		}
	}
	return NULL;
}
void *Parallel_assign_ad_mode(void *arg)
{
	struct Parallel_para_being * p=(struct Parallel_para_being *)arg;
	char* p_seq=p->p_seq;
	uint8_t* p_hashtable=p->p_hashtable;
	uint64_t seq_length=p->seq_length;
	uint64_t* p_hashValueTable=p->p_hashValueTable;
	uint32_t kmer_len=p->kmer_len;
	uint64_t thread_id=p->thread_id;
	uint64_t thread_num=p->thread_num;
	uint64_t start=p->start;
	uint32_t len=p->len;

	if(start==0)
	{
		if(p_hashValueTable[0]%thread_num==thread_id-1)
		{
			Transition('0',p_seq[start+kmer_len],(p_hashtable+p_hashValueTable[0]));
		}
		for(uint64_t i=1;i<len;i++)
		{
			if(p_hashValueTable[i]%thread_num==thread_id-1)
			{
				Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
			}
		}
	}
	else if(start+len-1==seq_length-kmer_len)
	{
		for(uint64_t i=0;i<len-1;i++)
		{
			if(p_hashValueTable[i]%thread_num==thread_id-1)
			{
				Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
			}
		}
		if(p_hashValueTable[len-1]%thread_num==thread_id-1)
		{
			Transition(p_seq[start+len-1-1],'0',(p_hashtable+p_hashValueTable[len-1]));
		}
	}
	else
	{
		for(uint64_t i=0;i<len;i++)
		{
			if(p_hashValueTable[i]%thread_num==thread_id-1)
			{
				Transition(p_seq[start+i-1],p_seq[start+i+kmer_len],(p_hashtable+p_hashValueTable[i]));
			}
		}
	}

	return NULL;
}
void Find_Short_Branched_Kmer(struct para_being tmp)
{
	char* seq;
	char* dataset=tmp.seq;
	uint64_t seq_length;
//	uint32_t refNumber=tmp.seq_length;
	uint64_t kmer_len_being=tmp.kmer_len_being;
	char *outputFileName=tmp.outputFileName;
	uint32_t thread_num=tmp.thread_num;
	uint32_t Task_Size=tmp.Task_Size;

	Task_Size=pow(2,Task_Size);
	uint64_t kmer_length;
	kmer_length=kmer_len_being;
	cout << "the k-mer length is: " << kmer_length << endl;

	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);

	uint64_t hash_table_length=pow(4,kmer_length);
	uint8_t *p_hash_table=(uint8_t *)malloc(hash_table_length*sizeof(uint8_t));
	memset(p_hash_table,0,hash_table_length);

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);


	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		cout << p_ref_path.pRefFilePath[ref_i]<< endl;
		if(thread_num>1)
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			uint64_t *p_hashvalueTable;
			uint64_t p_hashvalueTable_length=Task_Size;
			p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);

			uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
			uint64_t start_loop=0;
			uint64_t end_loop=0;
			struct Parallel_para_being *para_each;
			para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
			for(uint32_t i=0;i<loop_num;i++)
			{
				cout << i << ":" << loop_num << endl;
				if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
				{
					end_loop=start_loop+p_hashvalueTable_length-1;
				}
				else
				{
					end_loop=seq_length-kmer_length;
				}
				//1.1)
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].kmer_len=kmer_length;
					para_each[i].p_seq=seq;
					para_each[i].seq_length=seq_length;
					para_each[i].thread_id=i+1;
					para_each[i].thread_num=thread_num;
					para_each[i].start=start_loop;
					para_each[i].len=end_loop-start_loop+1;
					para_each[i].p_hashValueTable=p_hashvalueTable;
					if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].p_hashtable=p_hash_table;
					if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				start_loop=end_loop+1;
			}
			free(para_each);
			free(t);
			free(p_hashvalueTable);
		}
		else
		{
			uint64_t current;
			current=cal_hash_value_directly(seq,kmer_length);
			Transition('0',seq[kmer_length],(p_hash_table+current));

			for(uint64_t i=1;i<seq_length-kmer_length;i++)
			{
				current=cal_hash_value_indirectly(seq+i,current,kmer_length);
				Transition(seq[i-1],seq[i+kmer_length],(p_hash_table+current));
			}
			current=cal_hash_value_indirectly(seq+seq_length-kmer_length,current,kmer_length);
			Transition(seq[seq_length-kmer_length-1],'0',(p_hash_table+current));
		}
		free(seq);
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"HashTable construct time is: "<<span << endl;

	struct timeval tvs_check,tve_check;
	gettimeofday(&tvs_check,NULL);


	uint64_t sum_being=0,sum_in=0,sum_out=0;
	FILE *fp_in,*fp_out;
	char* inc;
	inc=(char*)malloc(sizeof(char)*10);
	char* outc;
	outc=(char*)malloc(sizeof(char)*10);
	char in_add[3]="in";
	char out_add[4]="out";
	strcpy(inc,outputFileName);
	strcpy(inc+6,in_add);
	inc[8]='\0';
	strcpy(outc,outputFileName);
	strcpy(outc+6,out_add);
	inc[9]='\0';
	fp_in=fopen(inc,"w+");
	fp_out=fopen(outc,"w+");
	cout <<"start checking hash table for branched k-mer."<< endl;

	if(thread_num>1)
	{
		uint64_t **reslt_in;
		uint64_t **reslt_out;
		reslt_in=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_out=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

		uint64_t * reslt_size_in;
		reslt_size_in=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);
		uint64_t * reslt_size_out;
		reslt_size_out=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

		struct reslt_para_being* tmp;
		tmp=(struct reslt_para_being*)malloc(sizeof(struct reslt_para_being)*thread_num);

		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		uint64_t unit_size=Task_Size;
		uint32_t loop_num=hash_table_length/unit_size+hash_table_length%unit_size;
		for(uint32_t i=0;i<loop_num;i++)
		{
			uint64_t start_loop=unit_size*i;
			uint64_t end_loop;
			if(i==loop_num-1)
			{
				end_loop=hash_table_length-1;
			}
			else
			{
				end_loop=unit_size*(i+1)-1;
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=start_loop;
				tmp[j].task_len=end_loop-start_loop+1;
				tmp[j].thread_id=j;
				tmp[j].thread_num=thread_num;
				tmp[j].reslt_in=&(reslt_in[j]);
				tmp[j].reslt_out=&(reslt_out[j]);
				tmp[j].reslt_size_in=&(reslt_size_in[j]);
				tmp[j].reslt_size_out=&(reslt_size_out[j]);
				tmp[j].p_hash_table=&p_hash_table;

				if(pthread_create(t+j, NULL, judge_thread_being, (void*)(tmp+j))!=0)
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
				for(uint64_t k=0;k<reslt_size_in[j];k++)
				{
					fprintf(fp_in,"%ld\n", reslt_in[j][k]);
					sum_in++;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				for(uint64_t k=0;k<reslt_size_out[j];k++)
				{
					fprintf(fp_out,"%ld\n",reslt_out[j][k]);
					sum_out++;
				}
			}
		}

		free(t);
		free(reslt_size_in);
		free(reslt_size_out);
		for(uint32_t i=0;i<thread_num;i++)
		{
			free(reslt_in[i]);
			free(reslt_out[i]);
		}
		free(reslt_in);
		free(reslt_out);
		free(tmp);
	}
	else
	{
		for(uint64_t i=0;i<hash_table_length;i++)
		{
			uint32_t tempkind;
			tempkind=Judge(p_hash_table[i]);
			if(tempkind==1||tempkind==3){
				fprintf(fp_in,"%ld\n",i);
				sum_in++;
			}
			if(tempkind==2||tempkind==3)
			{
				fprintf(fp_out,"%ld\n",i);
				sum_out++;
			}
			if(p_hash_table[i]!=0){
				sum_being++;
			}
		}
	}
	cout <<"the # of being k-mer is: ."<<sum_being << endl;
	cout <<"the # of in-branched k-mer is: "<<sum_in << endl;
	cout <<"the # of out-branched k-mer is: "<<sum_out << endl;
//	remove(inc);
//	remove(outc);
	free(inc);
	free(outc);
	fclose(fp_in);
	fclose(fp_out);

	free(p_hash_table);
	gettimeofday(&tve_check,NULL);
	double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
	cout <<"check time is: "<<span_check<< endl;
	gettimeofday(&tve_total,NULL);
	double span_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total time is: "<<span_total <<endl;
}
void Find_Short_Branched_Kmer_bit(struct para_being tmp)
{
	uint8_t* seq[4];
	uint64_t seq_length;
	char* dataset=tmp.seq;
//	uint32_t refNumber=tmp.seq_length;
	uint64_t kmer_len_being=tmp.kmer_len_being;
	char *outputFileName=tmp.outputFileName;
	uint32_t thread_num=tmp.thread_num;
	uint32_t Task_Size=tmp.Task_Size;

	struct bitKmerPara para;
	para.kmer1Len=kmer_len_being*2;
	cout << para.kmer1Len<<endl;
	cout << para.kmer1Len/8<<endl;

	para.remainer1to8=para.kmer1Len%8;
	cout << para.remainer1to8<<endl;
	para.kmer8Len=(para.kmer1Len/8)+(para.remainer1to8?1:0);
	cout << para.kmer8Len<<endl;
	para.remainer1to64=para.kmer1Len%64;
	cout << para.remainer1to64<<endl;
	para.kmer64Len=(para.kmer1Len/64)+(para.remainer1to64?1:0);
	cout << para.kmer64Len<<endl;
	para.remainer8to64=para.kmer8Len%8;
	cout << para.remainer8to64<<endl;

	Task_Size=pow(2,Task_Size);
	uint64_t kmer_length;
	kmer_length=kmer_len_being;
	cout << "the k-mer length is: " << kmer_length << endl;

	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);

	uint64_t hash_table_length=pow(4,kmer_length);
	uint8_t *p_hash_table=(uint8_t *)malloc(hash_table_length*sizeof(uint8_t));
	memset(p_hash_table,0,hash_table_length);

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);

	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq_bit(seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		cout << p_ref_path.pRefFilePath[ref_i]<< endl;
//		if(thread_num>1)
//		{
//			pthread_t* t;
//			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
//			uint64_t *p_hashvalueTable;
//			uint64_t p_hashvalueTable_length=Task_Size;
//			p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);
//
//			uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
//			uint64_t start_loop=0;
//			uint64_t end_loop=0;
//			struct Parallel_para_being *para_each;
//			para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
//			for(uint32_t i=0;i<loop_num;i++)
//			{
//				cout << i << ":" << loop_num << endl;
//				if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
//				{
//					end_loop=start_loop+p_hashvalueTable_length-1;
//				}
//				else
//				{
//					end_loop=seq_length-kmer_length;
//				}
//				//1.1)
//				for(uint32_t i=0;i<thread_num;i++)
//				{
//					para_each[i].kmer_len=kmer_length;
//					para_each[i].p_seq=seq;
//					para_each[i].seq_length=seq_length;
//					para_each[i].thread_id=i+1;
//					para_each[i].thread_num=thread_num;
//					para_each[i].start=start_loop;
//					para_each[i].len=end_loop-start_loop+1;
//					para_each[i].p_hashValueTable=p_hashvalueTable;
//					if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
//					{
//						cout << "error!" << endl;
//					}
//				}
//				for(uint32_t i=0;i<thread_num;i++)
//				{
//					pthread_join(t[i], NULL);
//				}
//				for(uint32_t i=0;i<thread_num;i++)
//				{
//					para_each[i].p_hashtable=p_hash_table;
//					if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
//					{
//						cout << "error!" << endl;
//					}
//				}
//				for(uint32_t i=0;i<thread_num;i++)
//				{
//					pthread_join(t[i], NULL);
//				}
//				start_loop=end_loop+1;
//			}
//			free(para_each);
//			free(t);
//			free(p_hashvalueTable);
//		}
//		else
//		{
			uint8_t current[32];
			uint64_t current_hv[4];
			Fetch_kmer_on_bit_sequence(seq,0,para,current);
			calMulDimHashValue(current,current_hv,para);
			uint64_t posr=0+kmer_length;
			uint64_t shiftr=posr%4;
			uint64_t idr=posr/4;
			uint8_t tmpr=seq[shiftr][idr];
			tmpr=tmpr>>6;
			uint64_t posl;
			uint64_t shiftl;
			uint64_t idl;
			uint8_t tmpl;

			Transition_bitchar(4,tmpr,p_hash_table+current_hv[0]);

			for(uint64_t i=1;i<seq_length-kmer_length;i++)
			{
				Fetch_kmer_on_bit_sequence(seq,i,para,current);
				calMulDimHashValue(current,current_hv,para);
				posr=i+kmer_length;
				shiftr=posr%4;
				idr=posr/4;
				tmpr=seq[shiftr][idr];
				tmpr=tmpr>>6;
				posl=i-1;
				shiftl=posl%4;
				idl=posl/4;
				tmpl=seq[shiftl][idl];
				tmpl=tmpl>>6;
				Transition_bitchar(tmpl,tmpr,p_hash_table+current_hv[0]);
			}
			Fetch_kmer_on_bit_sequence(seq,seq_length-kmer_length,para,current);
			calMulDimHashValue(current,current_hv,para);
			posl=seq_length-kmer_length-1;
			shiftl=posl%4;
			idl=posl/4;
			tmpl=seq[shiftl][idl];
			tmpl=tmpl>>6;
			Transition_bitchar(tmpl,4,p_hash_table+current_hv[0]);
//		}
		free(seq[0]);
		free(seq[1]);
		free(seq[2]);
		free(seq[3]);
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"HashTable construct time is: "<<span << endl;

	struct timeval tvs_check,tve_check;
	gettimeofday(&tvs_check,NULL);


	uint64_t sum_being=0,sum_in=0,sum_out=0;
	FILE *fp_in,*fp_out;
	char* inc;
	inc=(char*)malloc(sizeof(char)*10);
	char* outc;
	outc=(char*)malloc(sizeof(char)*10);
	char in_add[3]="in";
	char out_add[4]="out";
	strcpy(inc,outputFileName);
	strcpy(inc+6,in_add);
	inc[8]='\0';
	strcpy(outc,outputFileName);
	strcpy(outc+6,out_add);
	inc[9]='\0';
	fp_in=fopen(inc,"w+");
	fp_out=fopen(outc,"w+");
	cout <<"start checking hash table for branched k-mer."<< endl;

	if(thread_num>1)
	{
		uint64_t **reslt_in;
		uint64_t **reslt_out;
		reslt_in=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_out=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

		uint64_t * reslt_size_in;
		reslt_size_in=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);
		uint64_t * reslt_size_out;
		reslt_size_out=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

		struct reslt_para_being* tmp;
		tmp=(struct reslt_para_being*)malloc(sizeof(struct reslt_para_being)*thread_num);

		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		uint64_t unit_size=Task_Size;
		uint32_t loop_num=hash_table_length/unit_size+hash_table_length%unit_size;
		for(uint32_t i=0;i<loop_num;i++)
		{
			uint64_t start_loop=unit_size*i;
			uint64_t end_loop;
			if(i==loop_num-1)
			{
				end_loop=hash_table_length-1;
			}
			else
			{
				end_loop=unit_size*(i+1)-1;
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=start_loop;
				tmp[j].task_len=end_loop-start_loop+1;
				tmp[j].thread_id=j;
				tmp[j].thread_num=thread_num;
				tmp[j].reslt_in=&(reslt_in[j]);
				tmp[j].reslt_out=&(reslt_out[j]);
				tmp[j].reslt_size_in=&(reslt_size_in[j]);
				tmp[j].reslt_size_out=&(reslt_size_out[j]);
				tmp[j].p_hash_table=&p_hash_table;

				if(pthread_create(t+j, NULL, judge_thread_being, (void*)(tmp+j))!=0)
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
				for(uint64_t k=0;k<reslt_size_in[j];k++)
				{
					fprintf(fp_in,"%ld\n", reslt_in[j][k]);
					sum_in++;
				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				for(uint64_t k=0;k<reslt_size_out[j];k++)
				{
					fprintf(fp_out,"%ld\n",reslt_out[j][k]);
					sum_out++;
				}
			}
		}

		free(t);
		free(reslt_size_in);
		free(reslt_size_out);
		for(uint32_t i=0;i<thread_num;i++)
		{
			free(reslt_in[i]);
			free(reslt_out[i]);
		}
		free(reslt_in);
		free(reslt_out);
		free(tmp);
	}
	else
	{
		for(uint64_t i=0;i<hash_table_length;i++)
		{
			uint32_t tempkind;
			tempkind=Judge(p_hash_table[i]);
			if(tempkind==1||tempkind==3){
				fprintf(fp_in,"%ld\n",i);
				sum_in++;
			}
			if(tempkind==2||tempkind==3)
			{
				fprintf(fp_out,"%ld\n",i);
				sum_out++;
			}
			if(p_hash_table[i]!=0){
				sum_being++;
			}
		}
	}
	cout <<"the # of being k-mer is: ."<<sum_being << endl;
	cout <<"the # of in-branched k-mer is: "<<sum_in << endl;
	cout <<"the # of out-branched k-mer is: "<<sum_out << endl;
//	remove(inc);
//	remove(outc);
	free(inc);
	free(outc);
	fclose(fp_in);
	fclose(fp_out);

	free(p_hash_table);
	gettimeofday(&tve_check,NULL);
	double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
	cout <<"check time is: "<<span_check<< endl;
	gettimeofday(&tve_total,NULL);
	double span_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total time is: "<<span_total <<endl;
}
void Find_Short_Branched_Kmer_256bit_novel(struct para_being tmp)
{
	char* seq;
	char* dataset=tmp.seq;
	uint64_t seq_length;
	uint64_t kmer_len_being=tmp.kmer_len_being;
	char *outputFileName=tmp.outputFileName;
	uint32_t thread_num=tmp.thread_num;
	uint32_t Task_Size=tmp.Task_Size;
	uint32_t ad_label=tmp.ad_label;//20200207
	uint32_t rc_label=tmp.rc_label;//20200207

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len_being*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);

	Task_Size=pow(2,Task_Size);
	uint64_t kmer_length;
	kmer_length=kmer_len_being;
	cout << "the k-mer length is: " << kmer_length << endl;

	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);

	uint64_t hash_table_length=pow(4,kmer_length);
	uint8_t *p_hash_table=(uint8_t *)malloc(hash_table_length*sizeof(uint8_t));
	memset(p_hash_table,0,hash_table_length);

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);


	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		cout << p_ref_path.pRefFilePath[ref_i]<< endl;
		if(thread_num>1)
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			uint64_t *p_hashvalueTable;
			uint64_t p_hashvalueTable_length=Task_Size;
			p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);

			uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
			uint64_t start_loop=0;
			uint64_t end_loop=0;
			struct Parallel_para_being *para_each;
			para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
			for(uint32_t i=0;i<loop_num;i++)
			{
//				cout << i << ":" << loop_num << endl;
				if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
				{
					end_loop=start_loop+p_hashvalueTable_length-1;
				}
				else
				{
					end_loop=seq_length-kmer_length;
				}
				//1.1)
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].kmer_len=kmer_length;
					para_each[i].p_seq=seq;
					para_each[i].seq_length=seq_length;
					para_each[i].thread_id=i+1;
					para_each[i].thread_num=thread_num;
					para_each[i].start=start_loop;
					para_each[i].len=end_loop-start_loop+1;
					para_each[i].p_hashValueTable=p_hashvalueTable;
					if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].p_hashtable=p_hash_table;
					if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				start_loop=end_loop+1;
			}
			free(para_each);
			free(t);
			free(p_hashvalueTable);
		}
		else
		{
			uint64_t current[4];
			current[0]=cal_hash_value_directly(seq,kmer_length);
			Transition('0',seq[kmer_length],(p_hash_table+current[0]));

			for(uint64_t i=1;i<seq_length-kmer_length;i++)
			{
				current[0]=cal_hash_value_indirectly(seq+i,current[0],kmer_length);
				Transition(seq[i-1],seq[i+kmer_length],(p_hash_table+current[0]));
			}
			current[0]=cal_hash_value_indirectly(seq+seq_length-kmer_length,current[0],kmer_length);
			Transition(seq[seq_length-kmer_length-1],'0',(p_hash_table+current[0]));
		}
		if(rc_label==1)//20200207
        {
            rc(&seq,seq_length);
            cout << "rc:" << p_ref_path.pRefFilePath[ref_i]<< endl;
            if(thread_num>1)
            {
                pthread_t* t;
                t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                uint64_t *p_hashvalueTable;
                uint64_t p_hashvalueTable_length=Task_Size;
                p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);

                uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
                uint64_t start_loop=0;
                uint64_t end_loop=0;
                struct Parallel_para_being *para_each;
                para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
                for(uint32_t i=0;i<loop_num;i++)
                {
    //				cout << i << ":" << loop_num << endl;
                    if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
                    {
                        end_loop=start_loop+p_hashvalueTable_length-1;
                    }
                    else
                    {
                        end_loop=seq_length-kmer_length;
                    }
                    //1.1)
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_each[i].kmer_len=kmer_length;
                        para_each[i].p_seq=seq;
                        para_each[i].seq_length=seq_length;
                        para_each[i].thread_id=i+1;
                        para_each[i].thread_num=thread_num;
                        para_each[i].start=start_loop;
                        para_each[i].len=end_loop-start_loop+1;
                        para_each[i].p_hashValueTable=p_hashvalueTable;
                        if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
                        {
                            cout << "error!" << endl;
                        }
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        pthread_join(t[i], NULL);
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_each[i].p_hashtable=p_hash_table;
                        if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
                        {
                            cout << "error!" << endl;
                        }
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        pthread_join(t[i], NULL);
                    }
                    start_loop=end_loop+1;
                }
                free(para_each);
                free(t);
                free(p_hashvalueTable);
            }
            else
            {
                uint64_t current[4];
                current[0]=cal_hash_value_directly(seq,kmer_length);
                Transition('0',seq[kmer_length],(p_hash_table+current[0]));

                for(uint64_t i=1;i<seq_length-kmer_length;i++)
                {
                    current[0]=cal_hash_value_indirectly(seq+i,current[0],kmer_length);
                    Transition(seq[i-1],seq[i+kmer_length],(p_hash_table+current[0]));
                }
                current[0]=cal_hash_value_indirectly(seq+seq_length-kmer_length,current[0],kmer_length);
                Transition(seq[seq_length-kmer_length-1],'0',(p_hash_table+current[0]));
            }
        }
		free(seq);
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"HashTable construct time is: "<<span << endl;

	struct timeval tvs_check,tve_check;
	gettimeofday(&tvs_check,NULL);


	uint64_t sum_being=0,sum_in=0,sum_out=0,sum_uni=0;
	FILE *fp_in,*fp_out,*fp_uni;

	char* outputFileName_in;
	outputFileName_in=(char*)malloc(sizeof(char)*11);
	struct para_getN tmp_outname;
	tmp_outname.kmerlen=kmer_len_being;
	tmp_outname.FileName=outputFileName_in;
	tmp_outname.InOut_label=0;
	tmp_outname.isMiddle=0;
	tmp_outname.isposition=0;
	getFileName(tmp_outname);
	fp_in=fopen(outputFileName_in,"wb+");

	char* outputFileName_out;
	outputFileName_out=(char*)malloc(sizeof(char)*11);
	tmp_outname.kmerlen=kmer_len_being;
	tmp_outname.FileName=outputFileName_out;
	tmp_outname.InOut_label=1;
	tmp_outname.isMiddle=0;
	tmp_outname.isposition=0;
	getFileName(tmp_outname);
	fp_out=fopen(outputFileName_out,"wb+");

	outputFileName_out[7]='u';
	outputFileName_out[8]='n';
	outputFileName_out[9]='i';

	fp_uni=fopen(outputFileName_out,"wb+");

	FILE *fp_in_ad,*fp_out_ad,*fp_uni_ad;
	char* outputFileName_ad_in;
    struct para_getN tmp_ad_in;
    if(ad_label==1)
    {
        outputFileName_ad_in=(char*)malloc(sizeof(char)*13);
        tmp_ad_in.kmerlen=kmer_len_being;
        tmp_ad_in.FileName=outputFileName_ad_in;
        tmp_ad_in.InOut_label=0;
        tmp_ad_in.isMiddle=2;
        tmp_ad_in.isposition=0;
        getFileName(tmp_ad_in);
        fp_in_ad=fopen(outputFileName_ad_in,"ab+");
        free(outputFileName_ad_in);
    }

    char* outputFileName_ad_out;
    struct para_getN tmp_ad_out;
    if(ad_label==1)
    {
        outputFileName_ad_out=(char*)malloc(sizeof(char)*13);
        tmp_ad_out.kmerlen=kmer_len_being;
        tmp_ad_out.FileName=outputFileName_ad_out;
        tmp_ad_out.InOut_label=1;
        tmp_ad_out.isMiddle=2;
        tmp_ad_out.isposition=0;
        getFileName(tmp_ad_out);
        fp_out_ad=fopen(outputFileName_ad_out,"ab+");

    	outputFileName_ad_out[7]='u';
    	outputFileName_ad_out[8]='n';
    	outputFileName_ad_out[9]='i';

    	fp_uni_ad=fopen(outputFileName_ad_out,"wb+");
        free(outputFileName_ad_out);

    }


	cout <<"start checking hash table for branched k-mer."<< endl;

	if(thread_num>1)
	{
		uint64_t **reslt_in;
		uint64_t **reslt_out;
		uint64_t **reslt_uni;
		reslt_in=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_out=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_uni=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

		uint8_t **reslt_in_ad;
		uint8_t **reslt_out_ad;
		uint8_t **reslt_uni_ad;

		if(ad_label==1)
        {
            reslt_in_ad=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);
            reslt_out_ad=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);
            reslt_uni_ad=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);
        }

		uint64_t * reslt_size_in;
		reslt_size_in=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);
		uint64_t * reslt_size_out;
		reslt_size_out=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);
		uint64_t * reslt_size_uni;
		reslt_size_uni=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

		struct reslt_para_being* tmp;
		tmp=(struct reslt_para_being*)malloc(sizeof(struct reslt_para_being)*thread_num);

		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		uint64_t unit_size=Task_Size;
		uint32_t loop_num=hash_table_length/unit_size+hash_table_length%unit_size;
		for(uint32_t i=0;i<loop_num;i++)
		{
			uint64_t start_loop=unit_size*i;
			uint64_t end_loop;
			if(i==loop_num-1)
			{
				end_loop=hash_table_length-1;
			}
			else
			{
				end_loop=unit_size*(i+1)-1;
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=start_loop;
				tmp[j].task_len=end_loop-start_loop+1;
				tmp[j].thread_id=j;
				tmp[j].thread_num=thread_num;
				tmp[j].reslt_in=&(reslt_in[j]);
				tmp[j].reslt_out=&(reslt_out[j]);
				tmp[j].reslt_uni=&(reslt_uni[j]);
				if(ad_label==1)
                {
                    tmp[j].reslt_in_ad=&(reslt_in_ad[j]);
                    tmp[j].reslt_out_ad=&(reslt_out_ad[j]);
                    tmp[j].reslt_uni_ad=&(reslt_uni_ad[j]);
                }
				tmp[j].reslt_size_in=&(reslt_size_in[j]);
				tmp[j].reslt_size_out=&(reslt_size_out[j]);
				tmp[j].reslt_size_uni=&(reslt_size_uni[j]);
				tmp[j].p_hash_table=&p_hash_table;
				tmp[j].ad_label=ad_label;

				if(pthread_create(t+j, NULL, judge_thread_being_all, (void*)(tmp+j))!=0)
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
                fwrite(reslt_in[j],sizeof(uint64_t),reslt_size_in[j],fp_in);
                if(ad_label==1)
                {
                    fwrite(reslt_in_ad[j],sizeof(uint8_t),reslt_size_in[j],fp_in_ad);
                }
                sum_in+=reslt_size_in[j];
//				for(uint64_t k=0;k<reslt_size_in[j];k++)
//				{
//					fwrite(&(reslt_in[j][k]),sizeof(uint64_t),1,fp_in);
//					sum_in++;
//				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
			    fwrite(reslt_out[j],sizeof(uint64_t),reslt_size_out[j],fp_out);
                if(ad_label==1)
                {
                    fwrite(reslt_out_ad[j],sizeof(uint8_t),reslt_size_out[j],fp_out_ad);
                }
                sum_out+=reslt_size_out[j];
//				for(uint64_t k=0;k<reslt_size_out[j];k++)
//				{
//					fwrite(&(reslt_out[j][k]),sizeof(uint64_t),1,fp_out);
//					sum_out++;
//				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
                fwrite(reslt_uni[j],sizeof(uint64_t),reslt_size_uni[j],fp_uni);
                if(ad_label==1)
                {
                    fwrite(reslt_uni_ad[j],sizeof(uint8_t),reslt_size_uni[j],fp_uni_ad);
                }
                sum_uni+=reslt_size_uni[j];
//				for(uint64_t k=0;k<reslt_size_in[j];k++)
//				{
//					fwrite(&(reslt_in[j][k]),sizeof(uint64_t),1,fp_in);
//					sum_in++;
//				}
			}
		}

		free(t);
		free(reslt_size_in);
		free(reslt_size_out);
		free(reslt_size_uni);
		for(uint32_t i=0;i<thread_num;i++)
		{
			free(reslt_in[i]);
			free(reslt_out[i]);
			free(reslt_uni[i]);
		}
		free(reslt_in);
		free(reslt_out);
		free(reslt_uni);
		if(ad_label==1)
		{
			for(uint32_t i=0;i<thread_num;i++)
			{
				free(reslt_in_ad[i]);
				free(reslt_out_ad[i]);
				free(reslt_uni_ad[i]);
			}
			free(reslt_in_ad);
			free(reslt_out_ad);
			free(reslt_uni_ad);
		}
		free(tmp);
	}
	else
	{
		for(uint64_t i=0;i<hash_table_length;i++)
		{
			uint32_t tempkind;
			tempkind=Judge_novel(p_hash_table[i]);
			if(tempkind==0){
				fwrite(&i,sizeof(uint64_t),1,fp_uni);
				if(ad_label==1)
                {
                    uint8_t tmp;
                    tmp=p_hash_table[i];
                    fwrite(&tmp,sizeof(uint8_t),1,fp_uni_ad);
                }
				sum_uni++;
			}
			if(tempkind==1||tempkind==3){
				fwrite(&i,sizeof(uint64_t),1,fp_in);
				if(ad_label==1)
                {
                    uint8_t tmp;
                    tmp=(p_hash_table[i]>>4)&15;
                    fwrite(&tmp,sizeof(uint8_t),1,fp_in_ad);
                }
				sum_in++;
			}
			if(tempkind==2||tempkind==3)
			{
				fwrite(&i,sizeof(uint64_t),1,fp_out);
				if(ad_label==1)
                {
                    uint8_t tmp;
                    tmp=p_hash_table[i]&15;
                    fwrite(&tmp,sizeof(uint8_t),1,fp_out_ad);
                }
				sum_out++;
			}
//			if(p_hash_table[i]!=0){
//				sum_being++;
//			}
		}
	}
//	cout <<"the # of being k-mer is: "<<sum_being << endl;
	cout <<"the # of in-branched k-mer is: "<<sum_in << endl;
	cout <<"the # of out-branched k-mer is: "<<sum_out << endl;
	cout <<"the # of uni-branched k-mer is: "<<sum_out << endl;
	free(outputFileName_in);
	free(outputFileName_out);
	fclose(fp_in);
	fclose(fp_out);
	fclose(fp_uni);
	if(ad_label==1)
    {
        fclose(fp_in_ad);
        fclose(fp_out_ad);
        fclose(fp_uni_ad);
    }

	free(p_hash_table);
	gettimeofday(&tve_check,NULL);
	double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
	cout <<"check time is: "<<span_check<< endl;
	gettimeofday(&tve_total,NULL);
	double span_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total time is: "<<span_total <<endl;
}
void Find_Short_Branched_Kmer_256bit(struct para_being tmp)
{
	char* seq;
	char* dataset=tmp.seq;
	uint64_t seq_length;
	uint64_t kmer_len_being=tmp.kmer_len_being;
	char *outputFileName=tmp.outputFileName;
	uint32_t thread_num=tmp.thread_num;
	uint32_t Task_Size=tmp.Task_Size;
	uint32_t ad_label=tmp.ad_label;//20200207
	uint32_t rc_label=tmp.rc_label;//20200207

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len_being*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);

	Task_Size=pow(2,Task_Size);
	uint64_t kmer_length;
	kmer_length=kmer_len_being;
	cout << "the k-mer length is: " << kmer_length << endl;

	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);

	uint64_t hash_table_length=pow(4,kmer_length);
	uint8_t *p_hash_table=(uint8_t *)malloc(hash_table_length*sizeof(uint8_t));
	memset(p_hash_table,0,hash_table_length);

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	struct RefFilePath p_ref_path;
	getRefFilePathes(dataset, &p_ref_path);


	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		cout << p_ref_path.pRefFilePath[ref_i]<< endl;
		if(thread_num>1)
		{
			pthread_t* t;
			t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
			uint64_t *p_hashvalueTable;
			uint64_t p_hashvalueTable_length=Task_Size;
			p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);

			uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
			uint64_t start_loop=0;
			uint64_t end_loop=0;
			struct Parallel_para_being *para_each;
			para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
			for(uint32_t i=0;i<loop_num;i++)
			{
//				cout << i << ":" << loop_num << endl;
				if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
				{
					end_loop=start_loop+p_hashvalueTable_length-1;
				}
				else
				{
					end_loop=seq_length-kmer_length;
				}
				//1.1)
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].kmer_len=kmer_length;
					para_each[i].p_seq=seq;
					para_each[i].seq_length=seq_length;
					para_each[i].thread_id=i+1;
					para_each[i].thread_num=thread_num;
					para_each[i].start=start_loop;
					para_each[i].len=end_loop-start_loop+1;
					para_each[i].p_hashValueTable=p_hashvalueTable;
					if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					para_each[i].p_hashtable=p_hash_table;
					if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
					{
						cout << "error!" << endl;
					}
				}
				for(uint32_t i=0;i<thread_num;i++)
				{
					pthread_join(t[i], NULL);
				}
				start_loop=end_loop+1;
			}
			free(para_each);
			free(t);
			free(p_hashvalueTable);
		}
		else
		{
			uint64_t current[4];
			current[0]=cal_hash_value_directly(seq,kmer_length);
			Transition('0',seq[kmer_length],(p_hash_table+current[0]));

			for(uint64_t i=1;i<seq_length-kmer_length;i++)
			{
				current[0]=cal_hash_value_indirectly(seq+i,current[0],kmer_length);
				Transition(seq[i-1],seq[i+kmer_length],(p_hash_table+current[0]));
			}
			current[0]=cal_hash_value_indirectly(seq+seq_length-kmer_length,current[0],kmer_length);
			Transition(seq[seq_length-kmer_length-1],'0',(p_hash_table+current[0]));
		}
		if(rc_label==1)//20200207
        {
            rc(&seq,seq_length);
            cout << "rc:" << p_ref_path.pRefFilePath[ref_i]<< endl;
            if(thread_num>1)
            {
                pthread_t* t;
                t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
                uint64_t *p_hashvalueTable;
                uint64_t p_hashvalueTable_length=Task_Size;
                p_hashvalueTable=(uint64_t *)malloc(sizeof(uint64_t)*p_hashvalueTable_length);

                uint32_t loop_num=seq_length/p_hashvalueTable_length+1;
                uint64_t start_loop=0;
                uint64_t end_loop=0;
                struct Parallel_para_being *para_each;
                para_each=(struct Parallel_para_being *)malloc(sizeof(struct Parallel_para_being)*thread_num);
                for(uint32_t i=0;i<loop_num;i++)
                {
    //				cout << i << ":" << loop_num << endl;
                    if(start_loop+p_hashvalueTable_length<=seq_length-kmer_length)
                    {
                        end_loop=start_loop+p_hashvalueTable_length-1;
                    }
                    else
                    {
                        end_loop=seq_length-kmer_length;
                    }
                    //1.1)
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_each[i].kmer_len=kmer_length;
                        para_each[i].p_seq=seq;
                        para_each[i].seq_length=seq_length;
                        para_each[i].thread_id=i+1;
                        para_each[i].thread_num=thread_num;
                        para_each[i].start=start_loop;
                        para_each[i].len=end_loop-start_loop+1;
                        para_each[i].p_hashValueTable=p_hashvalueTable;
                        if(pthread_create(t+i, NULL, Parallel_cal_hashValue, (void*)(para_each+i))!=0)
                        {
                            cout << "error!" << endl;
                        }
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        pthread_join(t[i], NULL);
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        para_each[i].p_hashtable=p_hash_table;
                        if(pthread_create(t+i, NULL, Parallel_assign_ad_mode, (void*)(para_each+i))!=0)
                        {
                            cout << "error!" << endl;
                        }
                    }
                    for(uint32_t i=0;i<thread_num;i++)
                    {
                        pthread_join(t[i], NULL);
                    }
                    start_loop=end_loop+1;
                }
                free(para_each);
                free(t);
                free(p_hashvalueTable);
            }
            else
            {
                uint64_t current[4];
                current[0]=cal_hash_value_directly(seq,kmer_length);
                Transition('0',seq[kmer_length],(p_hash_table+current[0]));

                for(uint64_t i=1;i<seq_length-kmer_length;i++)
                {
                    current[0]=cal_hash_value_indirectly(seq+i,current[0],kmer_length);
                    Transition(seq[i-1],seq[i+kmer_length],(p_hash_table+current[0]));
                }
                current[0]=cal_hash_value_indirectly(seq+seq_length-kmer_length,current[0],kmer_length);
                Transition(seq[seq_length-kmer_length-1],'0',(p_hash_table+current[0]));
            }
        }
		free(seq);
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"HashTable construct time is: "<<span << endl;

	struct timeval tvs_check,tve_check;
	gettimeofday(&tvs_check,NULL);


	uint64_t sum_being=0,sum_in=0,sum_out=0;
	FILE *fp_in,*fp_out;

	char* outputFileName_in;
	outputFileName_in=(char*)malloc(sizeof(char)*11);
	struct para_getN tmp_outname;
	tmp_outname.kmerlen=kmer_len_being;
	tmp_outname.FileName=outputFileName_in;
	tmp_outname.InOut_label=0;
	tmp_outname.isMiddle=0;
	tmp_outname.isposition=0;
	getFileName(tmp_outname);
	fp_in=fopen(outputFileName_in,"wb+");

	char* outputFileName_out;
	outputFileName_out=(char*)malloc(sizeof(char)*11);
	tmp_outname.kmerlen=kmer_len_being;
	tmp_outname.FileName=outputFileName_out;
	tmp_outname.InOut_label=1;
	tmp_outname.isMiddle=0;
	tmp_outname.isposition=0;
	getFileName(tmp_outname);
	fp_out=fopen(outputFileName_out,"wb+");

	FILE *fp_in_ad,*fp_out_ad;
	char* outputFileName_ad_in;
    struct para_getN tmp_ad_in;
    if(ad_label==1)
    {
        outputFileName_ad_in=(char*)malloc(sizeof(char)*13);
        tmp_ad_in.kmerlen=kmer_len_being;
        tmp_ad_in.FileName=outputFileName_ad_in;
        tmp_ad_in.InOut_label=0;
        tmp_ad_in.isMiddle=2;
        tmp_ad_in.isposition=0;
        getFileName(tmp_ad_in);
        fp_in_ad=fopen(outputFileName_ad_in,"ab+");
        free(outputFileName_ad_in);
    }

    char* outputFileName_ad_out;
    struct para_getN tmp_ad_out;
    if(ad_label==1)
    {
        outputFileName_ad_out=(char*)malloc(sizeof(char)*13);
        tmp_ad_out.kmerlen=kmer_len_being;
        tmp_ad_out.FileName=outputFileName_ad_out;
        tmp_ad_out.InOut_label=0;
        tmp_ad_out.isMiddle=2;
        tmp_ad_out.isposition=0;
        getFileName(tmp_ad_out);
        fp_in_ad=fopen(outputFileName_ad_out,"ab+");
        free(outputFileName_ad_out);
    }


	cout <<"start checking hash table for branched k-mer."<< endl;

	if(thread_num>1)
	{
		uint64_t **reslt_in;
		uint64_t **reslt_out;
		reslt_in=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);
		reslt_out=(uint64_t **)malloc(sizeof(uint64_t *)*thread_num);

		uint8_t **reslt_in_ad;
		uint8_t **reslt_out_ad;

		if(ad_label==1)
        {
            reslt_in_ad=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);
            reslt_out_ad=(uint8_t **)malloc(sizeof(uint8_t *)*thread_num);
        }

		uint64_t * reslt_size_in;
		reslt_size_in=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);
		uint64_t * reslt_size_out;
		reslt_size_out=(uint64_t *)malloc(sizeof(uint64_t)*thread_num);

		struct reslt_para_being* tmp;
		tmp=(struct reslt_para_being*)malloc(sizeof(struct reslt_para_being)*thread_num);

		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		uint64_t unit_size=Task_Size;
		uint32_t loop_num=hash_table_length/unit_size+hash_table_length%unit_size;
		for(uint32_t i=0;i<loop_num;i++)
		{
			uint64_t start_loop=unit_size*i;
			uint64_t end_loop;
			if(i==loop_num-1)
			{
				end_loop=hash_table_length-1;
			}
			else
			{
				end_loop=unit_size*(i+1)-1;
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
				tmp[j].start=start_loop;
				tmp[j].task_len=end_loop-start_loop+1;
				tmp[j].thread_id=j;
				tmp[j].thread_num=thread_num;
				tmp[j].reslt_in=&(reslt_in[j]);
				tmp[j].reslt_out=&(reslt_out[j]);
				if(ad_label==1)
                {
                    tmp[j].reslt_in_ad=&(reslt_in_ad[j]);
                    tmp[j].reslt_out_ad=&(reslt_out_ad[j]);
                }
				tmp[j].reslt_size_in=&(reslt_size_in[j]);
				tmp[j].reslt_size_out=&(reslt_size_out[j]);
				tmp[j].p_hash_table=&p_hash_table;

				if(pthread_create(t+j, NULL, judge_thread_being, (void*)(tmp+j))!=0)
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
                fwrite(reslt_in[j],sizeof(uint64_t),reslt_size_in[j],fp_in);
                if(ad_label==1)
                {
                    fwrite(reslt_in_ad[j],sizeof(uint8_t),reslt_size_in[j],fp_in_ad);
                }
                sum_in+=reslt_size_in[j];
//				for(uint64_t k=0;k<reslt_size_in[j];k++)
//				{
//					fwrite(&(reslt_in[j][k]),sizeof(uint64_t),1,fp_in);
//					sum_in++;
//				}
			}
			for(uint32_t j=0;j<thread_num;j++)
			{
			    fwrite(reslt_out[j],sizeof(uint64_t),reslt_size_out[j],fp_out);
                if(ad_label==1)
                {
                    fwrite(reslt_out_ad[j],sizeof(uint8_t),reslt_size_out[j],fp_out_ad);
                }
                sum_out+=reslt_size_out[j];
//				for(uint64_t k=0;k<reslt_size_out[j];k++)
//				{
//					fwrite(&(reslt_out[j][k]),sizeof(uint64_t),1,fp_out);
//					sum_out++;
//				}
			}
		}

		free(t);
		free(reslt_size_in);
		free(reslt_size_out);
		for(uint32_t i=0;i<thread_num;i++)
		{
			free(reslt_in[i]);
			free(reslt_out[i]);
		}
		free(reslt_in);
		free(reslt_out);
		free(tmp);
	}
	else
	{
		for(uint64_t i=0;i<hash_table_length;i++)
		{
			uint32_t tempkind;
			tempkind=Judge(p_hash_table[i]);
			if(tempkind==1||tempkind==3){
				fwrite(&i,sizeof(uint64_t),1,fp_in);
				if(ad_label==1)
                {
                    uint8_t tmp;
                    tmp=(p_hash_table[i]>>4)&15;
                    fwrite(&tmp,sizeof(uint8_t),1,fp_in_ad);
                }
				sum_in++;
			}
			if(tempkind==2||tempkind==3)
			{
				fwrite(&i,sizeof(uint64_t),1,fp_out);
				if(ad_label==1)
                {
                    uint8_t tmp;
                    tmp=p_hash_table[i]&15;
                    fwrite(&tmp,sizeof(uint8_t),1,fp_out_ad);
                }
				sum_out++;
			}
			if(p_hash_table[i]!=0){
				sum_being++;
			}
		}
	}
	cout <<"the # of being k-mer is: "<<sum_being << endl;
	cout <<"the # of in-branched k-mer is: "<<sum_in << endl;
	cout <<"the # of out-branched k-mer is: "<<sum_out << endl;
	free(outputFileName_in);
	free(outputFileName_out);
	fclose(fp_in);
	fclose(fp_out);
	if(ad_label==1)
    {
        fclose(fp_in_ad);
        fclose(fp_out_ad);
    }

	free(p_hash_table);
	gettimeofday(&tve_check,NULL);
	double span_check = tve_check.tv_sec-tvs_check.tv_sec + (tve_check.tv_usec-tvs_check.tv_usec)/1000000.0;
	cout <<"check time is: "<<span_check<< endl;
	gettimeofday(&tve_total,NULL);
	double span_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total time is: "<<span_total <<endl;
}

