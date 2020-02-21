/*
 * Fast_c_dbg.cpp
 *
 *  Created on: Apr 14, 2018
 *      Author: bio
 */

#include "Basic.h"

void getFileName(struct para_getN tmp)
{
	uint32_t kmerlen=tmp.kmerlen;
	char* FileName=tmp.FileName;
	uint32_t label=tmp.InOut_label;
	uint32_t isMiddle=tmp.isMiddle;
	uint32_t isposition=tmp.isposition;

	if(isposition==0)
	{
		char init[5]="kmer";
		strcpy(FileName,init);
	}
	else if(isposition==10000)
	{
		char init[5]="unip";
		strcpy(FileName,init);
	}
	else
	{
		char init[5]="posi";
		strcpy(FileName,init);
	}
	char swap[4];
	sprintf(swap,"%d",kmerlen);

	char in_add[3]="in";
	char out_add[4]="out";
	char ad[3]="ad";

	if(kmerlen<100)
	{
		FileName[4]='0';
		strcpy(FileName+5,swap);
		if(isMiddle==1)
		{
			FileName[7]='\0';
		}
		else if(isMiddle==3)
		{
			strcpy(FileName+7,ad);
			FileName[9]='\0';
		}
		else
		{
			if(label==0)
			{
				strcpy(FileName+7,in_add);
				if(isMiddle==2)
				{
					strcpy(FileName+9,ad);
					FileName[11]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[9]='0';
						FileName[10]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[12]='\0';
					}
					else if(isposition<101)
					{
						FileName[9]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[12]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+9,swap);
						FileName[12]='\0';
					}
				}
				else
				{
					FileName[9]='\0';
				}
			}
			else
			{
				strcpy(FileName+7,out_add);
				if(isMiddle==2)
				{
					strcpy(FileName+10,ad);
					FileName[12]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[10]='0';
						FileName[11]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+12,swap);
						FileName[13]='\0';
					}
					else if(isposition<101)
					{
						FileName[10]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[13]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[13]='\0';
					}
				}
				else
				{
					FileName[10]='\0';
				}
			}
		}
	}
	else
	{
		strcpy(FileName+4,swap);
		if(isMiddle==1)
		{
			FileName[7]='\0';
		}
		else if(isMiddle==3)
		{
			strcpy(FileName+7,ad);
			FileName[9]='\0';
		}
		else
		{
			if(label==0)
			{
				strcpy(FileName+7,in_add);
				if(isMiddle==2)
				{
					strcpy(FileName+9,ad);
					FileName[11]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[9]='0';
						FileName[10]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[12]='\0';
					}
					else if(isposition<101)
					{
						FileName[9]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[12]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+9,swap);
						FileName[12]='\0';
					}
				}
				else
				{
					FileName[9]='\0';
				}
			}
			else
			{
				strcpy(FileName+7,out_add);
				if(isMiddle==2)
				{
					strcpy(FileName+10,ad);
					FileName[12]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[10]='0';
						FileName[11]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+12,swap);
						FileName[13]='\0';
					}
					else if(isposition<101)
					{
						FileName[10]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[13]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[13]='\0';
					}
				}
				else
				{
					FileName[10]='\0';
				}
			}
		}
	}
}

void ReadSeq(char **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t buffer_size=256;
	char buffer_line[256];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
		return;
	}

	uint64_t total_size=0;
	fseek(fp,0,2);
	total_size=ftell(fp);

	char *seq;
	seq=(char*) malloc (sizeof(char)*total_size);

	fseek(fp,0,0);

	uint64_t len=0;
	while (fgets(buffer_line,buffer_size-1,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			for(uint32_t i=0;i<buffer_size;i++)
			{
				if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
				{
					break;
				}
				if(buffer_line[i]>='a')
				{
					buffer_line[i]-=32;
				}
				seq[len]=buffer_line[i];
				len++;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len;
	*seq1=seq;
	cout << "the length of seq is: " << len << endl;
}

void rc(char **seq_or, uint64_t seq_length)
{
    char* seq;
    seq=*seq_or;
    char* rcseq;
    char* tmp;
    rcseq=(char*)malloc(sizeof(char)*seq_length);

    tmp=seq+seq_length-1;
    for(uint64_t i=0;i<seq_length;i++)
    {
        switch(*tmp)
		{
			case 'A':
				rcseq[i]='T';
				break;
			case 'C':
				rcseq[i]='G';
				break;
			case 'G':
				rcseq[i]='C';
				break;
			case 'T':
				rcseq[i]='A';
				break;
			default:
				rcseq[i]='T';
				break;
		}
		tmp--;
    }

    free(seq);
    *seq_or=rcseq;
}

void ReadSeq_identical(char* p_ref)
{
	FILE * fp_out;
	fp_out=fopen("a.fa","a+");

	uint32_t buffer_size=256;
	char buffer_line[256];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint64_t len=0;
	uint32_t cycle=0;
	while (fgets(buffer_line,buffer_size-1,fp)!=NULL)
	{
		for(uint32_t i=0;i<buffer_size;i++)
		{
			fprintf(fp_out,"%c",buffer_line[i]);
			if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
			{
				fprintf(fp_out,"\n");
				break;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	cout << "the length of seq is: " << len << endl;
}
void getRefFilePathes(char* pathFile, struct RefFilePath* p)
{
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	ifstream int_input;

	cout << pathFile << endl;
	int_input.open(pathFile,ios::in);
	uint64_t total_hash_number=0;
	char hv_tmp[64];
	while(int_input.getline(hv_tmp,63))
	{
		total_hash_number++;
	}
	int_input.close();
	cout << "total ref number is: " <<total_hash_number << endl;

	p->pRefFilePath=(char**)malloc(sizeof(char*)*total_hash_number);
	for(uint32_t i=0;i<total_hash_number;i++)
	{
		p->pRefFilePath[i]=(char*)malloc(sizeof(char)*256);
	}

	uint32_t line_ref_path=0;
	int_input.open(pathFile,ios::in);
	while(int_input.getline(hv_tmp,63))
	{
		strcpy(p->pRefFilePath[line_ref_path],hv_tmp);
		line_ref_path++;
	}
	int_input.close();
	p->NumberOfPathes=line_ref_path;

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
}

/***************************ForBitCharKmer*********************************/
void fill_char_with_four_char(uint8_t* current,char*p)
{
	for(uint32_t i=0;i<3;i++)
	{
		switch(p[i])
		{
			case 'A':
				*current=*current<<2;
				break;
			case 'C':
				*current=*current|1;
				*current=*current<<2;
				break;
			case 'G':
				*current=*current|2;
				*current=*current<<2;
				break;
			case 'T':
				*current=*current|3;
				*current=*current<<2;
				break;
			default:
				*current=*current<<2;
				break;
		}
	}
	switch(p[3])
	{
		case 'A':
			break;
		case 'C':
			*current=*current|1;
			break;
		case 'G':
			*current=*current|2;
			break;
		case 'T':
			*current=*current|3;
			break;
		default:
			break;
	}
}
void ReadSeq_bit(uint8_t **seq1,uint64_t *seq_length,char* p_ref)
{
	uint64_t cursize;
	uint64_t maxsize;
	uint64_t addsize;

	maxsize=pow(2,29);
	addsize=pow(2,29);

	uint8_t *seq;//=seq1;
	seq=(uint8_t*) malloc (sizeof(uint8_t)*maxsize);
	cursize=maxsize;


	uint64_t * p;
	p=(uint64_t *)malloc(sizeof(uint64_t)*5);


	uint32_t buffer_size=_BufferSize_;
	char buffer_line[_BufferSize_];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint64_t len=0;
	uint32_t cycle=0;
	uint32_t number_left=0;
	uint32_t number_buffer_char;
	char array_left[4];
	while (1)
	{
		if(fgets(buffer_line,buffer_size-number_left,fp)==NULL)
		{
			break;
		}
		if(buffer_line[0]=='>')
		{
			continue;
		}
		else
		{
//			cout << buffer_line << endl;
//			cout << strlen(buffer_line) << endl;
			char buffer_line_tmp[_BufferSize_];
			if(number_left!=0)
			{
				strcpy(buffer_line_tmp,array_left);
				strcpy(buffer_line_tmp+number_left,buffer_line);
			}
			else
			{
				strcpy(buffer_line_tmp,buffer_line);
			}

			uint32_t buffer_line_char_number=strlen(buffer_line_tmp);
			if(buffer_line_tmp[buffer_line_char_number-1]=='\n')
			{
				buffer_line_char_number--;
			}
//			cout <<buffer_line_char_number << endl;
			number_left=buffer_line_char_number%4;
			number_buffer_char=buffer_line_char_number/4;
			if(number_left!=0)
			{
				for(uint32_t i=0;i<number_left;i++)
				{
					array_left[i]=buffer_line_tmp[buffer_line_char_number-(number_left-i)];
				}
				array_left[number_left]='\0';
			}

			if(len+number_buffer_char<cursize)
			{
				for(uint32_t i=0;i<number_buffer_char;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
			}
			else
			{
				seq=(uint8_t*) realloc (seq,sizeof(uint8_t)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_line_char_number/4;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
				cout <<"add 1024*1024*1024 byte for seq: " << cycle++ <<endl;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len*4+number_left;
	if(number_left!=0)
	{
		for(uint32_t i=number_left;i<4;i++)
		array_left[i]='A';
	}
	uint8_t tmp;
	fill_char_with_four_char(&tmp,array_left);
	seq[len]=tmp;
	len++;
	seq1[0]=seq;

	uint8_t* seq_shift1=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift1[i]=seq[i]<<2;
		seq_shift1[i]=seq_shift1[i]|(seq[i+1]>>6);
	}
	seq_shift1[len-1]=seq[len-1]<<2;
	seq1[1]=seq_shift1;

	uint8_t* seq_shift2=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift2[i]=seq_shift1[i]<<2;
		seq_shift2[i]=seq_shift2[i]|(seq_shift1[i+1]>>6);
	}
	seq_shift2[len-1]=seq_shift1[len-1]<<2;
	seq1[2]=seq_shift2;

	uint8_t* seq_shift3=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift3[i]=seq_shift2[i]<<2;
		seq_shift3[i]=seq_shift3[i]|(seq_shift2[i+1]>>6);
	}
	seq_shift3[len-1]=seq_shift2[len-1]<<2;
	seq1[3]=seq_shift3;

	cout << "the length of seq is: " << *seq_length << endl;
}
void Fetch_kmer_on_bit_sequence(uint8_t **seq,uint64_t start, struct bitKmerPara para,uint8_t *kmer)
{
	uint64_t s_id=start/4;
	uint32_t s_shift=start%4;
	for(uint32_t i=0;i<para.kmer8Len;i++)
	{
		kmer[i]=seq[s_shift][s_id+i];
	}
}
void calMulDimHashValue(uint8_t *kmer,uint64_t *hv64,struct bitKmerPara para)
{
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		hv64[i]=0;
		if(i==para.kmer64Len-1)
		{
			for(uint32_t j=0;j<para.remainer8to64-1;j++)
			{
				hv64[i]=hv64[i]|kmer[8*i+j];
				hv64[i]=hv64[i]<<8;
			}
			hv64[i]=hv64[i]|kmer[8*i+para.remainer8to64-1];
			if(para.remainer1to8!=0)
			{
				hv64[i]=hv64[i]>>(8-para.remainer1to8);
			}
		}
		else
		{
			for(uint32_t j=0;j<7;j++)
			{
				hv64[i]=hv64[i]|kmer[8*i+j];
				hv64[i]=hv64[i]<<8;
			}
			hv64[i]=hv64[i]|kmer[8*i+7];
		}
	}
}
void Transition_bitchar(uint8_t start,uint8_t end,uint8_t* temp_res){
	uint8_t temp_b;
	uint8_t *res;
	res=temp_res;
	switch((uint32_t)start)
	{
		case 0:
			temp_b=16;
			break;
		case 1:
			temp_b=32;
			break;
		case 2:
			temp_b=64;
			break;
		case 3:
			temp_b=128;
			break;
		case 4:
			temp_b=0;
			break;
		default:
			temp_b=0;
			break;
	}
	(*res)=(*res)|temp_b;
	switch((uint32_t)end)
	{
		case 0:
			temp_b=1;
			break;
		case 1:
			temp_b=2;
			break;
		case 2:
			temp_b=4;
			break;
		case 3:
			temp_b=8;
			break;
		case 4:
			temp_b=0;
			break;
		default:
			temp_b=0;
			break;
	}
	(*res)=(*res)|temp_b;
}
/***************************ForBitCharKmer*********************************/

/***************************For256BitCharKmer******************************/
void cal_hash_value_directly_256bit(char *seq,uint64_t * current,\
		struct bit256KmerPara para)
{
	uint64_t tmp;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		char* loop_tmp=k_mer_temp+32*i;
		if(i==para.kmer64Len-1&&para.remainer1to64!=0)
		{
			tmp=0;
			for(uint32_t j=0;j<para.remainer1to64/2-1;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[para.remainer1to64/2-1])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
		else
		{
			tmp=0;
			for(uint32_t j=0;j<31;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[31])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
	}
}
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para)
{
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len-1;i++)
	{
		if(i==para.kmer64Len-2&&para.remainer1to64!=0)
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>(para.remainer1to64-2));
		}
		else
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>62);
		}
	}
	if(para.remainer1to64==0)
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
	else
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		current[para.kmer64Len-1]=current[para.kmer64Len-1]&para.codefor1;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
}
void cal_newExtendedKmerIn(uint64_t* original,uint64_t* newKmer,uint64_t j,uint32_t exLen,struct bit256KmerPara para)
{
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		newKmer[i]=original[i];
	}
	uint32_t exBitLen=2*exLen;
	if(para.remainer1to64==0)
	{
		newKmer[para.kmer64Len]=j;
	}
	else
	{
		if(exBitLen+para.remainer1to64<=64)
		{
			newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]<<exBitLen;
			newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|j;
		}
		else
		{
			newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]<<(64-para.remainer1to64);
			newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|(j>>(exBitLen-(64-para.remainer1to64)));
			newKmer[para.kmer64Len]=j<<(128-(exBitLen+para.remainer1to64));
			newKmer[para.kmer64Len]=newKmer[para.kmer64Len]>>(128-(exBitLen+para.remainer1to64));
		}
	}
}
void cal_newExtendedKmerOut(uint64_t* original,uint64_t* newKmer,uint64_t j,uint32_t exLen,struct bit256KmerPara para)
{
	uint32_t exBitLen=2*exLen;
	if(para.remainer1to64==0)
	{
		for(uint32_t i=0;i<para.kmer64Len;i++)
		{
			newKmer[i]=0;
			newKmer[i]=original[i]>>exBitLen;
		}
		newKmer[0]=newKmer[0]|(j<<(64-exBitLen));
		for(uint32_t i=1;i<para.kmer64Len;i++)
		{
			newKmer[i]=newKmer[i]|(original[i-1]<<(64-exBitLen));
		}
		newKmer[para.kmer64Len]=original[para.kmer64Len-1];
		newKmer[para.kmer64Len]=(newKmer[para.kmer64Len]<<(64-exBitLen))>>((64-exBitLen));
	}
	else
	{
		if(para.kmer64Len>1)
		{
			for(uint32_t i=0;i<para.kmer64Len-1;i++)
			{
				newKmer[i]=0;
				newKmer[i]=original[i]>>exBitLen;
			}
			newKmer[0]=newKmer[0]|(j<<(64-exBitLen));
			for(uint32_t i=1;i<para.kmer64Len-1;i++)
			{
				newKmer[i]=newKmer[i]|(original[i-1]<<(64-exBitLen));
			}
			if(exBitLen+para.remainer1to64<=64)
			{
				newKmer[para.kmer64Len-1]=0;
				newKmer[para.kmer64Len-1]=original[para.kmer64Len-1];
				newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|((original[para.kmer64Len-2]<<(64-exBitLen))>>(64-exBitLen-para.remainer1to64));
			}
			else
			{
				newKmer[para.kmer64Len-1]=0;
				newKmer[para.kmer64Len-1]=original[para.kmer64Len-1]>>(exBitLen+para.remainer1to64-64);
				newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|(original[para.kmer64Len-2]<<(64-exBitLen));

				newKmer[para.kmer64Len]=original[para.kmer64Len-1];
				newKmer[para.kmer64Len]=newKmer[para.kmer64Len]<<(128-exBitLen-para.remainer1to64);
				newKmer[para.kmer64Len]=newKmer[para.kmer64Len]>>(128-exBitLen-para.remainer1to64);
			}
		}
		else
		{
			if(exBitLen+para.remainer1to64<=64)
			{
				newKmer[para.kmer64Len-1]=0;
				newKmer[para.kmer64Len-1]=original[para.kmer64Len-1];
				newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|(j<<para.remainer1to64);
			}
			else
			{
				newKmer[para.kmer64Len-1]=0;
				newKmer[para.kmer64Len-1]=original[para.kmer64Len-1]>>(exBitLen+para.remainer1to64-64);
				newKmer[para.kmer64Len-1]=newKmer[para.kmer64Len-1]|(j<<(64-exBitLen));

				newKmer[para.kmer64Len]=original[para.kmer64Len-1];
				newKmer[para.kmer64Len]=newKmer[para.kmer64Len]<<(128-exBitLen-para.remainer1to64);
				newKmer[para.kmer64Len]=newKmer[para.kmer64Len]>>(128-exBitLen-para.remainer1to64);
			}
		}
	}

}
uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else
		{
			if(a[i]>b[i])
			{
				r=1;
				break;
			}
		}
	}
	return r;
}
/***************************For256BitCharKmer******************************/

void Transition(char start,char end,uint8_t* temp_res){
	uint8_t temp_b;
	uint8_t *res;
	res=temp_res;
	switch(start)
	{
		case 'A':
			temp_b=16;
			break;
		case 'C':
			temp_b=32;
			break;
		case 'G':
			temp_b=64;
			break;
		case 'T':
			temp_b=128;
			break;
		case '0':
			temp_b=0;
			break;
		default:
			temp_b=16;
			break;
	}
	(*res)=(*res)|temp_b;
	switch(end)
	{
		case 'A':
			temp_b=1;
			break;
		case 'C':
			temp_b=2;
			break;
		case 'G':
			temp_b=4;
			break;
		case 'T':
			temp_b=8;
			break;
		case '0':
			temp_b=0;
			break;
		default:
			temp_b=1;
			break;
	}
	(*res)=(*res)|temp_b;
}
uint32_t Judge(uint8_t temp_a){
	uint32_t res;
	uint8_t temp=15;
	temp=temp_a&temp;
	if(!(temp==0||temp==1||temp==2||temp==4||temp==8)){
		res=2;
	}
	else{
		res=0;
	}
	temp=temp_a>>4;
	temp=temp&15;
	if(!(temp==0||temp==1||temp==2||temp==4||temp==8)){
		if(res==0){
			res=1;
		}
		else{
			res=3;
		}
	}
	return res;
}
uint32_t Judge_novel(uint8_t temp_a){
	if(temp_a == 0)
	{
		return 4;
	}
	uint32_t res;
	uint8_t temp=15;
	temp=temp_a&temp;
	if(!(temp==0||temp==1||temp==2||temp==4||temp==8)){
		res=2;
	}
	else{
		res=0;
	}
	temp=temp_a>>4;
	if(!(temp==0||temp==1||temp==2||temp==4||temp==8)){
		if(res==0){
			res=1;
		}
		else{
			res=3;
		}
	}

	return res;
}

void Transition_4bit(char c_ad,uint8_t* temp_res,uint32_t label_4bit)
{
	//label_4bit=0,deal with the higher 4bit, otherwise the lower 4bit.
	uint8_t temp_b;
	uint8_t *res;
	res=temp_res;
	switch(c_ad)
	{
		case 'A':
			temp_b=1;
			break;
		case 'C':
			temp_b=2;
			break;
		case 'G':
			temp_b=4;
			break;
		case 'T':
			temp_b=8;
			break;
		case '0':
			temp_b=0;
			break;
		default:
			temp_b=1;
			break;
	}

	if(label_4bit==0)
	{
		temp_b=temp_b<<4;
		(*res)=(*res)|temp_b;
	}
	else
	{
		(*res)=(*res)|temp_b;
	}
}
uint32_t Judge_4bit(uint8_t temp_a,uint32_t label_4bit)
{
	uint32_t res;
	char temp;
	if(label_4bit==0)
	{
		temp=temp_a>>4;
		temp=temp&15;
		if(!(temp==0||temp==1||temp==2||temp==4||temp==8))
		{
			res=1;
		}
		else
		{
			res=0;
		}
	}
	else
	{
		temp=temp_a&15;
		if(!(temp==0||temp==1||temp==2||temp==4||temp==8))
		{
			res=1;
		}
		else{
			res=0;
		}
	}
	return res;
}

uint64_t cal_hash_value_directly(char *seq,uint32_t len)
{
	uint64_t current=0;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<len;i++)
	{

		switch(k_mer_temp[i])
		{
			case 'A':
				current=current<<2;
				break;
			case 'C':
				current=current|1;
				current=current<<2;
				break;
			case 'G':
				current=current|2;
				current=current<<2;
				break;
			case 'T':
				current=current|3;
				current=current<<2;
				break;
			default:
				current=current<<2;
				break;
		}
	}
	current=current>>2;
	return current;
}
uint64_t cal_hash_value_indirectly(char *seq,uint64_t current,uint32_t len)
{

	char *k_mer_temp=seq;
	uint64_t high=0;
	for(uint32_t i=0;i<len*2-3;i++){
		high=high|1;
		high=high<<1;
	}
	high=high|1;
//	high=high<<(HashSize-2);
	current=current&high;//0000 0000 0000 0011 1111 1111 1111 1111
	current=current<<2;
	switch(k_mer_temp[len-1])
	{
		case 'A':
			//current=current|0;
			break;
		case 'C':
			current=current|1;
			//current=current<<2;
			break;
		case 'G':
			current=current|2;
			//current=current<<2;
			break;
		case 'T':
			current=current|3 ;
			//current=current<<2;
			break;
		default:
			break;
	}
	return current;
}
