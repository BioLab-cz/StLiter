#include "Combine_Kmers.h"

void uion_Kmers_bit(uint64_t kmer_len,uint32_t ad_label)
{
	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);

	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);

	struct para_getN tmp;
	tmp.kmerlen=kmer_len;
	tmp.FileName=inputFileName;
	tmp.InOut_label=0;
	tmp.isMiddle=0;
	tmp.isposition=0;
	getFileName(tmp);

	//read in-branched kmer file
	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	uint64_t *h0;
	uint64_t x;
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
	h0=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*para.kmer64Len);
	fseek(int_input,0,0);
	x=fread(h0,sizeof(uint64_t),total_hash_number*para.kmer64Len,int_input);

	FILE* int_input_ad;
	uint8_t *h0_ad;
	if(ad_label==1)
	{
		h0_ad=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number);
		tmp.FileName=inputFileName;
		tmp.InOut_label=0;
		tmp.isMiddle=2;
		getFileName(tmp);
		int_input_ad=fopen(inputFileName,"rb");
		x=fread(h0_ad,sizeof(uint8_t),total_hash_number,int_input_ad);
		for(uint64_t i=0;i<total_hash_number;i++)
		{
			h0_ad[i]=h0_ad[i]<<4;
		}
	}

	//read out-branched kmer file
	tmp.kmerlen=kmer_len;
	tmp.FileName=inputFileName;
	tmp.InOut_label=1;
	tmp.isMiddle=0;
	tmp.isposition=0;
	tmp.isposition=0;
	getFileName(tmp);

	FILE* int_input1;
	int_input1=fopen(inputFileName,"rb");
	uint64_t total_hash_number1=0;
	fseek(int_input1,0,2);
	uint64_t *h1;
	uint64_t x1;
	total_hash_number1=ftell(int_input1)/(sizeof(uint64_t)*para.kmer64Len);
	h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number1*para.kmer64Len);
	fseek(int_input1,0,0);
	x1=fread(h1,sizeof(uint64_t),total_hash_number1*para.kmer64Len,int_input1);

	FILE* int_input1_ad;
	uint8_t *h1_ad;
	if(ad_label==1)
	{
		h1_ad=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number1);
		tmp.FileName=inputFileName;
		tmp.InOut_label=1;
		tmp.isMiddle=2;
		tmp.isposition=0;
		getFileName(tmp);
		int_input1_ad=fopen(inputFileName,"rb");
		x=fread(h1_ad,sizeof(uint8_t),total_hash_number,int_input1_ad);
	}
	cout << "read files over!" << endl;
	fclose(int_input);
	fclose(int_input1);
	if(ad_label==1)
	{
		fclose(int_input_ad);
		fclose(int_input1_ad);
	}

	uint64_t nin=0;
	uint64_t nout=0;
	tmp.kmerlen=kmer_len;
	tmp.FileName=inputFileName;
	tmp.InOut_label=1;
	tmp.isMiddle=1;
	getFileName(tmp);
	tmp.isposition=0;
	FILE* int_input2;
	int_input2=fopen(inputFileName,"wb");

	tmp.isMiddle=3;
	tmp.isposition=0;
	getFileName(tmp);
	FILE* int_input2_ad;
	if(ad_label==1)
	{
		int_input2_ad=fopen(inputFileName,"wb");
	}
	while(1)
	{
		if(nin==total_hash_number&&nout==total_hash_number1)
		{
			break;
		}
		else if(nin==total_hash_number)
		{
			fwrite(h1+nout*para.kmer64Len,sizeof(uint64_t),para.kmer64Len*(total_hash_number1-nout),int_input2);
			if(ad_label==1)
			{
				fwrite(h1_ad+nout,sizeof(uint8_t),(total_hash_number1-nout),int_input2_ad);
			}
			break;
		}
		else if(nout==total_hash_number1)
		{
			fwrite(h0+nin*para.kmer64Len,sizeof(uint64_t),para.kmer64Len*(total_hash_number-nin),int_input2);
			if(ad_label==1)
			{
				fwrite(h0_ad+nin,sizeof(uint8_t),(total_hash_number-nin),int_input2_ad);
			}
			break;
		}
		else
		{
			if(ad_label==1)
			{
				if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==1)
				{
					fwrite(h1+nout*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					fwrite(h1_ad+nout,sizeof(uint8_t),1,int_input2_ad);
					nout++;
				}
				else if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==0)
				{
					fwrite(h0+nin*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					fwrite(h0_ad+nin,sizeof(uint8_t),1,int_input2_ad);
					nin++;
				}
				else if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==2)
				{
					fwrite(h0+nin*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					uint8_t tmp=0;
					tmp=h0_ad[nin];
					tmp=tmp|h1_ad[nout];
					fwrite(&tmp,sizeof(uint8_t),1,int_input2_ad);
					nin++;
					nout++;
				}
			}
			else
			{
				if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==1)
				{
					fwrite(h1+nout*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					nout++;
				}
				else if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==0)
				{
					fwrite(h0+nin*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					nin++;
				}
				else if(cmp256BitKmer(h0+nin*para.kmer64Len,h1+nout*para.kmer64Len,para.kmer64Len)==2)
				{
					fwrite(h0+nin*para.kmer64Len,sizeof(uint64_t),para.kmer64Len,int_input2);
					nin++;
					nout++;
				}
			}
		}
	}
	fclose(int_input2);
	if(ad_label==1)
	{
		fclose(int_input2_ad);
	}
	free(inputFileName);

	cout << "done!" << endl;
}
