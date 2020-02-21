/*
 * litool.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: bio
 */
#include "litool.h"
void transKmerBinaryToText(char *p1,char*p2,uint32_t n)
{
	FILE* fp_in;
	FILE* fp_out;

	fp_in=fopen(p1,"rb");
	fp_out=fopen(p2,"wb+");

	uint64_t x[4];
	while(fread(&x,sizeof(uint64_t),n,fp_in)==n)
	{
		for(uint32_t i=0;i<n-1;i++)
		{
			fprintf(fp_out,"%llu	",x[i]);
		}
		fprintf(fp_out,"%llu",x[n-1]);
		fprintf(fp_out,"\n");
	}
	fclose(fp_in);
	fclose(fp_out);
}
void transCountKmerBinaryToText(char *p1,char*p2,uint32_t n)
{
	FILE* fp_in;
	FILE* fp_out;

	fp_in=fopen(p1,"rb");
	fp_out=fopen(p2,"wb+");

	uint64_t x[5];
	while(fread(&x,sizeof(uint64_t),n+1,fp_in)==n+1)
	{
		for(uint32_t i=0;i<n;i++)
		{
			fprintf(fp_out,"%llu	",x[i]);
		}
		fprintf(fp_out,"%llu",x[n]);
		fprintf(fp_out,"\n");
	}
	fclose(fp_in);
	fclose(fp_out);
}
void CheckResult_num(char *p1)
{
	ifstream f1;
	f1.open(p1,ios::in);
	uint64_t n1=0;

	char buffer1[64];
	while(f1.getline(buffer1,63))
	{
		n1++;
	}
	cout << n1 << endl;

	f1.close();
}
void CheckResult_content(char *p1,char *p2)
{
	ifstream f1;
	ifstream f2;
	f1.open(p1,ios::in);
	f2.open(p2,ios::in);
	uint64_t n1=0;
	uint64_t n2=0;

	char buffer1[64];
	char buffer2[64];
	while(f1.getline(buffer1,63)&&f2.getline(buffer2,63))
	{
		uint64_t x1;
		uint64_t x2;
		x1=atoi(buffer1);
		x2=atoi(buffer2);
		if(x1!=x2)
		{
			cout << "error: a!=b!" << endl;
			return ;
		}
		n1++;
		n2++;
	}
	f1.close();
	f2.close();

}
void Check_bitBybit(char *p1,char *p2)
{
	FILE * f1;
	FILE * f2;
	f1=fopen(p1,"rb");
	f2=fopen(p2,"rb");

	uint64_t n1;
	uint64_t n2;

	fseek(f1,0,2);
	n1=ftell(f1);
	fseek(f2,0,2);
	n2=ftell(f2);
	cout << n1 << "	" << n2 <<endl;
	if(n1==n2)
	{
		fseek(f1,0,0);
		fseek(f2,0,0);
		uint8_t * cp1=(uint8_t *)malloc(sizeof(uint8_t)*n1);
		fread(cp1,sizeof(uint8_t),n1,f1);
		uint8_t * cp2=(uint8_t *)malloc(sizeof(uint8_t)*n2);
		fread(cp2,sizeof(uint8_t),n2,f2);
		char x;
		for(uint64_t i=0;i<n1;i++)
		{
			if(cp1[i]!=cp2[i])
			{
				cout << "content not identical!" << endl;
				cout << (uint32_t)cp1[i] << " " <<(uint32_t)cp2[i] << endl;
				cin >> x;
			}
		}
	}
	else
	{
		cout << "size not identical!" <<endl;
	}
}
void SortOutKmersUsingBplusTreebinarytobinary(char*inputFileName,char*outputFileName,uint32_t len)
{
	struct bit256KmerPara para;

	para.kmer64Len=len;

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*len);
	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*len);

	fseek(int_input,0,0);

	uint64_t x;
	x=fread(h1,sizeof(uint64_t),total_hash_number*len,int_input);
	cout << x << endl;
	cout << total_hash_number << endl;

	struct NodeBit * p_root;
	p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	struct nodeBit c_tmp;
	c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*len);
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		for(uint64_t j=0;j<len;j++)
		{
			c_tmp.hashValue[j]=h1[len*i+j];
		}
		c_tmp.arrayID=i+1;
		Insert_Value_bit(&p_root,c_tmp,para);
	}
	fclose(int_input);
	cout << "Constructing BplusTree is over!" <<endl;

	//mapping the first hash value to the first dimensional coordinate of HashAd
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		uint64_t ID=MappingHashValueToID_bit(p_root,h1+len*i,para);
		if((i+1)!=ID)
		{
			cout << i << endl;
			cout << h1[len*i] << "	" << h1[len*i+1] << endl;
			cout << "error bplustree!" << endl;
			char y;
			cin >> y;
		}
	}

	struct NodeBit * p_leaf=Find_Minimal_Node_bit(p_root);

	FILE* outputFile;
	outputFile=fopen(outputFileName,"wb");

	struct NodeBit *p_tmp=p_leaf;
	while(p_tmp!=NULL)
	{
		for(uint32_t i=0;i<p_tmp->Node_Size;i++)
		{
			fwrite(p_tmp->data[i].hashValue,sizeof(uint64_t),len,outputFile);
		}
		p_tmp=p_tmp->brother;
	}
	fclose(outputFile);
}
void SortKmers(char*inputFileName,char*outputFileName,char*inputFileName_ad,char*outputFileName_ad,uint32_t len)
{
	struct bit256KmerPara para;

	para.kmer64Len=len;

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*len);
	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*len);

	fseek(int_input,0,0);

	uint64_t x;
	x=fread(h1,sizeof(uint64_t),total_hash_number*len,int_input);
	cout << x << endl;
	cout << total_hash_number << endl;

	FILE* int_input_ad;
	uint8_t *h1_ad;
	if(inputFileName_ad!=NULL)
	{
		int_input_ad=fopen(inputFileName_ad,"rb");
		h1_ad=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number);
		x=fread(h1_ad,sizeof(uint8_t),total_hash_number,int_input_ad);
		cout << x << endl;
		cout << total_hash_number << endl;
	}

	struct NodeBit * p_root;
	p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	struct nodeBit c_tmp;
	c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*len);
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		for(uint64_t j=0;j<len;j++)
		{
			c_tmp.hashValue[j]=h1[len*i+j];
		}
		if(inputFileName_ad!=NULL)
		{
			c_tmp.ad=h1_ad[i];
		}

		c_tmp.arrayID=i+1;
//		Insert_Value_bit(&p_root,c_tmp,para);
		Insert_ad_256bitKmer(&p_root,c_tmp,para);
	}
	fclose(int_input);
	if(inputFileName_ad!=NULL)
	{
		fclose(int_input_ad);
	}
	cout << "Constructing BplusTree is over!" <<endl;

	//mapping the first hash value to the first dimensional coordinate of HashAd
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		uint64_t ID=MappingHashValueToID_bit(p_root,h1+len*i,para);
		if((i+1)!=ID)
		{
			cout << i << endl;
			cout << h1[len*i] << "	" << h1[len*i+1] << endl;
			cout << "error bplustree!" << endl;
			char y;
			cin >> y;
		}
	}

	struct NodeBit * p_leaf=Find_Minimal_Node_bit(p_root);

	FILE* outputFile;
	outputFile=fopen(outputFileName,"wb");

	FILE* outputFile_ad;
	if(inputFileName_ad!=NULL)
	{
		outputFile_ad=fopen(outputFileName_ad,"wb");
	}

	struct NodeBit *p_tmp=p_leaf;
	while(p_tmp!=NULL)
	{
		for(uint32_t i=0;i<p_tmp->Node_Size;i++)
		{
			fwrite(p_tmp->data[i].hashValue,sizeof(uint64_t),len,outputFile);
			if(inputFileName_ad!=NULL)
			{
				fwrite(&(p_tmp->data[i].ad),sizeof(uint8_t),1,outputFile_ad);
			}
		}
		p_tmp=p_tmp->brother;
	}
	fclose(outputFile);
	if(inputFileName_ad!=NULL)
	{
		fclose(outputFile_ad);
	}
}
void SortOutCountKmersUsingBplusTreebinarytobinary(char*inputFileName,char*outputFileName,uint32_t len)
{
	struct bit256KmerPara para;
	para.kmer64Len=len;

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*(len+1));
	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*(len+1));

	fseek(int_input,0,0);

	uint64_t x;
	x=fread(h1,sizeof(uint64_t),total_hash_number*(len+1),int_input);

	struct NodeBit * p_root;
	p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	struct nodeBit c_tmp;
	c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*len);
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		for(uint64_t j=0;j<len;j++)
		{
			c_tmp.hashValue[j]=h1[(len+1)*i+j];
		}
		c_tmp.arrayID=h1[(i+1)*(len+1)-1];
		Insert_Value_bit(&p_root,c_tmp,para);
	}
	fclose(int_input);
	cout << "Constructing BplusTree is over!" <<endl;

	//mapping the first hash value to the first dimensional coordinate of HashAd
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		uint64_t ID=MappingHashValueToID_bit(p_root,h1+i*(len+1),para);
		if((h1[(i+1)*(len+1)-1])!=ID)
		{
			cout << "error bplustree!" << endl;
			char y;
			cin >> y;
		}
	}

	struct NodeBit * p_leaf=Find_Minimal_Node_bit(p_root);

	FILE* outputFile;
	outputFile=fopen(outputFileName,"wb");

	struct NodeBit *p_tmp=p_leaf;
	while(p_tmp!=NULL)
	{
		for(uint32_t i=0;i<p_tmp->Node_Size;i++)
		{
			fwrite(p_tmp->data[i].hashValue,sizeof(uint64_t),len,outputFile);
			fwrite(&(p_tmp->data[i].arrayID),sizeof(uint64_t),1,outputFile);
		}
		p_tmp=p_tmp->brother;
	}
	fclose(outputFile);
	destory_tree_bit(p_root,para);
}
void SortOutCountKmersSingle(char*inputFileName,char*outputFileName,uint32_t len)
{
	struct bit256KmerPara para;
	para.kmer64Len=len;

	FILE* int_input;
	int_input=fopen(inputFileName,"rb");
	uint64_t total_hash_number=0;
	fseek(int_input,0,2);
	total_hash_number=ftell(int_input)/(sizeof(uint64_t)*(len+1));
	uint64_t *h1=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number*(len+1));

	fseek(int_input,0,0);

	uint64_t x;
	x=fread(h1,sizeof(uint64_t),total_hash_number*(len+1),int_input);
	if(x!=total_hash_number*(len+1))
	{
		cout << "error: kmer number wrong!" << endl;
	}
	remove(inputFileName);

	struct NodeBit * p_root;
	p_root=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	struct nodeBit c_tmp;
	c_tmp.hashValue=(uint64_t *)malloc(sizeof(uint64_t)*len);
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		for(uint64_t j=0;j<len;j++)
		{
			c_tmp.hashValue[j]=h1[(len+1)*i+j];
		}
		c_tmp.arrayID=h1[(i+1)*(len+1)-1];
		Insert_Value_bit(&p_root,c_tmp,para);
	}
	fclose(int_input);
	cout << "Constructing BplusTree is over!" <<endl;

	//mapping the first hash value to the first dimensional coordinate of HashAd
	for(uint64_t i=0;i<total_hash_number;i++)
	{
		uint64_t ID=MappingHashValueToID_bit(p_root,h1+i*(len+1),para);
		if((h1[(i+1)*(len+1)-1])!=ID)
		{
			cout << "error bplustree!" << endl;
			char y;
			cin >> y;
		}
	}

	struct NodeBit * p_leaf=Find_Minimal_Node_bit(p_root);

	FILE* outputFile;
	outputFile=fopen(outputFileName,"ab+");

	struct NodeBit *p_tmp=p_leaf;
	while(p_tmp!=NULL)
	{
		for(uint32_t i=0;i<p_tmp->Node_Size;i++)
		{
			fwrite(p_tmp->data[i].hashValue,sizeof(uint64_t),len,outputFile);
			fwrite(&(p_tmp->data[i].arrayID),sizeof(uint64_t),1,outputFile);
		}
		p_tmp=p_tmp->brother;
	}
	fclose(outputFile);
	destory_tree_bit(p_root,para);
}
void binaryCoutInt(uint64_t a)
{
	uint32_t b[64];
	uint64_t x=1;
	for(uint32_t i=0;i<64;i++)
	{
		b[i]=(x&a)>>i;
		x=x<<1;
	}
	for(uint32_t i=0;i<64;i++)
	{
		uint32_t j=63-i;
		cout << b[j];
	}
	cout << endl;
}
void binaryCoutChar(uint8_t a)
{
	uint32_t b[8];
	uint64_t x=1;
	for(uint32_t i=0;i<8;i++)
	{
		b[i]=(x&a)>>i;
		x=x<<1;
	}
	for(uint32_t i=0;i<8;i++)
	{
		uint32_t j=7-i;
		cout << b[j];
	}
	cout << endl;
}
void Coutsubstring(char* p,uint32_t l)
{
	for(uint32_t i=0;i<l;i++)
	{
		cout << p[i];
	}
	cout << endl;
}

