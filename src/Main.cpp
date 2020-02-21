#include "Basic.h"
#include "litool.h"
#include "Combine_Kmers.h"
#include "Extend_Kmer.h"
#include "Being_Short.h"
#include "Builder.h"

//int main(int argc, char** argv)
//{
////	char * p1;
////	char * p2;
////	uint32_t len;
////
////	p1=argv[1];
////	p2=argv[2];
////	len=atoi(argv[3]);
//	SortKmers(argv[1],argv[2],argv[3],argv[4],atoi(argv[5]));
//
////	SortOutKmersUsingBplusTreebinarytobinary(p1,p2,len);
//
//}

int main(int argc, char** argv)
{
	struct CDBG_para CDBG_para_tmp;
	CDBG_para_tmp.Task_Size_Single_Thread=24;
	CDBG_para_tmp.ad_label=0;
	CDBG_para_tmp.extend_len=5;
	CDBG_para_tmp.memory=5;
	CDBG_para_tmp.thread_num=1;
	CDBG_para_tmp.isFinal_label=0;
	CDBG_para_tmp.inout_label=2;
	CDBG_para_tmp.kmer_len=0;
	CDBG_para_tmp.p_ref=NULL;
	CDBG_para_tmp.method='N';
	CDBG_para_tmp.pos_label=0;
	CDBG_para_tmp.rc_label=0;
	CDBG_para_tmp.sort_label=1;
	CDBG_para_tmp.pl_label=0;
	if(argc<4)
	{
		cout << "the input parameter is wrong!" <<endl;
	}
	else
	{
		//-M : method
		//-L : kmer len
		//-l : extend len
		//-r : reference path
		//-a : ad_label
		//-i : inout_label
		//-t : thread num
		//-m : memory size
		//-s : single task size
		//-p : bit vector label
		//-S : sort result label
		for(uint32_t i=1;i<argc;i=i+2)
		{
			if(argv[i][0]=='-'&&argv[i][1]=='M')//Method
			{
				CDBG_para_tmp.method=argv[i+1][0];
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='L')//Len_kmer
			{
				CDBG_para_tmp.kmer_len=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='r')//reference
			{
				CDBG_para_tmp.p_ref=argv[i+1];
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='l')//len_extend
			{
				CDBG_para_tmp.extend_len=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='a')//ad_label
			{
				CDBG_para_tmp.ad_label=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='i')//inout_label
			{
				CDBG_para_tmp.inout_label=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='t')//thread
			{
				CDBG_para_tmp.thread_num=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='m')//memory
			{
				CDBG_para_tmp.memory=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='s')//single_task_size
			{
				CDBG_para_tmp.Task_Size_Single_Thread=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='p')//pos_label
			{
				//0:unuse pos list and unsave new pos list
				//1:unuse pos list and save new pos list
				//2:use pos list and unsave new pos list
				//3:use pos list and save new pos list

				CDBG_para_tmp.pos_label=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='S')//sorting results
			{
				CDBG_para_tmp.sort_label=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='c')//sorting results
			{
				CDBG_para_tmp.rc_label=atoi(argv[i+1]);
			}
			else if(argv[i][0]=='-'&&argv[i][1]=='P')//sorting results
			{
				CDBG_para_tmp.pl_label=atoi(argv[i+1]);
			}
		}
	}
	if(CDBG_para_tmp.method=='N')
	{
		cout << "error: method is without setting!"<< endl;
		return -1;
	}
	else if(CDBG_para_tmp.method!='B'&& CDBG_para_tmp.method!='U'&&\
			CDBG_para_tmp.method!='E'&& CDBG_para_tmp.method!='F'&&\
			CDBG_para_tmp.method!='P')
	{
		cout << "error from main: the method is setting wrong!" << endl;
	}
	else if(CDBG_para_tmp.method=='B'||CDBG_para_tmp.method=='C'||\
			CDBG_para_tmp.method=='E'||CDBG_para_tmp.method=='F'||\
			CDBG_para_tmp.method=='P')
	{
		if(CDBG_para_tmp.p_ref==NULL)
		{
			cout << "error: p_ref is without setting!"<< endl;
			return -1;
		}
		else
		{
			cout << "start read the paths of reference sequences!"<< endl;
			ifstream int_input;
			int_input.open(CDBG_para_tmp.p_ref,ios::in);
			uint64_t ref_number=0;
			char hv_tmp[64];
			while(int_input.getline(hv_tmp,63))
			{
				ref_number++;
			}
			int_input.close();
			cout << "the paths are in the file : "<<CDBG_para_tmp.p_ref << endl;
		}
	}
	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);
	switch(CDBG_para_tmp.method)
	{
	    //U : union the in and out branched kmer
		//B : being short
		//C : count kmer
		//E : extend
		//F : final extend
		//P : unipath
		/******************default setting********************************
			CDBG_para_tmp.Task_Size_Single_Thread=24;
			CDBG_para_tmp.ad_label=2;
			CDBG_para_tmp.extend_len=0;
			CDBG_para_tmp.memory=0;
			CDBG_para_tmp.thread_num=1;
			CDBG_para_tmp.inout_label=2;
			CDBG_para_tmp.kmer_len=0;
			CDBG_para_tmp.p_ref=NULL;
			CDBG_para_tmp.method='N';
			CDBG_para_tmp.kmer_path=NULL;
	 	 ****************************************************************/
		case 'U':
		{
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ U method calling: kmer length is without setting!"<< endl;
			}
			else if(CDBG_para_tmp.ad_label==2)
			{
				cout << "error from main.cpp @ U method calling: ad label is without setting!"<< endl;
			}
			else
			{
				struct timeval tvs,tve;
				gettimeofday(&tvs,NULL);
				cout <<"start..."<<endl;
				uion_Kmers_bit(CDBG_para_tmp.kmer_len,CDBG_para_tmp.ad_label);
				cout << "end..."<< endl;
				gettimeofday(&tve,NULL);
				double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
				cout <<"union kmers time is: "<<span<<endl;
			}
			break;
		}
		case 'B':
		{
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ B method calling: kmer_len is without setting!"<< endl;
				break;
			}
			else if(CDBG_para_tmp.p_ref==NULL)
			{
				cout << "error from main.cpp @ B method calling: p_ref is without setting!"<< endl;
				break;
			}
			else if(CDBG_para_tmp.ad_label==2)
			{
				cout << "error from main.cpp @ B method calling: ad label is without setting!"<< endl;
			}
			uint32_t thread_num=CDBG_para_tmp.thread_num;
			uint64_t kmer_len_being=CDBG_para_tmp.kmer_len;
			char* inputFileName;
			inputFileName=(char*)malloc(sizeof(char)*11);
			struct para_getN tmp_name;
			tmp_name.kmerlen=kmer_len_being;
			tmp_name.FileName=inputFileName;
			tmp_name.isMiddle=1;
			getFileName(tmp_name);

			struct para_being tmp_being;
			tmp_being.seq=CDBG_para_tmp.p_ref;
			tmp_being.kmer_len_being=kmer_len_being;
			tmp_being.outputFileName=inputFileName;
			tmp_being.thread_num=thread_num;
			tmp_being.Task_Size=CDBG_para_tmp.Task_Size_Single_Thread;
			tmp_being.ad_label=CDBG_para_tmp.ad_label;
			tmp_being.rc_label=CDBG_para_tmp.rc_label;
			struct timeval tvs,tve;
			gettimeofday(&tvs,NULL);
			cout <<"start..."<<endl;
//			Find_Short_Branched_Kmer_256bit(tmp_being);
			Find_Short_Branched_Kmer_256bit_novel(tmp_being);
			cout << "end..."<< endl;
			free(inputFileName);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"being short branched kmers finding time is: "<<span<<endl;
			break;
		}
		case 'C':
		{
			if(CDBG_para_tmp.memory==0)
			{
				cout << "error from main.cpp @ C method calling: memory size is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.inout_label==2)
			{
				cout << "error from main.cpp @ C method calling: in_out label is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ C method calling: kmer_len is without setting!"<< endl;
				return -1;
			}
			struct para_calTask tmp;
			tmp.InOut_label=CDBG_para_tmp.inout_label;
			tmp.pos_label=CDBG_para_tmp.pos_label;
			tmp.seq=CDBG_para_tmp.p_ref;
			tmp.original_kmer_length=CDBG_para_tmp.kmer_len;
			tmp.isFinal=1;
			tmp.memory_size=CDBG_para_tmp.memory;
			tmp.thread_num=CDBG_para_tmp.thread_num;
			tmp.Task_Size_Single_Thread=CDBG_para_tmp.Task_Size_Single_Thread;

			struct timeval tvs,tve;
			gettimeofday(&tvs,NULL);
			cout <<"start..."<<endl;
			Cal_Task(tmp);
			cout << "end..."<< endl;
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"counting kmers time is: "<<span<<endl;
			break;
		}
		case 'E':
		{
			if(CDBG_para_tmp.ad_label==2)
			{
				cout << "error from main.cpp @ E method calling: ad label is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.extend_len==0)
			{
				cout << "error from main.cpp @ E method calling: extend_len is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.memory==0)
			{
				cout << "error from main.cpp @ E method calling: memory size is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.inout_label==2)
			{
				cout << "error from main.cpp @ E method calling: inout_label is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ E method calling: kmer_len is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.p_ref==NULL)
			{
				cout << "error from main.cpp @ E method calling: p_ref is without setting!"<< endl;
				return -1;
			}

			struct para_calTask tmp;
			tmp.pos_label=CDBG_para_tmp.pos_label;
			tmp.seq=CDBG_para_tmp.p_ref;
			tmp.extend_length=CDBG_para_tmp.extend_len;
			tmp.original_kmer_length=CDBG_para_tmp.kmer_len;
			tmp.InOut_label=CDBG_para_tmp.inout_label;
			tmp.isFinal=0;
			tmp.ad_label=CDBG_para_tmp.ad_label;
			tmp.memory_size=CDBG_para_tmp.memory;
			tmp.thread_num=CDBG_para_tmp.thread_num;
			tmp.Task_Size_Single_Thread=CDBG_para_tmp.Task_Size_Single_Thread;
			tmp.sort_label=CDBG_para_tmp.sort_label;
			tmp.rc_label=CDBG_para_tmp.rc_label;
			struct timeval tvs,tve;
			gettimeofday(&tvs,NULL);
			cout <<"start..."<<endl;
			Cal_Task(tmp);
			cout << "end..."<< endl;
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"Normal extending kmers time is: "<<span<<endl;
			break;
		}
		case 'F':
		{
			if(CDBG_para_tmp.ad_label==2)
			{
				cout << "error from main.cpp @ F method calling: ad label is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.extend_len==0)
			{
				cout << "error from main.cpp @ F method calling: extend_len is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.memory==0)
			{
				cout << "error from main.cpp @ F method calling: memory size is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.inout_label==2)
			{
				cout << "error from main.cpp @ F method calling: inout_label is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ F method calling: kmer_len is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.p_ref==NULL)
			{
				cout << "error from main.cpp @ F method calling: p_ref is without setting!"<< endl;
				return -1;
			}
			struct para_calTask tmp;
			tmp.pos_label=CDBG_para_tmp.pos_label;
			tmp.seq=CDBG_para_tmp.p_ref;
			tmp.extend_length=CDBG_para_tmp.extend_len;
			tmp.original_kmer_length=CDBG_para_tmp.kmer_len;
			tmp.InOut_label=CDBG_para_tmp.inout_label;
			tmp.isFinal=2;
			tmp.ad_label=CDBG_para_tmp.ad_label;
			tmp.memory_size=CDBG_para_tmp.memory;
			tmp.thread_num=CDBG_para_tmp.thread_num;
			tmp.Task_Size_Single_Thread=CDBG_para_tmp.Task_Size_Single_Thread;
			tmp.sort_label=CDBG_para_tmp.sort_label;
			tmp.rc_label=CDBG_para_tmp.rc_label;
			struct timeval tvs,tve;
			gettimeofday(&tvs,NULL);
			cout <<"start..."<<endl;
			tmp.isFinal=2;
			Cal_Task(tmp);
			cout << "end..."<< endl;
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"Final extending kmers time is: "<<span<<endl;
			break;
		}
		case 'P':
		{
			if(CDBG_para_tmp.kmer_len==0)
			{
				cout << "error from main.cpp @ P method calling: kmer_len is without setting!"<< endl;
				return -1;
			}
			if(CDBG_para_tmp.p_ref==NULL)
			{
				cout << "error from main.cpp @ P method calling: p_ref is without setting!"<< endl;
				return -1;
			}
			struct timeval tvs,tve;
			gettimeofday(&tvs,NULL);
			cout <<"start..."<<endl;
			struct para_findUnipath tmps;
			tmps.seq=CDBG_para_tmp.p_ref;
			tmps.pos_label=CDBG_para_tmp.pos_label;
			tmps.original_kmer_length=CDBG_para_tmp.kmer_len;
			tmps.thread_num=CDBG_para_tmp.thread_num;
			tmps.rc_label=CDBG_para_tmp.rc_label;
			tmps.memory_size=CDBG_para_tmp.memory;
			tmps.pl_label=CDBG_para_tmp.pl_label;
			FindUniPath_Novel(tmps);
			cout << "end..."<< endl;
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"time cost for finding unipath is: "<<span<<endl;
			break;
		}
	}

	cout << "finished!"<< endl;
	gettimeofday(&tve_total,NULL);
	double span_check_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total cost time is: "<<span_check_total<<endl;

}


