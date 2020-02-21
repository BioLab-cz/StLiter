/*
 * BplusTreeBit.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: bio
 */
#include "BplusTreeBit.h"

void Node_initial_bit(struct NodeBit * p)
{
	p->Node_Size=0;
	p->leaf_label=1;
	p->numAsChild=-1;

	p->parent=NULL;
	p->brother=NULL;
	for(uint32_t i=0;i<M;i++)
	{
		p->p_child[i]=NULL;
		p->data[i].hashValue=NULL;
		p->data[i].ad=0;
		p->data[i].p_posList=NULL;
		p->data[i].p_posList_length=0;
	}
}
void Insert_Value_bit(struct NodeBit** p_root,struct nodeBit insert_data,struct bit256KmerPara para)
{

	struct NodeBit*p_root_tmp= * p_root;
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
			cout << "error from bplustree insert!" << endl;
			return;
		}
		else if(pos_label==0)
		{
			pos=p_root_tmp->Node_Size;
		}
	}

//	cout << "pos is:" << pos << endl;

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
			if(p_root_tmp->leaf_label==1)//ycy20191028 for unipath counter
			{
				p_root_tmp->data[l+1].p_posList_length=p_root_tmp->data[l].p_posList_length;
			}
		}
		p_root_tmp->data[pos].arrayID=insert_data.arrayID;
		p_root_tmp->data[pos].p_posList_length=0;//ycy20191028 for unipath counter
		if(p_root_tmp->data[pos].hashValue==NULL)
		{
			p_root_tmp->data[pos].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
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

		p_root_tmp->Node_Size++;

		if(pos==0&&p_root_tmp->parent!=NULL)
		{
			NodeBit* p_parent_update=p_root_tmp->parent;
			p_parent_update->data[p_root_tmp->numAsChild].arrayID=\
					p_root_tmp->data[0].arrayID;
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_parent_update->data[p_root_tmp->numAsChild].hashValue[j]=\
					p_root_tmp->data[0].hashValue[j];
			}
			NodeBit* p1=p_root_tmp;
			NodeBit* p2=p_parent_update;
			while(p1->numAsChild==0&&p2->parent!=NULL)
			{
				p2->parent->data[p2->numAsChild].arrayID=\
						p2->data[0].arrayID;
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p2->parent->data[p2->numAsChild].hashValue[j]=\
							p2->data[0].hashValue[j];
				}
				p1=p2;
				p2=p2->parent;
			}
		}

	}

}
void Insert_Node_bit(struct NodeBit** p_root,struct NodeBit* p_inserted_node,\
		struct nodeBit insert_data,struct NodeBit* p_middle,struct bit256KmerPara para)
{
	struct NodeBit*p_root_tmp= p_inserted_node;

	uint32_t pos;
	uint32_t pos_label=0;
	for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
	{
		if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
		{
			pos=i;
			pos_label=1;
			break;
		}
	}
	if(pos_label==0)
	{
		pos=p_root_tmp->Node_Size;
	}
	if(p_root_tmp->Node_Size==M)
	{
		Divide_Node_bit(p_root,p_root_tmp,insert_data,pos,p_middle,0,0,para);
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

			p_root_tmp->p_child[l+1]=p_root_tmp->p_child[l];
			p_root_tmp->p_child[l+1]->numAsChild++;
		}
		p_root_tmp->data[pos].arrayID=insert_data.arrayID;
		if(p_root_tmp->data[pos].hashValue==NULL)
		{
			p_root_tmp->data[pos].hashValue=\
					(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
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
		p_root_tmp->p_child[pos]=p_middle;
		p_middle->numAsChild=pos;
		p_middle->parent=p_root_tmp;

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
							(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
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
void Divide_Node_bit(struct NodeBit** p_root,struct NodeBit* p_divide_node,struct nodeBit insert_data,\
		uint32_t pos,struct NodeBit* insert_child,uint32_t ad_label,uint32_t poslist_label,struct bit256KmerPara para)
{
	struct NodeBit * parent_tmp;
	parent_tmp=p_divide_node->parent;

	struct NodeBit * p_add=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_add);
	p_add->leaf_label=p_divide_node->leaf_label;
	p_add->Node_Size=M/2+1;
	if(p_divide_node->leaf_label==1)
	{
		p_add->brother=p_divide_node->brother;
		p_divide_node->brother=p_add;
	}

	if(pos<=p_divide_node->Node_Size-M/2-1)
	{
		uint32_t shift=(p_divide_node->Node_Size-M/2-1);
		for(uint32_t i=shift;i<p_divide_node->Node_Size;i++)
		{
			if(p_add->data[i-shift].hashValue==NULL)
			{
				p_add->data[i-shift].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift].hashValue[j]=\
							p_divide_node->data[i].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift].hashValue[j]=\
							p_divide_node->data[i].hashValue[j];
				}
			}
			p_add->data[i-shift].arrayID=\
					p_divide_node->data[i].arrayID;
//////////////////////////////////////////////////////////////////_
			if(ad_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift].ad=\
					p_divide_node->data[i].ad;
			}
			if(poslist_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift].p_posList_length=\
					p_divide_node->data[i].p_posList_length;
//				p_add->data[i-shift].p_posList=\
//					p_divide_node->data[i].p_posList;
//
//				p_divide_node->data[i].p_posList=NULL;
				p_divide_node->data[i].p_posList_length=0;
			}
//////////////////////////////////////////////////////////////////
			p_add->p_child[i-shift]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift;
			}
			p_divide_node->p_child[i]=NULL;
		}


		for(uint32_t i=pos;i<=shift-1;i++)
		{
			uint32_t l=shift-1-(i-pos);
			if(p_divide_node->data[l+1].hashValue==NULL)
			{
				p_divide_node->data[l+1].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_divide_node->data[l+1].hashValue[j]=p_divide_node->data[l].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_divide_node->data[l+1].hashValue[j]=p_divide_node->data[l].hashValue[j];
				}
			}
			p_divide_node->data[l+1].arrayID=p_divide_node->data[l].arrayID;
			p_divide_node->p_child[l+1]=p_divide_node->p_child[l];
//////////////////////////////////////////////////////////////////
			if(ad_label==1&&p_add->leaf_label==1)
			{
				p_divide_node->data[l+1].ad=p_divide_node->data[l].ad;
			}
			if(poslist_label==1&&p_add->leaf_label==1)
			{
				p_divide_node->data[l+1].p_posList_length=\
						p_divide_node->data[l].p_posList_length;
//				p_divide_node->data[l+1].p_posList=\
//						p_divide_node->data[l].p_posList;
			}
//////////////////////////////////////////////////////////////////
			if(p_divide_node->p_child[l+1]!=NULL)
			{
				p_divide_node->p_child[l+1]->numAsChild=l+1;
			}
		}

//////////////////////////////////////////////////////////////////
		if(poslist_label==1&&p_add->leaf_label==1)
		{
//			p_divide_node->data[pos].p_posList=NULL;
			p_divide_node->data[pos].p_posList_length=0;
		}
//////////////////////////////////////////////////////////////////

		p_divide_node->data[pos].arrayID=insert_data.arrayID;
		if(p_divide_node->data[pos].hashValue==NULL)
		{
			p_divide_node->data[pos].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_divide_node->data[pos].hashValue[j]=insert_data.hashValue[j];
			}
		}
		else
		{
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_divide_node->data[pos].hashValue[j]=insert_data.hashValue[j];
			}
		}
//////////////////////////////////////////////////////////////////
		if(ad_label==1&&p_add->leaf_label==1)
		{
			p_divide_node->data[pos].ad=insert_data.ad;
		}

		if(poslist_label==1&&p_add->leaf_label==1)
		{
			p_divide_node->data[pos].ad=insert_data.ad;//ycy20191028 for unipath counter
//			if(p_divide_node->data[pos].p_posList_length==0)
//			{
//				uint32_t tmp_len=p_divide_node->data[pos].p_posList_length;
//				p_divide_node->data[pos].p_posList=(uint64_t*)malloc(sizeof(uint64_t)*(tmp_len+1));
//				p_divide_node->data[pos].p_posList_length++;
//
//				p_divide_node->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
//			}
//			else
//			{
//				uint32_t tmp_len=p_divide_node->data[pos].p_posList_length;
//				p_divide_node->data[pos].p_posList=(uint64_t*)realloc(p_divide_node->data[pos].p_posList,sizeof(uint64_t)*(tmp_len+1));
//				p_divide_node->data[pos].p_posList_length++;
//
//				p_divide_node->data[pos].p_posList[tmp_len]=insert_data.p_posList[0];
//			}
		}
//////////////////////////////////////////////////////////////////
		p_divide_node->p_child[pos]=insert_child;

		if(insert_child!=NULL)
		{
			insert_child->parent=p_divide_node;
			insert_child->numAsChild=pos;
		}
		p_divide_node->Node_Size=p_divide_node->Node_Size-M/2;

		if(pos==0&&p_divide_node->parent!=NULL)
		{
			NodeBit* p_parent_update=p_divide_node->parent;
			p_parent_update->data[p_divide_node->numAsChild].arrayID=\
					p_divide_node->data[0].arrayID;
			if(p_parent_update->data[p_divide_node->numAsChild].hashValue==NULL)
			{
				p_parent_update->data[p_divide_node->numAsChild].hashValue=\
						(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_parent_update->data[p_divide_node->numAsChild].hashValue[j]=\
							p_divide_node->data[0].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_parent_update->data[p_divide_node->numAsChild].hashValue[j]=\
							p_divide_node->data[0].hashValue[j];
				}
			}
			NodeBit* p1=p_divide_node;
			NodeBit* p2=p_divide_node->parent;
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

		if(parent_tmp==NULL)
		{
			update_root_bit(p_root,p_add,para);
		}
		else
		{
			Insert_Node_bit(p_root,parent_tmp,p_add->data[0],p_add,para);
		}
	}
	else
	{
		uint32_t shift=(p_divide_node->Node_Size-M/2);
		for(uint32_t i=p_divide_node->Node_Size-M/2;i<=pos-1;i++)
		{
			if(p_add->data[i-shift].hashValue==NULL)
			{
				p_add->data[i-shift].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift].hashValue[j]=p_divide_node->data[i].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift].hashValue[j]=p_divide_node->data[i].hashValue[j];
				}
			}
			p_add->data[i-shift].arrayID=p_divide_node->data[i].arrayID;
//////////////////////////////////////////////////////////////////
			if(ad_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift].ad=p_divide_node->data[i].ad;
			}
			if(poslist_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift].p_posList_length=p_divide_node->data[i].p_posList_length;
//				p_add->data[i-shift].p_posList=p_divide_node->data[i].p_posList;

				p_divide_node->data[i].p_posList_length=0;
//				p_divide_node->data[i].p_posList=NULL;
			}
//////////////////////////////////////////////////////////////////
			p_add->p_child[i-shift]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift;
			}
			p_divide_node->p_child[i]=NULL;
		}

		p_add->data[pos-shift].arrayID=insert_data.arrayID;
		if(p_add->data[pos-shift].hashValue==NULL)
		{
			p_add->data[pos-shift].hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_add->data[pos-shift].hashValue[j]=insert_data.hashValue[j];
			}
		}
		else
		{
			for(uint32_t j=0;j<para.kmer64Len;j++)
			{
				p_add->data[pos-shift].hashValue[j]=insert_data.hashValue[j];
			}
		}
//////////////////////////////////////////////////////////////////
		if(ad_label==1&&p_add->leaf_label==1)
		{
			p_add->data[pos-shift].ad=insert_data.ad;
		}
		if(poslist_label==1&&p_add->leaf_label==1)
		{
			p_add->data[pos-shift].p_posList_length=insert_data.p_posList_length;//ycy20191028 for unipath counter
//			if(p_add->data[pos-shift].p_posList_length==0)
//			{
//				uint32_t tmp_len=p_add->data[pos-shift].p_posList_length;
//				p_add->data[pos-shift].p_posList=(uint64_t*)malloc(sizeof(uint64_t)*(tmp_len+1));
//				p_add->data[pos-shift].p_posList_length++;
//
//				p_add->data[pos-shift].p_posList[tmp_len]=insert_data.p_posList[0];
//			}
//			else
//			{
//				uint32_t tmp_len=p_add->data[pos-shift].p_posList_length;
//				p_add->data[pos-shift].p_posList=(uint64_t*)realloc(p_add->data[pos-shift].p_posList,sizeof(uint64_t)*(tmp_len+1));
//				p_add->data[pos-shift].p_posList_length++;
//
//				p_add->data[pos-shift].p_posList[tmp_len]=insert_data.p_posList[0];
//			}
		}
//////////////////////////////////////////////////////////////////
		p_add->p_child[pos-shift]=insert_child;
		if(insert_child!=NULL)
		{
			insert_child->parent=p_add;
			insert_child->numAsChild=pos-shift;
		}
		for(uint32_t i=pos;i<=p_divide_node->Node_Size-1;i++)
		{
			if(p_add->data[i-shift+1].hashValue==NULL)
			{
				p_add->data[i-shift+1].hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift+1].hashValue[j]=p_divide_node->data[i].hashValue[j];
				}
			}
			else
			{
				for(uint32_t j=0;j<para.kmer64Len;j++)
				{
					p_add->data[i-shift+1].hashValue[j]=p_divide_node->data[i].hashValue[j];
				}
			}
			p_add->data[i-shift+1].arrayID=p_divide_node->data[i].arrayID;
//////////////////////////////////////////////////////////////////
			if(ad_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift+1].ad=p_divide_node->data[i].ad;
			}
			if(poslist_label==1&&p_add->leaf_label==1)
			{
				p_add->data[i-shift+1].p_posList_length=p_divide_node->data[i].p_posList_length;
//				p_add->data[i-shift+1].p_posList=p_divide_node->data[i].p_posList;

				p_divide_node->data[i].p_posList_length=0;
//				p_divide_node->data[i].p_posList=NULL;
			}
//////////////////////////////////////////////////////////////////
			p_add->p_child[i-shift+1]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift+1;
			}
			p_divide_node->p_child[i]=NULL;
		}
		p_divide_node->Node_Size=p_divide_node->Node_Size-M/2;

		if(parent_tmp==NULL)
		{
			update_root_bit(p_root,p_add,para);
		}
		else
		{
			Insert_Node_bit(p_root,parent_tmp,p_add->data[0],p_add,para);
		}
	}
}
void update_root_bit(struct NodeBit** p_root,struct NodeBit* p2,struct bit256KmerPara para)
{
	struct NodeBit *root_tmp=*p_root;
	struct NodeBit * p_tmp=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_tmp);

	p_tmp->Node_Size=2;
	p_tmp->leaf_label=0;

	if(p_tmp->data[0].hashValue==NULL)
	{
		p_tmp->data[0].hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			p_tmp->data[0].hashValue[j]=root_tmp->data[0].hashValue[j];
		}
	}
	else
	{
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			p_tmp->data[0].hashValue[j]=root_tmp->data[0].hashValue[j];
		}
	}
	p_tmp->data[0].arrayID=root_tmp->data[0].arrayID;
	p_tmp->p_child[0]=root_tmp;
	root_tmp->numAsChild=0;
	root_tmp->parent=p_tmp;

	if(p_tmp->data[1].hashValue==NULL)
	{
		p_tmp->data[1].hashValue=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			p_tmp->data[1].hashValue[j]=p2->data[0].hashValue[j];
		}
	}
	else
	{
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			p_tmp->data[1].hashValue[j]=p2->data[0].hashValue[j];
		}
	}
	p_tmp->data[1].arrayID=p2->data[0].arrayID;
	p_tmp->p_child[1]=p2;
	p2->numAsChild=1;
	p2->parent=p_tmp;

	*p_root=p_tmp;
}
void createBplusTree_bit(struct NodeBit** p_root,struct nodeBit* p_data,uint32_t p_data_size,struct bit256KmerPara para)
{
	for(uint32_t i=0;i<p_data_size;i++)
	{
		Insert_Value_bit(p_root,p_data[i],para);
	}
}
void destory_tree_bit(struct NodeBit* p_root,struct bit256KmerPara para)
{
	if(p_root!=NULL)
	{
		for(uint32_t i=0;i<p_root->Node_Size;i++)
		{
			destory_tree_bit(p_root->p_child[i],para);
			free(p_root->data[i].hashValue);
		}
		free(p_root);
	}
}
int64_t addAD1_bit(struct NodeBit** p_root,uint64_t* hashValue,uint8_t c, struct bit256KmerPara para)
{
	struct NodeBit* p_root_tmp = *p_root;
	while(1)
	{
		if(p_root_tmp->leaf_label==1)
		{
			break;
		}
		uint32_t child_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
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

	if(p_root_tmp->Node_Size==0)
	{
		return -2;
	}
	else
	{
		int32_t pos_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==2)
			{
				p_root_tmp->data[i].ad=p_root_tmp->data[i].ad|c;
				return 1;
			}
			else if(cmp256BitKmer(hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				return -1;
			}
		}
	}
	return -3;
}
int64_t addAD2_bit(struct NodeBit* p_root,uint64_t* hashValue,uint8_t c,struct bit256KmerPara para)
{
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(i!=0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				return -1;
//				p_tmp=NULL;
//				label=1;
//				break;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp==NULL)
	{
		return -2;
	}
	else
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==1)
			{
				return -3;
			}
			else if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==2)
			{
				p_tmp->data[i].ad=p_tmp->data[i].ad|c;
				return 1;
			}
		}
		return -4;
	}
}
int64_t addAD_bit(struct NodeBit* p_root,uint64_t* hashValue,uint8_t c, struct bit256KmerPara para)
{
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(i!=0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				return -1;
//				p_tmp=NULL;
//				label=1;
//				break;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp==NULL)
	{
		return -2;
	}
	else
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==1)
			{
				return -3;
			}
			else if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==2)
			{
				p_tmp->data[i].ad=p_tmp->data[i].ad|c;
				return 1;
			}
		}
		return -4;
	}
	return -5;
}
struct NodeBit* Find_Minimal_Node_bit(struct NodeBit* p_root)
{
	struct NodeBit* p_tmp=p_root;
	while(p_tmp->leaf_label!=1)
	{
		p_tmp=p_tmp->p_child[0];
	}
	return p_tmp;
}
int64_t MappingHashValueToID_bit(struct NodeBit* p_root,uint64_t* hashValue,struct bit256KmerPara para)
{
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			uint32_t lp=cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len);
			if(i!=0&&lp==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&lp==0)
			{
				return -1;
//				p_tmp=NULL;
//				label=1;
//				break;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp==NULL)
	{
		return -1;
	}
	else
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			uint32_t lp;
			lp=cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len);
			if(lp==1)
			{
				return -1;
			}
			else if(lp==2)
			{
				p_tmp->data[i].p_posList_length++;
				if(p_tmp->data[i].arrayID==-1)
				{
					cout << "error: from mapping inner!" << endl;
				}
				return p_tmp->data[i].arrayID;
			}
		}
		return -1;
	}
}
void count_hashValue_bit(struct NodeBit* p_root,uint64_t* hashValue,struct bit256KmerPara para)
{
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
 		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(i!=0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=NULL;
				label=1;
				break;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp!=NULL)
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==1)
			{
				break;
			}
			else if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==2)
			{
				p_tmp->data[i].arrayID++;
				break;
			}
		}
	}
}
int64_t MappingHashValueToID_bit_unipath_count(struct NodeBit* p_root,uint64_t* hashValue,struct bit256KmerPara para)
{
	//2019.10.23
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(i!=0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				return -1;
//				p_tmp=NULL;
//				label=1;
//				break;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp==NULL)
	{
		return -1;
	}
	else
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==1)
			{
				return -1;
			}
			else if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==2)
			{
				p_tmp->data[i].p_posList_length++;
				if(p_tmp->data[i].arrayID==-1)
				{
					cout << "error: from mapping inner!" << endl;
				}
				return p_tmp->data[i].arrayID;
			}
		}
		return -1;
	}
}
