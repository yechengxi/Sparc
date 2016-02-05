#ifndef __BASIC_DATA_STRUCTURE_H
#define __BASIC_DATA_STRUCTURE_H

#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"
using namespace std;
//typedef unsigned __int64 uint64_t;
//typedef __int64 int64_t;
//typedef unsigned __int8 uint8_t;
//typedef unsigned __int16 uint16_t;
//typedef unsigned __int32 uint32_t;
//typedef __int32 int32_t;



// These are the structures to save the k-mers.32 bases, 64 bases, 96 bases, 128 bases.
struct kmer_t
{
	uint64_t kmer;
};

struct kmer_t2
{
	uint64_t kmer[2];
};

struct kmer_t3
{
	uint64_t kmer[3];
};

struct kmer_t4
{
	uint64_t kmer[4];
};




//read structure

struct read_t
{
	char tag[1000];
	bool error_nt[1000];
	char c_seq[10000];//char representation
	//char *c_seq;//char representation
	
	//uint64_t read_bits[10000];//bit representation
	uint64_t *read_bits;
	//char read[1000];//char representation
	int readLen;// read length
	int read_idx;
};



struct ref_read_t
{
	char tag[1000];
	uint64_t *read_bits;//bit representation
	size_t read_idx;
	int alloc_sz;
	int readLen;// read length
	int contig_no;

};


struct align_profile
{
	vector<int> match_vec, mismatch_vec, deletion_vec, insertion_vec;
};


struct query_info
{
	string qName, tName, qAlignedSeq, matchPattern, tAlignedSeq;
	int qLength, qStart, qEnd, tLength, tStart, tEnd, score, numMatch, numMismatch, numIns, numDel, mapQV;
	char qStrand, tStrand;
	size_t read_idx;
	int report_b, report_e;
	int n_exist;
	int n_new;
	bool Patch,Fill;
	int Patch_K;
	int Patch_D; 
	int Patch_G;
};




struct consensus_edge_node
{
	int edge_cov;
	struct consensus_node *node_ptr;
	struct consensus_edge_node *nxt_edge;
};


struct sparse_consensus_edge_node
{
	uint32_t edge;
	int32_t edge_cov: 24, len : 8;
	struct sparse_consensus_node *node_ptr;
	struct sparse_consensus_edge_node *nxt_edge;
};





struct consensus_node
{
	uint32_t kmer;
	uint32_t coord;
	int64_t score;
	uint32_t cns_coord; 
	uint32_t cov : 28, used : 1, in_backbone : 1, in_cns_backbone : 1;
	consensus_edge_node *left;
	consensus_edge_node *right;
	consensus_node * last_node;
};

struct sparse_consensus_node
{
	uint32_t kmer;
	int64_t score; 
	uint32_t coord;
	uint32_t cns_coord;
	uint32_t cov : 28, used : 1, in_backbone : 1, in_cns_backbone : 1;
	sparse_consensus_edge_node *left;
	sparse_consensus_edge_node *right;
	sparse_consensus_node * last_node;
};



struct backbone_info
{
	vector<consensus_node *> node_vec;
	vector<sparse_consensus_node *> sparse_node_vec;
	vector<uint64_t> cov_vec;
	vector<uint64_t> cnt_vec;
	string backbone;
	int64_t n_nodes;
	int64_t n_edges;
	int CovTh;
	int ScoringMethod;
	size_t ref_matched, ref_mismatch;
	int boost;
	int gap;
	double threshold;
	string ContigPrefix;
};

struct read_index
{
	vector< map<uint64_t,bool> > repeat_maps;
	uint64_t repeat_cnt;
	int MaxMatch;
	vector<int> read_len_vt;
};


struct reads_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerBlock;
	int current_block;
	int current_byte;
	int current_read;
	map<int,uint64_t *> read_ptr;
	vector<int> read_len_vt;
};

//contig graph

struct contigs_info
{
	int total_contigs;
	int K_size;
	vector<int> contig_sz_vt,kmer_cnt_vt,comp_vt;
	vector<int> contigs_hp_b,contigs_hp_e;
	vector<string> contigs_str;
	map<int, vector<int> > Length_ID;
	map<int, vector<int> > Cov_Length;
	//map<int, vector<int> > ctg_in_scf;
	//vector<vector<int> > scaffolds;
	//vector<vector<int> > gaps_in_scaffolds;

	vector < vector<int>::iterator> LengthRank;
	vector<int> cov_vt;
	vector<struct c_info> c_info_vt;
	vector< map<int,struct scaffold_contig_info> > scaffold_adjacency_left,scaffold_adjacency_right;
	vector< map<int,struct adjacent_contig_info> > contig_adjacency_left,contig_adjacency_right;	
};




//path information in the BFS bubble removal

bool get_a_fasta_read(ifstream & fasta_in, string &tag, string &str, string & n_tag)
{
	
	ifstream tmp_ifstream;
	string temp;
	if(!getline(fasta_in,temp))
	{return 0;}
	if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
	{temp.resize(temp.size()-1);}

	str.clear();
	if(temp[0]=='>')
	{tag=temp;}
	else
	{
		tag=n_tag;
		str=temp;
	}

	
	while(getline(fasta_in,temp))
	{
		
		if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
		{temp.resize(temp.size()-1);}

		if((temp.size()>0&&(temp[0]=='>'||temp[0]=='\n'||temp[0]=='\r')))
		{
			n_tag=temp;
			return 1;
		}
		else
		{
			str+=temp;
			
		}
		
	}
	return 1;
}


bool get_a_fastq_read(ifstream & fastq_in, string &tag, string &seq, string & quality)
{
	
	ifstream tmp_ifstream;
	string temp;
	if(!getline(fastq_in,temp))
	{return 0;}
	seq.clear();
	if(temp[0]=='@')
	{tag=temp;}
	else
	{
		return 0;
	}
	getline(fastq_in,seq);//seq
	if(seq[seq.size()-1]=='\n'||seq[seq.size()-1]=='\r')
	{seq.resize(seq.size()-1);}
	getline(fastq_in,temp);//'+'
	getline(fastq_in,quality);
	if(quality[quality.size()-1]=='\n'||quality[quality.size()-1]=='\r')
	{quality.resize(quality.size()-1);}

	return 1;
}

//left shift and right shift of shift_sz bits of the whole bit array, arr_sz is the array length
static inline void L_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{

	uint64_t temp_arr[100];

	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/

	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=0;i<arr_sz;++i)
	{
		if(i+jmp+1<arr_sz)
		{

			uint64_t tt=0;
			if(offset==0)
			{
				tt=0;
			}
			else
			{
				tt=(bitsarr[i+jmp+1]>>(64-offset));
			}
			temp_arr[i]=((bitsarr[i+jmp]<<offset)|tt);
		}
		if(i+jmp+1==arr_sz)
		{temp_arr[i]=bitsarr[i+jmp]<<offset;}
		if(i+jmp+1>arr_sz)
		{temp_arr[i]=0;}

	}

	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/
}

static inline void R_shift_NB(uint64_t * bitsarr, int shift_sz,int arr_sz)
{
	uint64_t temp_arr[100];
	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/
	int jmp=shift_sz/64;
	int offset=shift_sz%64;

	for (int i=arr_sz-1;i>=0;--i)
	{
		if(i-jmp>0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset)|(bitsarr[i-jmp-1]<<(64-offset));}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp==0)
		{
			if(offset>0)
			{temp_arr[i]=(bitsarr[i-jmp]>>offset);}
			else
			{temp_arr[i]=bitsarr[i-jmp];}
		}
		if (i-jmp<0)
		{temp_arr[i]=0;}

	}
	memcpy(bitsarr,temp_arr,sizeof(uint64_t)*arr_sz);
	/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/

}

// get reverse complement of a k-mer.
static inline uint64_t get_rev_comp_seq(uint64_t seq, int seq_size)
{
	seq =~seq;

	seq = ((seq & 0x3333333333333333 )<< 2) | ((seq & 0xCCCCCCCCCCCCCCCC )>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0 )>> 4);
	seq = ((seq & 0x00FF00FF00FF00FF )<< 8) | ((seq & 0xFF00FF00FF00FF00 )>> 8);
	seq = ((seq & 0x0000FFFF0000FFFF )<<16) | ((seq & 0xFFFF0000FFFF0000 )>>16);
	seq = ((seq & 0x00000000FFFFFFFF )<<32) | ((seq & 0xFFFFFFFF00000000 )>>32);

	return seq >> (64 - (seq_size*2));
}

static inline uint64_t* get_rev_comp_seq_arr(uint64_t *seq_arr, int seq_size,int arr_sz)
{


	int tot_bits=arr_sz*64;
	for(int i=0;i<arr_sz;++i)
	{
		seq_arr[i]=~seq_arr[i];
		seq_arr[i] = ((seq_arr[i] & 0x3333333333333333 )<< 2) | ((seq_arr[i] & 0xCCCCCCCCCCCCCCCC )>> 2);
		seq_arr[i] = ((seq_arr[i] & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq_arr[i] & 0xF0F0F0F0F0F0F0F0 )>> 4);
		seq_arr[i] = ((seq_arr[i] & 0x00FF00FF00FF00FF )<< 8) | ((seq_arr[i] & 0xFF00FF00FF00FF00 )>> 8);
		seq_arr[i] = ((seq_arr[i] & 0x0000FFFF0000FFFF )<<16) | ((seq_arr[i] & 0xFFFF0000FFFF0000 )>>16);
		seq_arr[i] = ((seq_arr[i] & 0x00000000FFFFFFFF )<<32) | ((seq_arr[i] & 0xFFFFFFFF00000000 )>>32);
	}

	int j=0,k=arr_sz-1;
	for (;j<k;++j,--k)
	{
		uint64_t temp;
		temp=seq_arr[j];
		seq_arr[j]=seq_arr[k];
		seq_arr[k]=temp;
	}

	R_shift_NB(seq_arr,tot_bits-(seq_size*2),arr_sz);
	return seq_arr;
	//return seq >> (64 - (seq_size*2));
}

// get sub bit array from a bit array.
 inline void get_sub_arr(uint64_t * bitsarr_in,int bitsarr_len,int begin_pos,int sub_sz,uint64_t * bitsarr_out)
{
	if(bitsarr_len<sub_sz)
	{cout<<"Error! Input kmer too short."<<bitsarr_len <<" "<<sub_sz<<endl;return;}
	int arr_sz_in=bitsarr_len/32+1;
	int rem=bitsarr_len%32;
	if(rem==0)
	{arr_sz_in--;}

	int arr_sz_out=sub_sz/32+1;
	if(sub_sz%32==0)
	{arr_sz_out--;}

	uint64_t temp_arr[10];
	memset(temp_arr,0,sizeof(temp_arr));

	memset(bitsarr_out,0,sizeof(uint64_t)*arr_sz_out);

	int rem2=(32-rem+begin_pos)%32;
	int block_beg=(32-rem+begin_pos)/32;
	if(rem==0)
	{block_beg--;}

	int rem3=(32-rem+begin_pos+sub_sz)%32;
	int block_end=(32-rem+begin_pos+sub_sz)/32;
	if(rem3!=0)
	{rem3=32-rem3;}
	else
	{
		block_end--;
	}
	if(rem==0)
	{block_end--;}

	int orig_sz=(block_end-block_beg+1);
	memcpy(temp_arr,&bitsarr_in[block_beg],orig_sz*sizeof(uint64_t));
	L_shift_NB(temp_arr,rem2*2,orig_sz);
	R_shift_NB(temp_arr,(rem2+rem3)%32*2,arr_sz_out);
	memcpy(bitsarr_out,temp_arr,arr_sz_out*sizeof(uint64_t));


}

 uint64_t* str2bitsarr(const char * c_str,int len, uint64_t* b_str,int arr_sz )
{

	for (int  k=0;k<arr_sz;++k)
	{
		b_str[k]=0;
	}
	int arr_sz_needed=len/32+1;
	int rem=len%32;
	if(rem==0)
	{arr_sz_needed--;}

	int beg_arr_idx=arr_sz-arr_sz_needed;
	if(rem==0&&arr_sz_needed>0)
	{rem=32;}
	for (int k=0;k<len;k++)
	{
		if(rem==0)
		{beg_arr_idx++;rem=32;}


		switch(c_str[k])
		{
		case ('A'):case ('a'):case ('0'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			rem--;
			//b_str<<=2;
			break;


		case ('C'):case ('c'):case ('1'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			rem--;
//			++(b_str<<=2);
			break;


		case 'G':case 'g':case '2':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			rem--;//(++(b_str<<=1))<<=1;
			break;

		case 'T':case 't':case '3':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			rem--;
			//++((++(b_str<<=1))<<=1);
			break;
		default:
			return b_str;
		}

	//	cout<<b_str<<endl;
	}
	return b_str;
}

 //hash functions
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed )
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	unsigned int h1 = seed ^ len;
	unsigned int h2 = 0;

	const unsigned int * data = (const unsigned int *)key;

	while(len >= 8)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;

		unsigned int k2 = *data++;
		k2 *= m; k2 ^= k2 >> r; k2 *= m;
		h2 *= m; h2 ^= k2;
		len -= 4;
	}

	if(len >= 4)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;
	}

	switch(len)
	{
	case 3: h2 ^= ((unsigned char*)data)[2] << 16;
	case 2: h2 ^= ((unsigned char*)data)[1] << 8;
	case 1: h2 ^= ((unsigned char*)data)[0];
			h2 *= m;
	};

	h1 ^= h2 >> 18; h1 *= m;
	h2 ^= h1 >> 22; h2 *= m;
	h1 ^= h2 >> 17; h1 *= m;
	h2 ^= h1 >> 19; h2 *= m;

	uint64_t h = h1;

	h = (h << 32) | h2;

	return h;
}

 //convert a string of nucleotide bases into bit array.
void Init_Read(string &seq,struct read_t & read)
{
	read.readLen=strlen(seq.c_str());
	int Read_arr_sz=read.readLen/32+1;
	int rem=read.readLen%32;
	if(rem==0)
	{Read_arr_sz--;}

	str2bitsarr(seq.c_str(),(int)seq.size(),read.read_bits,Read_arr_sz);

}




void Init_Ref_Read(string &seq, struct ref_read_t & read)
{
	read.readLen = strlen(seq.c_str());
	int Read_arr_sz = read.readLen / 32 + 1;
	int rem = read.readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	str2bitsarr(seq.c_str(), (int)seq.size(), read.read_bits, Read_arr_sz);

}

void reverse_complement_str(string & str)
{
	if (str.size() == 0)
	{
		return;
	}

	reverse(str.begin(), str.end());
	for (size_t i = 0; i != str.size(); ++i)
	{
		switch (str[i])
		{
		case 'A':case 'a':
			str[i] = 'T';
			break;
		case 'C':case 'c':
			str[i] = 'G';
			break;
		case 'G':case 'g':
			str[i] = 'C';
			break;
		case 'T':case 't':
			str[i] = 'A';
			break;
		case 'N':case 'n':
			str[i] = 'N';
			break;
		case '-':
			str[i] = '-';
			break;
		default:
			cout << "error: complement_str" << str[i] << endl;
			return;


		}

	}
}

//get the complement of a string of nucleotide bases
void complement_str(string & str)
{
	for (size_t i=0;i!=str.size();++i)
	{
		switch (str[i])
		{
		case 'A':case 'a':
			str[i]='T';
			break;
		case 'C':case 'c':
			str[i]='G';
			break;
		case 'G':case 'g':
			str[i]='C';
			break;
		case 'T':case 't':
			str[i]='A';
			break;
		case 'N':case 'n':
			str[i]='N';
			break;
		case '-':
			str[i]='-';
			break;
		default:
			cout<<"error: complement_str"<<str[i]<<endl;
			return;


		}

	}
}

//convert a bit array into a string
char * bitsarr2str(uint64_t* b_seq, int len,char * c_str,int arr_sz)
{

	int tot_bits=arr_sz*64;
	//char *c_str;
	//c_str=(char*) malloc(sizeof(char)*(len+1));
	//#pragma omp parallel for
	for (int i=0;i<len;++i)
	{
		uint64_t temp,temp2[100];/////////////////////////
		for (int k=0;k<arr_sz;++k)
		{
			temp2[k]=b_seq[k];
		}
		L_shift_NB(temp2,tot_bits-(len-i)*2,arr_sz);
		R_shift_NB(temp2,tot_bits-2,arr_sz);
		//uint64_t temp=(b_seq<<(64-(len-i)*2))>>62;
		temp=temp2[arr_sz-1];
		switch(temp)
		{
		case 0:
			c_str[i]='A';
			break;
		case 1:
			c_str[i]='C';
			break;
		case 2:
			c_str[i]='G';
			break;
		case 3:
			c_str[i]='T';
			break;

		}
	}
	c_str[len]='\0';
	return c_str;
}



#endif
