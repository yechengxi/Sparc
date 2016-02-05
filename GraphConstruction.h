#ifndef __GRAPH_CONSTRUCTION_H
#define __GRAPH_CONSTRUCTION_H


#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"
#include "Align.h"
#include "BasicDataStructure.h"
using namespace std;


void Consensus_Kmer_Graph_Construction(struct ref_read_t *read, struct backbone_info *backbone_info, int K_size)
{
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;

	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	int Kmer_arr_sz = K_size / 32 + 1;
	int rem1 = K_size % 32;
	if (rem1 == 0)
	{
		Kmer_arr_sz--;
	}

	int tot_bits = Read_arr_sz * 64;

	uint64_t seq;
	//check the read to see if there is a saved kmer in the hashtable or bloom filter
	consensus_node *previous_node = NULL, *current_node = NULL;
	backbone_info->node_vec.resize(OverlappingKmers);

	for (int j = 0; j < OverlappingKmers; j++)
	{

		previous_node = current_node;
		get_sub_arr(read->read_bits, read->readLen, j, K_size, &seq);

		current_node = (consensus_node*)malloc(sizeof(consensus_node));
		backbone_info->node_vec[j] = current_node;
		memset(current_node, 0, sizeof(consensus_node));
		memcpy(&(current_node->kmer), &seq, sizeof(uint64_t));
		current_node->in_backbone = 1;
		current_node->cov++;
		backbone_info->n_nodes++;
		if (j >= 1)
		{
			//left edge,right edge
			previous_node->right = (consensus_edge_node*)malloc(sizeof(consensus_edge_node));
			memset(previous_node->right, 0, sizeof(consensus_edge_node));
			previous_node->right->node_ptr = current_node;
			//previous_node->right->edge_cov++;
			previous_node->coord = j-1;
			backbone_info->n_edges++;
			current_node->left = (consensus_edge_node*)malloc(sizeof(consensus_edge_node));
			memset(current_node->left, 0, sizeof(consensus_edge_node));
			current_node->left->node_ptr = previous_node;
			//current_node->left->edge_cov++;
			backbone_info->n_edges++;
		}


	}




}


void Consensus_Sparse_Kmer_Graph_Construction(struct ref_read_t *read, struct backbone_info *backbone_info, int K_size)
{
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;
	int gap = backbone_info->gap;
	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}
	int Kmer_arr_sz = K_size / 32 + 1;
	int rem1 = K_size % 32;
	if (rem1 == 0)
	{
		Kmer_arr_sz--;
	}

	int tot_bits = Read_arr_sz * 64;
	uint64_t seq;


	sparse_consensus_node *previous_node = NULL, *current_node = NULL;
	backbone_info->sparse_node_vec.resize(OverlappingKmers);
	for (int j = 0; j < OverlappingKmers; j++)
	{
		backbone_info->sparse_node_vec[j] = NULL;
	}
	int last_position = -1;
	for (int j = 0; j < OverlappingKmers; j++)
	{
		if (j%gap == 0)
		{
			
			previous_node = current_node;
			get_sub_arr(read->read_bits, read->readLen, j, K_size, &seq);
			//char test[100];
			//bitsarr2str(&seq, K_size, test, 1);
			current_node = (sparse_consensus_node*)malloc(sizeof(sparse_consensus_node));
			backbone_info->sparse_node_vec[j] = current_node;
			memset(current_node, 0, sizeof(sparse_consensus_node));
			current_node->kmer = (uint32_t)seq;
			current_node->in_backbone = 1;
			current_node->cov++;
			current_node->coord = j;
			backbone_info->n_nodes++;

			//left edge,right edge
			uint64_t edge_bits;
			if (j >= 1)
			{
				get_sub_arr(read->read_bits, read->readLen, last_position + K_size, gap, &edge_bits);
				previous_node->right = (sparse_consensus_edge_node*)malloc(sizeof(sparse_consensus_edge_node));
				memset(previous_node->right, 0, sizeof(sparse_consensus_edge_node));
				previous_node->right->node_ptr = current_node;
				previous_node->right->edge = edge_bits;
				previous_node->right->len = gap;
				//previous_node->right->edge_cov++;//decide whether assign weight to the backbone edge
				backbone_info->n_edges++;

				get_sub_arr(read->read_bits, read->readLen, last_position, gap, &edge_bits);
				current_node->left = (sparse_consensus_edge_node*)malloc(sizeof(sparse_consensus_edge_node));
				memset(current_node->left, 0, sizeof(sparse_consensus_edge_node));
				current_node->coord = j;
				current_node->left->node_ptr = previous_node;
				//current_node->left->edge_cov++;//decide whether assign weight to the backbone edge
				current_node->left->len = gap;
				current_node->left->edge = edge_bits;
				backbone_info->n_edges++;
			}
			last_position = j;
		}




	}




}


void NormalizeAlignment(query_info *query_info)
{

	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2 * seq_sz);
	tAlignedSeq_new.resize(2 * seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int n_char = 0;
	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		if (query_info->qAlignedSeq[i] == query_info->tAlignedSeq[i])
		{
			//qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
			//tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
		}
		else
		{
			if (query_info->qAlignedSeq[i] != '-'&&query_info->tAlignedSeq[i] != '-')
			{
				//tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;

				bool replace = 0;
				for (int j = i + 1; j < query_info->tAlignedSeq.size(); ++j)
				{
					if (query_info->tAlignedSeq[j] != '-')
					{
						if (query_info->tAlignedSeq[j] == query_info->qAlignedSeq[i])
						{
							replace = 1;
							query_info->tAlignedSeq[j] = '-';
							//tAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
							tAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
						}
						break;
					}
				}
				if (replace == 0)
				{
					//tAlignedSeq_new.push_back('-');
					tAlignedSeq_new[n_char] = '-';
				}
				n_char--;

				//qAlignedSeq_new.push_back('-');
				//qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
				qAlignedSeq_new[n_char] = '-';
				n_char++;
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];

				n_char++;

			}
			else
			{
				if (query_info->qAlignedSeq[i] == '-')
				{

					for (int j = i + 1; j < query_info->qAlignedSeq.size(); ++j)
					{
						if (query_info->qAlignedSeq[j] != '-')
						{
							if (query_info->qAlignedSeq[j] == query_info->tAlignedSeq[i])
							{
								query_info->qAlignedSeq[i] = query_info->qAlignedSeq[j];
								query_info->qAlignedSeq[j] = '-';

							}
							break;
						}
					}

				}
				else
				{

					for (int j = i + 1; j < query_info->qAlignedSeq.size(); ++j)
					{
						if (query_info->tAlignedSeq[j] != '-')
						{
							if (query_info->tAlignedSeq[j] == query_info->qAlignedSeq[i])
							{
								query_info->tAlignedSeq[i] = query_info->tAlignedSeq[j];
								query_info->tAlignedSeq[j] = '-';

							}
							break;
						}
					}
				}

				//qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
				//tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;
			}
		}
	}
	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);
	query_info->qAlignedSeq = qAlignedSeq_new;
	query_info->tAlignedSeq = tAlignedSeq_new;




}


void PatchGaps(query_info *query_info)
{
	bool Align = 1;
	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2* seq_sz);
	tAlignedSeq_new.resize(2* seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int MinGapSize = 2, K_size = query_info->Patch_K, SearchDepth = query_info->Patch_D;
	int n_char = 0;
	
	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		bool Patch = 0;
			
		if (query_info->qAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j <i+ MinGapSize; ++j)
			{
				if ((query_info->qAlignedSeq[j] != '-') || (j+1==seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}
		
		if (!Patch&&Align)
		{
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
		}

		if (Patch)
		{
			int  shift = -1;
			map<string, int> target_position, query_position;
			vector<int> target_index, query_index;
			target_index.resize(300);
			query_index.resize(300);
			string target_crop, query_crop;
			int target_bases = 0, query_bases = 0;
			for (int j = i; j < seq_sz; ++j)
			{
				if (query_info->tAlignedSeq[j] != '-')
				{
					target_crop.push_back(query_info->tAlignedSeq[j]);
					target_index[target_bases] = j;
					target_bases++;
						
				}
				if (target_bases >= SearchDepth + K_size)
				{
					break;
				}
			}

			for (int j = i; j < seq_sz; ++j)
			{
				if (query_info->qAlignedSeq[j] != '-')
				{
					query_crop.push_back(query_info->qAlignedSeq[j]);
					query_index[query_bases] = j;
					query_bases++;
				}
				if (query_bases >= SearchDepth + K_size)
				{
					break;
				}
			}



			if (!Align)
			{
				shift = MinGapSize;
				for (int j = 0; j + K_size < target_crop.size(); ++j)
				{
					string kmer = target_crop.substr(j, K_size);
					if (target_position.count(kmer) == 0)
					{
						target_position[kmer] = j;
					}
				}
				int target_matched = -1, query_matched = -1;
				for (int j = 0; j + K_size < query_crop.size(); ++j)
				{
					string kmer = query_crop.substr(j, K_size);
					if (target_position.count(kmer))
					{
						target_matched = target_position[kmer];
						query_matched = j;
						break;
					}
				}


				if (target_matched >= 0)
				{
					bool Debug = 0;
					if (Debug)
					{
						cout << "before: " << endl;
						cout << query_info->tAlignedSeq.substr(i, 100) << endl;
						cout << query_info->qAlignedSeq.substr(i, 100) << endl;
					}
					//clear the query bases
					for (int j = 0; j < query_matched + K_size; ++j)
					{
						query_info->qAlignedSeq[query_index[j]] = '-';
					}
					for (int j = 0; j < K_size; ++j)
					{
						query_info->qAlignedSeq[target_index[target_matched + j]] = query_info->tAlignedSeq[target_index[target_matched + j]];

					}

					if (Debug)
					{
						cout << "after: " << endl;
						cout << query_info->tAlignedSeq.substr(i, 100) << endl;
						cout << query_info->qAlignedSeq.substr(i, 100) << endl;
					}


				}

			}
			else
			{

				int match = 100, mismatch = -100, gap_cost = -10, band_width = 50;
				string A_aln, B_aln;
				struct aln_t aln_t;
				int score = 0;
				char qry_char[300];
				char ref_char[300];
				strcpy(ref_char, target_crop.c_str());
				strcpy(qry_char, query_crop.c_str());

				score = GlobalAlign(&aln_t, qry_char, ref_char, match, mismatch, gap_cost, band_width);
				printAlign(&aln_t, qry_char, ref_char, A_aln, B_aln);

				bool Debug = 0;

				if (Debug)
				{
					cout << i << endl;
				
					if (i == 27)
					{
						cout << "";
					}
					cout << "before: " << endl;
					cout << query_info->tAlignedSeq.substr(i, 70) << endl;
					cout << query_info->qAlignedSeq.substr(i, 70) << endl;
				}
				if (Debug)
				{
					cout << "after: " << endl;
					
					cout << B_aln << endl;
					cout << A_aln << endl;
				}
				int A_cnt = 0, B_cnt = 0;
				for (int jj = 0; jj != A_aln.size(); ++jj)
				{
					if (A_aln[jj] != '-')
					{
						A_cnt++;
					}
				}

				for (int jj = 0; jj != B_aln.size(); ++jj)
				{
					if (B_aln[jj] != '-')
					{
						B_cnt++;
					}
				}


				int cnt = 0;
				int match_position = -1;
				for (int jj = (int)A_aln.size() - 1; jj >= 0; --jj)
				{
					if (A_aln[jj] == B_aln[jj])
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt == K_size)
					{
						match_position = jj + K_size - 1;
						break;
					}
				}
				//append align

				if (match_position == -1)
				{
					shift = MinGapSize - 1;
					qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
					tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
					n_char++;

					i=i + shift - 1;
					continue;
				}
				for (int jj = 0; jj <= match_position; ++jj)
				{
					qAlignedSeq_new[n_char] = A_aln[jj];
					tAlignedSeq_new[n_char] = B_aln[jj];
					n_char++;

				}
				


				int target_shift = 0, query_shift = 0;
				for (int jj = 0; jj <= match_position; ++jj)
				{
					if (B_aln[jj] != '-')
					{
						target_shift++;
					}
					if (A_aln[jj] != '-')
					{
						query_shift++;
					}
				}

				
				//clear the query bases
				int q_base = 0;
				for (int j = i; j < seq_sz; ++j)
				{
					if (query_info->qAlignedSeq[j]!='-')
					{
						q_base++;
						query_info->qAlignedSeq[j] = '-';
					}
					if (q_base >= query_shift)
					{
						break;
					}
				}
				
				int t_base = 0;
				for (int j = i; j < seq_sz; ++j)
				{

					if (query_info->tAlignedSeq[j] != '-')
					{
						t_base++;
						query_info->tAlignedSeq[j] = '-';
					}
					if (t_base >= target_shift)
					{
						shift = j-i;
						break;
					}
				}


			}

			i = i + shift - 1;


		}

		
	}

	if (!Align)
	{
		for (int i = 0; i < seq_sz; ++i)
		{
			if (query_info->tAlignedSeq[i] != '-' || query_info->qAlignedSeq[i] != '-')
			{
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;
			}

		}
	}
	

	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);
	query_info->qAlignedSeq = qAlignedSeq_new;

	query_info->tAlignedSeq = tAlignedSeq_new;
}



void FillGaps(query_info *query_info)
{
	bool Align = 1;
	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2 * seq_sz);
	tAlignedSeq_new.resize(2 * seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int MinGapSize = query_info->Patch_G, K_size = query_info->Patch_K, SearchDepth = query_info->Patch_D;
	int n_char = 0;

	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		bool Patch = 0;
		//detect gaps in query
		if (query_info->qAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j <i + MinGapSize; ++j)
			{
				if ((query_info->qAlignedSeq[j] != '-') || (j + 1 == seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}

		//detect gaps in target
		if (query_info->tAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j <i + MinGapSize; ++j)
			{
				if ((query_info->tAlignedSeq[j] != '-') || (j + 1 == seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}

		if (!Patch)
		{
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
			continue;
		}
		
		
		int  shift = -1;
		vector<int> target_index, query_index;
		target_index.resize(300);
		query_index.resize(300);
		string target_crop, query_crop;
		target_crop.resize(300);
		query_crop.resize(300);
		int target_bases = 0, query_bases = 0;
		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->tAlignedSeq[j] != '-')
			{
				target_crop[target_bases]=query_info->tAlignedSeq[j];
				target_index[target_bases] = j;
				target_bases++;

			}
			if (target_bases >= SearchDepth + K_size)
			{
				break;
			}
		}

		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				query_crop[query_bases]=query_info->qAlignedSeq[j];
				query_index[query_bases] = j;
				query_bases++;
			}
			if (query_bases >= SearchDepth + K_size)
			{
				break;
			}
		}
		target_crop.resize(target_bases);
		query_crop.resize(query_bases);

		int match = 100, mismatch = -100, gap_cost = -10, band_width = 50;
		string Query_aln, Target_aln;
		struct aln_t aln_t;
		int score = 0;
		char qry_char[300];
		char ref_char[300];
		strcpy(ref_char, target_crop.c_str());
		strcpy(qry_char, query_crop.c_str());

		score = GlobalAlign(&aln_t, qry_char, ref_char, match, mismatch, gap_cost, band_width);
		printAlign(&aln_t, qry_char, ref_char, Query_aln, Target_aln);

		bool Debug = 0;
		if (Debug)
		{
			cout << i << endl;

			
			cout << "before: " << endl;
			cout << query_info->tAlignedSeq.substr(i, 70) << endl;
			cout << query_info->qAlignedSeq.substr(i, 70) << endl;
		}
		if (Debug)
		{
			cout << "after: " << endl;

			cout << Target_aln << endl;
			cout << Query_aln << endl;
		}
		
		int cnt = 0;
		int match_position = -1;
		for (int jj = (int)Query_aln.size() - 1; jj >= 0; --jj)
		{
			if (Query_aln[jj] == Target_aln[jj])
			{
				cnt++;
			}
			else
			{
				cnt = 0;
			}
			if (cnt == K_size)
			{
				match_position = jj + K_size - 1;
				break;
			}
		}
		//append align
		
		int shift_end;
		shift_end = i + 1;
		if (query_crop.size() == target_crop.size() && target_crop.size() == SearchDepth + K_size)
		{
			shift_end = i + SearchDepth;// min(target_index[SearchDepth + K_size - 1], query_index[SearchDepth + K_size - 1]);
		}

		if (match_position == -1)
		{
			int j;
			for (j = i; j < shift_end ; ++j)
			{
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[j];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[j];
				n_char++;
				
			}
			i=j-1;

			continue;
		}



		for (int jj = 0; jj <= match_position; ++jj)
		{
			qAlignedSeq_new[n_char] = Query_aln[jj];
			tAlignedSeq_new[n_char] = Target_aln[jj];
			n_char++;

		}



		int target_shift = 0, query_shift = 0;
		for (int jj = 0; jj <= match_position; ++jj)
		{
			if (Target_aln[jj] != '-')
			{
				target_shift++;
			}
			if (Query_aln[jj] != '-')
			{
				query_shift++;
			}
		}


		//clear the aligned bases
		int q_base = 0;
		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				q_base++;
				query_info->qAlignedSeq[j] = '-';
			}
			if (q_base >= query_shift)
			{
				shift = query_index[q_base - 1];
				break;
				
			}
		}

		int t_base = 0;
		for (int j = i; j < seq_sz; ++j)
		{

			if (query_info->tAlignedSeq[j] != '-')
			{
				t_base++;
				query_info->tAlignedSeq[j] = '-';
			}
			if (t_base >= target_shift)
			{
				shift = min(shift,target_index[t_base - 1]);
				break;
			}
		}
		
		i = shift ;
	}



	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);
	query_info->qAlignedSeq = qAlignedSeq_new;
	query_info->tAlignedSeq = tAlignedSeq_new;
}


void Add_Path_To_Backbone(struct backbone_info *backbone_info, struct query_info *query_info, int K_size)
{
	ofstream o_report_align;
	bool boost_edges = 0;
	string ContigPrefix = backbone_info->ContigPrefix;
	if (ContigPrefix.size() > 0)
	{
		if (query_info->qName.substr(0, ContigPrefix.size()) == ContigPrefix)
		{
			boost_edges = 1;
		}
	}
	bool DEBUG = 0;
	if (query_info->tStrand == '-')
	{
		reverse_complement_str(query_info->qAlignedSeq);
		reverse_complement_str(query_info->tAlignedSeq);
		reverse(query_info->matchPattern.begin(), query_info->matchPattern.end());

	}

	NormalizeAlignment(query_info);
	if (query_info->Patch)
	{
		PatchGaps(query_info);

	}
	else
	{
		if (query_info->Fill)
		{
			FillGaps(query_info);

		}

	}
	//		FillGaps(query_info);


	for (int i = 0; i < query_info->tAlignedSeq.size(); ++i)
	{
		if (query_info->tAlignedSeq[i] != '-')
		{
			if (query_info->qAlignedSeq[i] != '-')
			{
				backbone_info->ref_matched++;
			}
			else
			{
				backbone_info->ref_mismatch++;
			}
		}
	}


	if ((query_info->report_e > query_info->report_b) && (query_info->tStart<query_info->report_b) && (query_info->tEnd>query_info->report_e))
	{
		o_report_align.open("subalign.txt", ios_base::app);
		int target_position = query_info->tStart;
		string target_crop, query_crop;
		for (int k = 0; k < query_info->tAlignedSeq.size(); ++k)
		{
			if (query_info->tAlignedSeq[k] != '-')
			{
				target_position++;
			}
			if (target_position>query_info->report_b)
			{
				target_crop.push_back(query_info->tAlignedSeq[k]);
				query_crop.push_back(query_info->qAlignedSeq[k]);
				//
			}

			if (target_position>query_info->report_e)
			{
				o_report_align << ">target_" << query_info->read_idx << endl << target_crop << endl;//
				o_report_align << ">query_" << query_info->read_idx << endl << query_crop << endl;//
				break;
			}
		}
	}


	string backbone_seg = backbone_info->backbone.substr(query_info->tStart, query_info->tEnd - query_info->tStart);
	//cout << backbone_seg << endl;
	//cout << target_seg << endl;
	int target_position = query_info->tStart;
	consensus_node *previous_node = NULL, *current_node = NULL;
	int MatchPosition = -100, PreviousMatchPosition = -100;

	for (int i = 0; i + K_size <= query_info->qAlignedSeq.size(); ++i)
	{

		if (query_info->qAlignedSeq[i] == '-')
		{
			MatchPosition = -1;
			PreviousMatchPosition = -1;
			if (query_info->tAlignedSeq[i] != '-')
			{
				target_position++;
			}
			continue;
		}
		PreviousMatchPosition = MatchPosition;
		MatchPosition = target_position;
		for (int j = 0; j < K_size; ++j)
		{
			if (query_info->qAlignedSeq[i + j] != query_info->tAlignedSeq[i + j])
			{
				MatchPosition = -1;
				break;
			}
		}

		string KmerStr;
		KmerStr.resize(K_size);
		int KmerSize = 0;
		for (int j = i; j < query_info->qAlignedSeq.size(); ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				KmerStr[KmerSize]=query_info->qAlignedSeq[j];
				KmerSize++;
			}
			
			if (KmerSize == K_size)
			{
				break;
			}
		}

		if (KmerSize != K_size)
		{
			break;
		}
		uint64_t seq = 0;
		str2bitsarr(KmerStr.c_str(), K_size, &seq, 1);


		if (MatchPosition >= 0)
		{
			query_info->n_exist++;
			current_node = backbone_info->node_vec[MatchPosition];

			if (previous_node != NULL)
			{
				current_node->cov++;
			}

			//link previous node to the current one

			if (previous_node != NULL)
			{

				consensus_edge_node **edge_p2p = &(previous_node->right);

				while (*edge_p2p != NULL)
				{

					if ((*edge_p2p)->node_ptr == current_node)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}


					edge_p2p = &((*edge_p2p)->nxt_edge);
				}

				if (*edge_p2p == NULL)
				{

					(*edge_p2p) = (consensus_edge_node *)malloc(sizeof(consensus_edge_node));
					memset((*edge_p2p), 0, sizeof(consensus_edge_node));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = current_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}

				edge_p2p = &(current_node->left);
				while (*edge_p2p != NULL)
				{
					if ((*edge_p2p)->node_ptr == previous_node)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}
				if (*edge_p2p == NULL)
				{
					(*edge_p2p) = (consensus_edge_node *)malloc(sizeof(consensus_edge_node));
					memset((*edge_p2p), 0, sizeof(consensus_edge_node));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = previous_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}


			}

		}
		else
		{

			//link previous node to the current one

			if (previous_node == NULL)
			{
				current_node = (consensus_node *)malloc(sizeof(consensus_node));
				memset(current_node, 0, sizeof(consensus_node));

				current_node->kmer = (uint32_t)seq;
				current_node->cov++;
				backbone_info->n_nodes++;
				previous_node = current_node;
				query_info->n_new++;
				//cout << i << ", ";
			}
			else
			{

				consensus_edge_node **edge_p2p = &(previous_node->right);
				while (*edge_p2p != NULL)
				{
					if ((*edge_p2p)->node_ptr->kmer == (uint32_t)seq && (*edge_p2p)->node_ptr->in_backbone == 0)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}
						current_node = (*edge_p2p)->node_ptr;
						current_node->cov++;

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}

				if (*edge_p2p == NULL)
				{
					(*edge_p2p) = (consensus_edge_node *)malloc(sizeof(consensus_edge_node));
					memset((*edge_p2p), 0, sizeof(consensus_edge_node));
					(*edge_p2p)->nxt_edge = NULL;
					current_node = (consensus_node *)malloc(sizeof(consensus_node));
					query_info->n_new++;
					memset(current_node, 0, sizeof(consensus_node));
					current_node->kmer = (uint32_t)seq;
					current_node->cov++;
					(*edge_p2p)->node_ptr = current_node;

					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
					backbone_info->n_nodes++;
					//cout << i <<", ";


				}
				//link current node to previous
				edge_p2p = &(current_node->left);
				while (*edge_p2p != NULL)
				{
					if ((*edge_p2p)->node_ptr->kmer == previous_node->kmer)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}
				if (*edge_p2p == NULL)
				{
					(*edge_p2p) = (consensus_edge_node *)malloc(sizeof(consensus_edge_node));
					memset((*edge_p2p), 0, sizeof(consensus_edge_node));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = previous_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}

			}
		}

		previous_node = current_node;

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

	}


}

void Add_Path_To_Backbone_Sparse(struct backbone_info *backbone_info, struct query_info *query_info, int K_size)
{
	//bool TEST = 1;
	int gap = backbone_info->gap;
	ofstream o_report_align;
	bool boost_edges = 0;
	string ContigPrefix = backbone_info->ContigPrefix;
	if (ContigPrefix.size() > 0)
	{
		if (query_info->qName.substr(0, ContigPrefix.size()) == ContigPrefix)
		{
			boost_edges = 1;
		}
	}

	if (query_info->tStrand == '-')
	{
		reverse_complement_str(query_info->qAlignedSeq);
		reverse_complement_str(query_info->tAlignedSeq);
		reverse(query_info->matchPattern.begin(), query_info->matchPattern.end());

	}
	if (query_info->read_idx == 327)
	{
		//cout << "16" << endl;
	}
	int bases_q1 = 0, bases_t1 = 0, bases_q2 = 0, bases_t2 = 0;
	string t1, q1;
	for (int i = 0; i < query_info->tAlignedSeq.size(); ++i)
	{
		if (query_info->qAlignedSeq[i] != '-')
		{
			bases_q1++;
		//	q1.push_back(query_info->qAlignedSeq[i]);
		}
		if (query_info->tAlignedSeq[i] != '-')
		{
			bases_t1++;
		//	t1.push_back(query_info->tAlignedSeq[i]);
		}
	}
	NormalizeAlignment(query_info);
	if (query_info->Patch)
	{
		PatchGaps(query_info);

	}
	else
	{
		if (query_info->Fill)
		{
			FillGaps(query_info);

		}

	}
	string t2, q2;
	for (int i = 0; i < query_info->tAlignedSeq.size(); ++i)
	{
		if (query_info->qAlignedSeq[i] != '-')
		{
			bases_q2++;
		//	q2.push_back(query_info->qAlignedSeq[i]);
		}
		if (query_info->tAlignedSeq[i] != '-')
		{
			bases_t2++;
		//	t2.push_back(query_info->tAlignedSeq[i]);
		}
	}
	if ((bases_q1 != bases_q2) || (bases_t1 != bases_t2))
	{
		//if (q1.size() != q2.size())
		//{
		//	cout << q1.substr(0,200) << endl;
		//	cout << q2.substr(0, 200) << endl;
		//}
		//cout << "Normalization error." << endl;
	}
	for (int i = 0; i < query_info->tAlignedSeq.size(); ++i)
	{
		if (query_info->tAlignedSeq[i] != '-')
		{
			if (query_info->qAlignedSeq[i] != '-')
			{
				backbone_info->ref_matched++;
			}
			else
			{
				backbone_info->ref_mismatch++;
			}
		}
	}

	if ((query_info->report_e > query_info->report_b) && (query_info->tStart<query_info->report_b) && (query_info->tEnd>query_info->report_e))
	{
		o_report_align.open("subalign.txt", ios_base::app);
		int target_position = query_info->tStart;
		string target_crop, query_crop;
		for (int k = 0; k < query_info->tAlignedSeq.size(); ++k)
		{
			if (query_info->tAlignedSeq[k] != '-')
			{
				target_position++;
			}
			if (target_position>query_info->report_b)
			{
				target_crop.push_back(query_info->tAlignedSeq[k]);
				query_crop.push_back(query_info->qAlignedSeq[k]);
				//
			}

			if (target_position>query_info->report_e)
			{
				o_report_align << ">target_" << query_info->read_idx << endl << target_crop << endl;//
				o_report_align << ">query_" << query_info->read_idx << endl << query_crop << endl;//
				break;
			}
		}
	}


	int target_position = query_info->tStart;
	sparse_consensus_node *previous_node = NULL, *current_node = NULL;
	int FirstTargetPosition = -100;
	int qStart = -1;
	vector<bool> matched;
	matched.resize(query_info->qAlignedSeq.size());
	vector<int> Align2Backbone, Align2Query, Query2Align;

	Align2Backbone.resize(query_info->qAlignedSeq.size());
	Align2Query.resize(query_info->qAlignedSeq.size());
	Query2Align.resize(query_info->qAlignedSeq.size());
	string processed_query, processed_target;
	processed_query.resize(query_info->qAlignedSeq.size());
	processed_target.resize(query_info->qAlignedSeq.size());
	int n_t = 0, n_q = 0;

	for (int i = 0; i < query_info->qAlignedSeq.size(); ++i)
	{
		if (query_info->qAlignedSeq[i] != '-')
		{
			processed_query[n_q] = (query_info->qAlignedSeq[i]);
			n_q++;
		}
		if (query_info->tAlignedSeq[i] != '-')
		{
			processed_target[n_t] = (query_info->tAlignedSeq[i]);
			n_t++;
		}
	}
	processed_query.resize(n_q);
	processed_target.resize(n_t);

	int base_cnt = 0;
	for (int i = 0; i < query_info->qAlignedSeq.size(); ++i)
	{
		Align2Query[i] = base_cnt;
		if (query_info->qAlignedSeq[i] != '-')
		{
			Query2Align[base_cnt] = i;
			base_cnt++;
		}
	}

	base_cnt = query_info->tStart;
	for (int i = 0; i + K_size <= query_info->qAlignedSeq.size(); ++i)
	{
		Align2Backbone[i] = base_cnt;
		if (query_info->tAlignedSeq[i] != '-')
		{
			base_cnt++;
			/*
			if (TEST)
			{
				if (backbone_info->sparse_node_vec[Align2Backbone[i]] != NULL)
				{
					uint64_t kmer = backbone_info->sparse_node_vec[Align2Backbone[i]]->kmer;
					char kmer_str[100];
					bitsarr2str(&kmer, K_size, kmer_str, 1);
					if (kmer_str[0] != query_info->tAlignedSeq[i])
					{
						cout << "Align2Backbone error." << endl;
					}
				}

			}
			*/
		}

		
	}


	for (int i = 0; i  < query_info->qAlignedSeq.size(); ++i)
	{
		matched[i] = 0;
	}

	for (int i = 0; i + K_size <= query_info->tAlignedSeq.size(); ++i)
	{
		bool match = 1;
		int match_bases = 0;
		for (int j = 0; j < query_info->tAlignedSeq.size(); ++j)
		{
			if (query_info->qAlignedSeq[i + j] == query_info->tAlignedSeq[i + j] )
			{
				if (query_info->qAlignedSeq[i + j] != '-')
				{
					match_bases++;
				}
				if (match_bases == K_size)
				{
					break;
				}
				
			}
			else
			{
				match = 0;
				break;
			}

		}

		matched[i] = match;
		if (qStart < 0 && match)
		{
			qStart = i;
		}
		
	}
	

	if (qStart < 0)
	{
		return;
	}

	if (backbone_info->sparse_node_vec[Align2Backbone[qStart]] != NULL)
	{
		current_node = backbone_info->sparse_node_vec[Align2Backbone[qStart]];
	}
	else
	{
		string KmerStr = processed_query.substr(Align2Query[qStart], K_size);
		
		if (KmerStr.size() != K_size)
		{
			return;
		}
		uint64_t seq = 0;
		str2bitsarr(KmerStr.c_str(), K_size, &seq, 1);
		current_node = (sparse_consensus_node *)malloc(sizeof(sparse_consensus_node));
		memset(current_node, 0, sizeof(sparse_consensus_node));
		current_node->kmer = (uint32_t)seq;
		backbone_info->n_nodes++;
		backbone_info->sparse_node_vec[Align2Backbone[qStart]] = current_node;
		current_node->in_backbone = 1;
	}

	for (int i = qStart; i < query_info->qAlignedSeq.size();)
	{
		//explore the edges of the current node
		
		if (query_info->read_idx == 8 && i >= 167)
		{
			cout << "";
		}
		if (query_info->qAlignedSeq[i] == '-')
		{
			i++;// must skip to prmise the correctness of the algorithm
			continue;
		}
		current_node->cov++;//increase the coverage
		string EdgeStr;
		EdgeStr = processed_query.substr(Align2Query[i] + K_size, gap);
		

		if (EdgeStr.size() == 0)
		{
			
			return;
		}
		struct sparse_consensus_edge_node* current_edge = current_node->right;
		bool edge_exists = 0;
		while (current_edge != NULL)
		{
			if (EdgeStr.size() < current_edge->len)
			{
				current_edge = current_edge->nxt_edge;
				continue;
			}
			uint64_t edge_bits;
			str2bitsarr(EdgeStr.c_str(), current_edge->len, &edge_bits, 1);
			
			if (edge_bits == current_edge->edge)
			{
				//cout << "Path exists.";
				
				bool edge_matched = 0;

				int align_idx = Query2Align[Align2Query[i] + current_edge->len];
				if (matched[align_idx] && backbone_info->sparse_node_vec[Align2Backbone[align_idx]] == current_edge->node_ptr)
				{
					edge_matched = 1;//backbone matched, so move on
				}
				if (!matched[align_idx] && !current_edge->node_ptr->in_backbone)
				{
					edge_matched = 1;//off-backbone branch matched, so move on
				}
				
				if (edge_matched)
				{
					current_edge->edge_cov++;
					if (boost_edges)
					{
						current_edge->edge_cov += backbone_info->boost;
					}
					
					edge_exists = 1;
					query_info->n_exist++;
					current_node = current_edge->node_ptr;
					
					string LeftEdgeStr = processed_query.substr(Align2Query[i], current_edge->len);
					
					str2bitsarr(LeftEdgeStr.c_str(), current_edge->len, &edge_bits, 1);
					struct sparse_consensus_edge_node* left_edge = current_node->left;
					while (left_edge != NULL)
					{
						if (edge_bits == left_edge->edge&& LeftEdgeStr.size() == current_edge->len)
						{
							break;
						}

						left_edge = left_edge->nxt_edge;
					}
					if (left_edge != NULL)
					{
						left_edge->edge_cov++;
					}
					else
					{
						cout << "Error." << endl;
						return;
						cout << "left edge error." << endl;

					}

					

					i = align_idx;
					break;

				}


			}

			current_edge = current_edge->nxt_edge;
		}

		if (edge_exists)
		{
			continue;
		}


		//There is a backbone match but was not recorded, or there is no match.

		/////non-matching, so pick the best stretch
		// case 1, within g bases, there is a match with the backbone -- put the node into the backbone nodes
		// otherwise, skip g bases, and create a new node.

		current_node->cov++;
		bool backbone_match = 0;

		int jj;//jj is the position in the processed query read
		int ii = Align2Query[i];// these two lines are very tricky, 
		int align_idx = -1;
		for (jj = ii + 1; jj <= ii + gap; ++jj)
		{
			//	cout << jj << endl;
			align_idx = Query2Align[jj];
			if (matched[align_idx] )
			{
				backbone_match = 1;
				break;
			}
			if (jj + K_size >= processed_query.size())
			{
				break;
			}
		}
		jj = min(jj, ii + gap);

		string temp_str;
		temp_str = processed_query.substr(ii, jj - ii + K_size);

		if (temp_str.size() != jj - ii + K_size)
		{
			cout << "bug. skipped" << endl;
			return;
		}

		sparse_consensus_node* next_node;
		if (backbone_match&& backbone_info->sparse_node_vec[Align2Backbone[align_idx]]!=NULL)
		{
			//add edges to the existing backbone node.
			next_node = backbone_info->sparse_node_vec[Align2Backbone[align_idx]];
			if (next_node == current_node)
			{
				cout << "bug caught, linked to itself." << endl;
				cout << query_info->read_idx << " " << i << endl;
				i++;
				return;
			}


		}
		else
		{
			//allocate a new node and new edges 
			string kmer_str = temp_str.substr(temp_str.size() - K_size, K_size);
			next_node = (sparse_consensus_node*)malloc(sizeof(sparse_consensus_node));
			memset(next_node, 0, sizeof(sparse_consensus_node));
			query_info->n_new++;
			uint64_t seq;
			str2bitsarr(kmer_str.c_str(), K_size, &seq, 1);
			next_node->kmer = (uint32_t)seq;
			if (backbone_match)
			{
				backbone_info->sparse_node_vec[Align2Backbone[align_idx]] = next_node;
				next_node->in_backbone = 1;
				next_node->coord = Align2Backbone[align_idx];
			}
			backbone_info->n_nodes++;
			//cout << i << ", ";
			//cout << Query2Align[Align2Query[i]+1] << ", ";
		}


		//append edges, quite tricky
		string edge_str;
		edge_str = temp_str.substr(K_size, temp_str.size());
		if (edge_str.size() > 1)
		{
			cout << "";
		}

		uint64_t edge_bits;
		str2bitsarr(edge_str.c_str(), edge_str.size(), &edge_bits, 1);

		sparse_consensus_edge_node **edge_p2p = &(current_node->right);

		while (*edge_p2p != NULL)
		{

			edge_p2p = &((*edge_p2p)->nxt_edge);
		}

		if (*edge_p2p == NULL)
		{

			(*edge_p2p) = (sparse_consensus_edge_node *)malloc(sizeof(sparse_consensus_edge_node));
			memset((*edge_p2p), 0, sizeof(sparse_consensus_edge_node));
			(*edge_p2p)->edge = (uint32_t)edge_bits;
			(*edge_p2p)->len = edge_str.size();
			(*edge_p2p)->nxt_edge = NULL;
			(*edge_p2p)->node_ptr = next_node;
			(*edge_p2p)->edge_cov++;
			if (boost_edges)
			{
				(*edge_p2p)->edge_cov += backbone_info->boost;
			}
			backbone_info->n_edges++;
		}



		edge_str = temp_str.substr(0, temp_str.size() - K_size);
		str2bitsarr(edge_str.c_str(), edge_str.size(), &edge_bits, 1);

		edge_p2p = &(next_node->left);

		while (*edge_p2p != NULL)
		{

			edge_p2p = &((*edge_p2p)->nxt_edge);
		}

		if (*edge_p2p == NULL)
		{

			(*edge_p2p) = (sparse_consensus_edge_node *)malloc(sizeof(sparse_consensus_edge_node));
			memset((*edge_p2p), 0, sizeof(sparse_consensus_edge_node));
			(*edge_p2p)->edge = (uint32_t)edge_bits;
			(*edge_p2p)->len = edge_str.size();
			(*edge_p2p)->nxt_edge = NULL;
			(*edge_p2p)->node_ptr = current_node;
			(*edge_p2p)->edge_cov++;
			if (boost_edges)
			{
				(*edge_p2p)->edge_cov += backbone_info->boost;
			}
			backbone_info->n_edges++;
		}
		uint64_t kmer_bits = current_node->kmer;

		char kmer_str[100];
		bitsarr2str(&kmer_bits, K_size, kmer_str, 1);
		if ((*edge_p2p)->len == K_size&&strcmp(kmer_str, edge_str.c_str()) != 0)
		{
			cout << "bug" << endl;
			cout << kmer_str << endl;
			cout << edge_str << endl;
		}

		current_node = next_node;
		i = align_idx;
		//i = i + 1;

	}

}




#endif
