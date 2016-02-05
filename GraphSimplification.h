#ifndef __GRAPH_SIMPLIFICATION_H
#define __GRAPH_SIMPLIFICATION_H
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
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
using namespace std;


void MergeNodes(struct backbone_info *backbone_info)
{
	//break the bubble links in the bfs.
	ofstream o_debug("debug_merge.txt");
	int n_merged = 0, n_replacements = 0;
	for (int i = 0; i + 1 <backbone_info->node_vec.size(); ++i)
	{
		//cout << i << " ";
		consensus_edge_node *edge_ptr = backbone_info->node_vec[i]->right;
		map<uint32_t, consensus_node *> node_map;
		map< consensus_node *, int> node_cov;
		uint32_t kmer = edge_ptr->node_ptr->kmer;
		int backbone_edge_cov = backbone_info->node_vec[i]->right->edge_cov;
		//edge_ptr = edge_ptr->nxt_edge;
		int max_cov = 0;
		while (edge_ptr != NULL)
		{
			if (node_map.count(edge_ptr->node_ptr->kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer]] += edge_ptr->edge_cov;
				n_merged++;
				if (edge_ptr->node_ptr->used)
				{
					//cout << "err" << endl;
				}
				edge_ptr->node_ptr->used = 1;

				
			}
			edge_ptr = edge_ptr->nxt_edge;
		}

		edge_ptr = backbone_info->node_vec[i]->right;
		

		map< consensus_node *, int>::iterator it;

		for (it = node_cov.begin(); it != node_cov.end(); ++it)
		{
			if (it->first->kmer == kmer)
			{
				backbone_edge_cov = node_cov[node_map[kmer]];
			}
			else
			{
				if (max_cov < it->second)
				{
					max_cov = it->second;
				}
			}

		}
		
		if (max_cov >= backbone_edge_cov)
		{
			n_replacements++;

			o_debug << "position:" << i << " " << endl;
			for (it = node_cov.begin(); it != node_cov.end(); ++it)
			{
				if (it->first->kmer == kmer)
				{
					o_debug << "bb " << it->first->kmer<<" "<<it->second << endl;
				}
				else
				{
					o_debug << it->first->kmer << " " << it->second << endl;;
				}

			}
		}
		

		

	}


	cout << backbone_info->n_nodes << " nodes." << endl;
	cout << backbone_info->n_edges << " edges." << endl;
	cout << n_merged << " merges." << endl;
	cout << n_replacements << " replacements." << endl;

}














void OutputPathsFromANode(consensus_node * begin_node, string filename, map<consensus_node *, bool> &Visited)
{


	ofstream o_graph(filename.c_str(), ios_base::app);
	map<uint64_t, consensus_node *> node_map;

	list<consensus_node *> node_list;
	node_list.push_back(begin_node);
	if (begin_node->coord == 14240)
	{
		cout << "";
	}
	//o_graph << "\"" << begin_node << "\"" << " [label=\"" << begin_node->kmer.kmer << "(" << begin_node->coord << ")" "\"];" << endl;
	while (node_list.size() > 0)
	{
		consensus_node * current_node = node_list.front();
		node_list.pop_front();
		consensus_edge_node * edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			
			
				
			o_graph << "\"" << current_node << "\"" << " -> " << "\"" << edge_ptr->node_ptr << "\"";
			o_graph << " [label=\"" << edge_ptr->edge_cov << "\"];";
			o_graph << endl;

			

			if ((!edge_ptr->node_ptr->in_backbone) )
			{
				if (Visited.count(edge_ptr->node_ptr) == 0)
				{
					Visited[edge_ptr->node_ptr]=1;
					node_list.push_back(edge_ptr->node_ptr);
				}
				if (edge_ptr->node_ptr->coord > 0)
				{
					o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [color=red,label=\"" << edge_ptr->node_ptr->kmer  << "(" << edge_ptr->node_ptr->coord << ")" "\"];" << endl;

				}
				else
				{
					o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [color=red,label=\"" << edge_ptr->node_ptr->kmer << "\"];" << endl;

				}
				
			}
			else
			{
				
				o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [label=\"" << edge_ptr->node_ptr->kmer << "[" << edge_ptr->node_ptr->cns_coord << "]" << "[" << edge_ptr->node_ptr->coord << "]" "\"];" << endl;
				
			}
				

		
			
			edge_ptr = edge_ptr->nxt_edge;
		}
	}

	o_graph.close();

}


void OutputPathsFromASparseNode(sparse_consensus_node * begin_node, string filename, map<sparse_consensus_node *, bool> &Visited)
{


	ofstream o_graph(filename.c_str(), ios_base::app);

	list<sparse_consensus_node *> node_list;
	node_list.push_back(begin_node);
	if (begin_node->coord == 14240)
	{
		cout << "";
	}
	//o_graph << "\"" << begin_node << "\"" << " [label=\"" << begin_node->kmer.kmer << "(" << begin_node->coord << ")" "\"];" << endl;
	while (node_list.size() > 0)
	{
		sparse_consensus_node * current_node = node_list.front();
		node_list.pop_front();
		sparse_consensus_edge_node * edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			uint64_t edge_bits = edge_ptr->edge;
			char edge_str[100];
			bitsarr2str(&edge_bits,edge_ptr->len,edge_str,1);
			o_graph << "\"" << current_node << "\"" << " -> " << "\"" << edge_ptr->node_ptr << "\"";
			o_graph << " [label=\"" <<edge_str << " (" << edge_ptr->edge_cov << ")\"];";
			o_graph << endl;



			if ((!edge_ptr->node_ptr->in_backbone))
			{
				if (Visited.count(edge_ptr->node_ptr) == 0)
				{
					Visited[edge_ptr->node_ptr]=1;
					node_list.push_back(edge_ptr->node_ptr);
				}
				if (edge_ptr->node_ptr->coord > 0)
				{
					o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [color=red,label=\"" << edge_ptr->node_ptr->kmer << "(" << edge_ptr->node_ptr->coord << ")" "\"];" << endl;

				}
				else
				{
					o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [color=red,label=\"" << edge_ptr->node_ptr->kmer << "\"];" << endl;

				}

			}
			else
			{

				o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [label=\"" << edge_ptr->node_ptr->kmer << "[" << edge_ptr->node_ptr->cns_coord << "]" << "[" << edge_ptr->node_ptr->coord << "]" "\"];" << endl;

			}




			edge_ptr = edge_ptr->nxt_edge;
		}
	}

	o_graph.close();

}



void OutputSubGraph(struct backbone_info *backbone_info, int begin, int end, string filename)
{


	ofstream o_graph(filename.c_str());
	o_graph << "digraph G {" << endl;
	o_graph << "\"" << backbone_info->node_vec[begin] << "\"" << " [label=\"" << backbone_info->node_vec[begin]->kmer << "(" << backbone_info->node_vec[begin]->coord << ")" "\"];" << endl;

	o_graph.close();
	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	map<consensus_node *, bool> Visited;
	for (int i = begin; i + 1 <end; ++i)
	{
		OutputPathsFromANode(backbone_info->node_vec[i], filename, Visited);

	}
	
	o_graph.open(filename.c_str(),ios_base::app);
	o_graph << "}" << endl;
	o_graph.close();
	o_graph.clear();
}





void OutputSparseSubGraph(struct backbone_info *backbone_info, int begin, int end, string filename)
{

	map<sparse_consensus_node *, bool> Visited;

	ofstream o_graph(filename.c_str());
	o_graph << "digraph G {" << endl;
	while (backbone_info->sparse_node_vec[begin] == NULL)
	{
		begin++;
	}
	o_graph << "\"" << backbone_info->sparse_node_vec[begin] << "\"" << " [label=\"" << backbone_info->sparse_node_vec[begin]->kmer << "(" << backbone_info->sparse_node_vec[begin]->coord << ")" "\"];" << endl;

	o_graph.close();
	uint32_t n_Rbranched = 0, n_Lbranched = 0;

	for (int i = begin; i + 1 <end; ++i)
	{
		if (backbone_info->sparse_node_vec[i] != NULL)
		{
			if (i == 82)
			{
				cout << "";
			}
			OutputPathsFromASparseNode(backbone_info->sparse_node_vec[i], filename,  Visited);
		}
		
	}

	o_graph.open(filename.c_str(), ios_base::app);
	o_graph << "}" << endl;
	o_graph.close();
	o_graph.clear();
}



char * multiply(const char *a_in, const char *b_in, char * mul)
{
	char c[500],a[500],b[500];
	strcpy(a, a_in);
	strcpy(b, b_in);
	char temp[500];
	int la, lb;
	int i, j, k = 0, x = 0, y;
	long int r = 0;
	long sum = 0;
	la = strlen(a) - 1;
	lb = strlen(b) - 1;


	for (i = 0; i <= la; i++){
		a[i] = a[i] - 48;
	}

	for (i = 0; i <= lb; i++){
		b[i] = b[i] - 48;
	}

	for (i = lb; i >= 0; i--){
		r = 0;
		for (j = la; j >= 0; j--){
			temp[k++] = (b[i] * a[j] + r) % 10;
			r = (b[i] * a[j] + r) / 10;
		}
		temp[k++] = r;
		x++;
		for (y = 0; y<x; y++){
			temp[k++] = 0;
		}
	}

	k = 0;
	r = 0;
	for (i = 0; i<la + lb + 2; i++){
		sum = 0;
		y = 0;
		for (j = 1; j <= lb + 1; j++){
			if (i <= la + j){
				sum = sum + temp[y + i];
			}
			y += j + la + 1;
		}
		c[k++] = (sum + r) % 10;
		r = (sum + r) / 10;
	}
	c[k] = r;
	j = 0;
	for (i = k - 1; i >= 0; i--){
		mul[j++] = c[i] + 48;
	}
	mul[j] = '\0';
	return mul;
}



void BFSFindBestPath(struct backbone_info *backbone_info,int node_idx)
{
	consensus_node * begin_node = backbone_info->node_vec[node_idx];
	list<consensus_node *> node_list;
	node_list.push_back(begin_node);
	int max_cov = backbone_info->cov_vec[node_idx];
	if (node_idx == 511)
	{
		cout << "";
	}
	while (node_list.size() > 0)
	{
		consensus_node * current_node = node_list.front();
		node_list.pop_front();
		consensus_edge_node * edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			//if (edge_ptr->node_ptr->score < (current_node->score + edge_ptr->edge_cov))
			bool update = 0;
			if (backbone_info->ScoringMethod == 1)
			{
				double score_org = 0.0;
				double score_new = 0.0;
				if (current_node->score == 0)
				{
					cout << "";
				}


				
				stringstream itoa_str1, itoa_str2, itoa_str3, itoa_str4;
				string score_str1, score_str2, degree_str1, degree_str2;

				itoa_str1 << edge_ptr->node_ptr->score;
				score_str1 = itoa_str1.str();
				itoa_str2 << edge_ptr->node_ptr->cov;
				degree_str1 = itoa_str2.str();

				itoa_str3 << current_node->score + edge_ptr->edge_cov;
				score_str2 = itoa_str3.str();
				itoa_str4 << current_node->cov + 1;
				degree_str2 = itoa_str4.str();
				char score_org_str[600], score_new_str[600];

				multiply(score_str1.c_str(), degree_str2.c_str(), score_org_str);
				multiply(score_str2.c_str(), degree_str1.c_str(), score_new_str);

				if (edge_ptr->node_ptr->score == 0)
				{
					update = 1;
				}
				else
				{
					string score_org2, score_new2;
					bool output = 0;
					for (int ii = 0; ii < strlen(score_org_str); ++ii)
					{
						if (score_org_str[ii] != '0')
						{
							output = 1;
						}
						if (output)
						{
							score_org2.push_back(score_org_str[ii]);
						}
					}
					output = 0;
					for (int ii = 0; ii < strlen(score_new_str); ++ii)
					{
						if (score_new_str[ii] != '0')
						{
							output = 1;
						}
						if (output)
						{
							score_new2.push_back(score_new_str[ii]);
						}
					}
					if (score_org2.size() < score_new2.size() || (score_org2.size() == score_new2.size() && score_org2<score_new2))
					{
						update = 1;
					}
				}
			}

			int new_score = 0;
			if (backbone_info->ScoringMethod == 2)
			{
				if (backbone_info->threshold < 0.0)
				{
					new_score = (current_node->score + max(edge_ptr->edge_cov - backbone_info->CovTh, -2));
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}
				else
				{
					int threshold = round((max_cov*backbone_info->threshold));
					if (threshold < backbone_info->CovTh)
					{
						threshold = backbone_info->CovTh;
					}
					new_score = (current_node->score + edge_ptr->edge_cov - threshold);
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}
				
			}
			
			
			

			//if (score_org <score_new)
			//
			if (update)
			{
				edge_ptr->node_ptr->cov = current_node->cov + 1;
				//edge_ptr->node_ptr->score = (current_node->score + max(edge_ptr->edge_cov - backbone_info->CovTh,-2));
					
				edge_ptr->node_ptr->score = new_score;

				edge_ptr->node_ptr->last_node = current_node;
				if (!edge_ptr->node_ptr->in_backbone)
				{
					node_list.push_back(edge_ptr->node_ptr);
				}

				
			}

			edge_ptr = edge_ptr->nxt_edge;
		}
	}

}


void BFSFindBestPathSparse(struct backbone_info *backbone_info, int node_idx)
{
	sparse_consensus_node * begin_node = backbone_info->sparse_node_vec[node_idx];
	list<sparse_consensus_node *> node_list;
	node_list.push_back(begin_node);
	int max_cov = backbone_info->cov_vec[node_idx];
	while (node_list.size() > 0)
	{
		sparse_consensus_node * current_node = node_list.front();
		node_list.pop_front();
		sparse_consensus_edge_node * edge_ptr = current_node->right;


		while (edge_ptr != NULL)
		{
			//if (edge_ptr->node_ptr->score < (current_node->score + edge_ptr->edge_cov))
			bool update = 0;

			int new_score = 0;
			if (backbone_info->ScoringMethod == 2)
			{
				


				if (backbone_info->threshold < 0.0)
				{
					new_score = current_node->score + edge_ptr->len*max((int)edge_ptr->edge_cov - backbone_info->CovTh, -2);
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}
				else
				{
					int threshold = round((max_cov*backbone_info->threshold));
					if (threshold < backbone_info->CovTh)
					{
						threshold = backbone_info->CovTh;
					}
					new_score = (current_node->score + edge_ptr->len*(edge_ptr->edge_cov - threshold));
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}

				

			}


			if (update)
			{
				edge_ptr->node_ptr->cov = current_node->cov + 1;
				edge_ptr->node_ptr->score = new_score;
			
				edge_ptr->node_ptr->last_node = current_node;
				if (!edge_ptr->node_ptr->in_backbone)
				{
					node_list.push_back(edge_ptr->node_ptr);
				}


			}

			edge_ptr = edge_ptr->nxt_edge;
		}
	}

}





void FindBestPath(struct backbone_info *backbone_info)
{
	
	uint32_t n_Rbranched = 0, n_Lbranched=0;
	backbone_info->node_vec[0]->cov = 1;//depth
	for (int i = 0; i + 1 <backbone_info->node_vec.size(); ++i)
	{
	
		BFSFindBestPath(backbone_info,i);

		//cout << i << " ";
		consensus_edge_node *edge_ptr = backbone_info->node_vec[i]->right;
		map<uint64_t, consensus_node *> node_map;
		map< consensus_node *, int> node_cov;
		int RBranch = 0;
		while (edge_ptr != NULL)
		{
			/*
			if (node_map.count(edge_ptr->node_ptr->kmer.kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer.kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] += edge_ptr->edge_cov;
			
			}
			*/
			RBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}
		if (RBranch > 1)
		{
			n_Rbranched++;
		}

		
		node_map.clear();
		node_cov.clear();

		edge_ptr = backbone_info->node_vec[i]->left;
		int LBranch = 0;
		while (edge_ptr != NULL)
		{
			/*
			if (node_map.count(edge_ptr->node_ptr->kmer.kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer.kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] += edge_ptr->edge_cov;

			}
			*/
			LBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}

		if (LBranch> 1)
		{
			n_Lbranched++;
		}

	}

	//cout << backbone_info->backbone.size() << " backbone size." << endl;
	//cout << backbone_info->n_nodes << " nodes." << endl;
	//cout << backbone_info->n_edges << " edges." << endl;
	//cout << n_Lbranched << " left branches." << endl;
	//cout << n_Rbranched << " right branches." << endl;

}




void FindBestPathSparse(struct backbone_info *backbone_info)
{

	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	backbone_info->sparse_node_vec[0]->cov = 1;//depth
	for (int i = 0; i + 1 <backbone_info->sparse_node_vec.size(); ++i)
	{
		if (backbone_info->sparse_node_vec[i] == NULL)
		{
			continue;
		}

		BFSFindBestPathSparse(backbone_info, i);

		//cout << i << " ";
		sparse_consensus_edge_node *edge_ptr = backbone_info->sparse_node_vec[i]->right;
		map<uint64_t, sparse_consensus_node *> node_map;
		map< sparse_consensus_node *, int> node_cov;
		int RBranch = 0;
		while (edge_ptr != NULL)
		{
			RBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}
		if (RBranch > 1)
		{
			n_Rbranched++;
		}


		node_map.clear();
		node_cov.clear();

		edge_ptr = backbone_info->sparse_node_vec[i]->left;
		int LBranch = 0;
		while (edge_ptr != NULL)
		{
			LBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}

		if (LBranch> 1)
		{
			n_Lbranched++;
		}

	}

	//cout << backbone_info->backbone.size() << " backbone size." << endl;
	//cout << backbone_info->n_nodes << " nodes." << endl;
	//cout << backbone_info->n_edges << " edges." << endl;
	//cout << n_Lbranched << " left branches." << endl;
	//cout << n_Rbranched << " right branches." << endl;

}




void BFSClear(struct backbone_info *backbone_info, int node_idx)
{
	consensus_node * begin_node = backbone_info->node_vec[node_idx];
	list<consensus_node *> node_list;
	node_list.push_back(begin_node);
	begin_node->in_backbone = 0;
	begin_node->coord = 0;
	while (node_list.size() > 0)
	{
		consensus_node * current_node = node_list.front();
		node_list.pop_front();
		consensus_edge_node * edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			edge_ptr->node_ptr->coord = 0;

			if (!edge_ptr->node_ptr->in_backbone)
			{
				node_list.push_back(edge_ptr->node_ptr);
			}
			

			edge_ptr = edge_ptr->nxt_edge;
		}
	}

}






void ClearInfo(struct backbone_info *backbone_info)
{

	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	backbone_info->node_vec[0]->cov = 1;//depth
	for (int i = 0; i + 1 <backbone_info->node_vec.size(); ++i)
	{

		BFSClear(backbone_info, i);

	}
	
	cout << "Nodes information clear." << endl;

}




#endif
