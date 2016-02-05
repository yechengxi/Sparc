#include "iostream"
#include "stdio.h"
#include "string"
#include "vector"
#include "cstdlib"
#include "bitset"
#include <map>
#include <math.h>
#include "memory"
#include <algorithm>
#include "fstream"
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
#include "stdint.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
using namespace std;
/*

take a blasr output and calculate consensus
19 fields

qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
string qName,tName,qAlignedSeq,matchPattern,tAlignedSeq;
int qLength,qStart,qEnd,tLength,tStart,tEnd,score,numMatch,numMismatch,numIns,numDel,mapQV;
char qStrand,tStrand;
m130404_014004_sidney_c100506902550000001823076808221337_s1_p0/1059/0_303/0_303 303 12 296 +  Backbone_2 77265 53349 53598 - -1025 246 0 38 3 254  ...

*/


int main(int argc, char* argv[])
{

	bool HELP = 0;
	string Filename, BackboneFile, OutputFilename="Consensus";
	int K_size = 2;
	int CovTh = 2;
	int ScoringMethod = 2;
	int subgraph_begin = 0, subgraph_end = 0;
	int cns_begin = 0, cns_end = 0;
	int report_begin = 0, report_end = 0;
	int boost = 1;
	bool Debug = 0;
	string ContigPrefix;
	int gap = 2;
	bool Patch = 0, Fill = 0;
	int Patch_K = 5;
	int Patch_D = 30;
	int Patch_G = 2;
	double threshold=-0.1;
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "h") == 0)
		{
			i++;
			HELP = 1;
			break;
		}
		if (strcmp(argv[i], "f") == 0 || strcmp(argv[i], "m") == 0)
		{
			i++;
			Filename = (argv[i]);
			continue;
		}
		if (strcmp(argv[i], "o") == 0 || strcmp(argv[i], "out") == 0)
		{
			i++;
			OutputFilename = (argv[i]);
			continue;
		}
		if (strcmp(argv[i], "b") == 0)
		{
			i++;
			BackboneFile = (argv[i]);
			continue;
		}
		if (strcmp(argv[i], "k") == 0)
		{
			i++;
			K_size = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "t") == 0)
		{
			i++;
			threshold = atof(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Patch") == 0)
		{
			i++;
			Patch = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Fill") == 0)
		{
			i++;
			Fill = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Patch_K") == 0)
		{
			i++;
			Patch_K = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Patch_D") == 0)
		{
			i++;
			Patch_D = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Patch_G") == 0)
		{
			i++;
			Patch_G = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Debug") == 0 || strcmp(argv[i], "DEBUG") == 0)
		{
			i++;
			Debug = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "g") == 0)
		{
			i++;
			gap = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "c") == 0)
		{
			i++;
			CovTh = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "sb") == 0)
		{
			i++;
			subgraph_begin = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "se") == 0)
		{
			i++;
			subgraph_end = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "rb") == 0)
		{
			i++;
			report_begin = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "re") == 0)
		{
			i++;
			report_end = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "cb") == 0)
		{
			i++;
			cns_begin = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "ce") == 0)
		{
			i++;
			cns_end = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "Score") == 0)
		{
			i++;
			ScoringMethod = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "HQ_Prefix") == 0 || strcmp(argv[i], "Contig_Prefix") == 0 || strcmp(argv[i], "HQPrefix") == 0 || strcmp(argv[i], "ContigPrefix") == 0)
		{
			i++;
			ContigPrefix = argv[i];
			continue;
		}
		if (strcmp(argv[i], "boost") == 0)
		{
			i++;
			boost = atoi(argv[i]);
			continue;
		}

	}
	cout << "For help: Sparc -h" << endl;
	if (HELP || argc<=2)
	{
		cout << " Example command: " << endl;

		cout << "Sparc b BackboneFile.fa m Mapped.m5 c 2 k 2 g 2 o ConsensusOutput" << endl;
		
		cout << "Parameters:" << endl;
		cout << "b: backbone file." << endl;
		cout << "m: the reads mapping files produced by blasr, using option -m 5." << endl;
		cout << "k: k-mer size. (range: [1,5])" << endl;
		cout << "c: coverage threshold. (range: [1,5])" << endl;
		//cout << "t: adaptive threshold. (range: [0.0,0.3]), default 0.0." << endl;
		cout << "g: skip size, the larger the value, the more memory efficient the algorithm is. (range: [1,5])" << endl;
		cout << "HQ_Prefix: Shared prefix of the high quality read names." << endl;
		cout << "boost: boosting weight for the high quality reads. (range: [1,5])" << endl;
		cout << "Author: Chengxi Ye cxy@umd.edu." << endl;
		cout << "last update: Jan 2, 2015." << endl;
		return 1;
	}


	ifstream backbone_in(BackboneFile.c_str());
	ifstream mapped_in(Filename.c_str());
	string str,backbone,tag,n_tag;
	get_a_fasta_read(backbone_in,tag,backbone,n_tag);

	align_profile backbone_align_profile;
	backbone_align_profile.match_vec.resize(backbone.size());
	backbone_align_profile.mismatch_vec.resize(backbone.size());
	backbone_align_profile.insertion_vec.resize(backbone.size());
	backbone_align_profile.deletion_vec.resize(backbone.size());

	
	ref_read_t ref;
	ref.read_bits = (uint64_t *)malloc((size_t)(backbone.size() / 4) + 100);
	ref.alloc_sz = (size_t)(backbone.size() / 4 + 100);
	Init_Ref_Read(backbone, ref);


	int64_t bucket_count=0, edge_cnt=0;
	struct backbone_info backbone_info;
	backbone_info.threshold = threshold;
	backbone_info.gap = gap;
	backbone_info.ContigPrefix = ContigPrefix;
	backbone_info.ScoringMethod = ScoringMethod;
	backbone_info.boost = boost;
	backbone_info.backbone = backbone;
	backbone_info.n_edges = 0;
	backbone_info.n_nodes = 0;
	backbone_info.CovTh = CovTh;
	//cout << "constructing backbone graph." << endl;
	struct backbone_info backbone_info_org = backbone_info;

	if (gap == 1)
	{
		Consensus_Kmer_Graph_Construction(&ref, &backbone_info_org, K_size);
	}
	
	if (Debug == 1 || gap != 1)
	{
		
		Consensus_Sparse_Kmer_Graph_Construction(&ref, &backbone_info, K_size);
	}
	
	//cout << "Nodes: " << backbone_info.n_nodes << " nodes." << endl;

	//cout << "Edges: " << backbone_info.n_edges << " edges." << endl;
	cout << "Backbone size: " << backbone_info.backbone.size() << endl;

	ofstream o_profile;
	if (Debug)
	{
		o_profile.open("align_profile.txt");

	}
	
	//cout << "adding query branches." << endl;

	size_t n_reads = 0;
	backbone_info.ref_matched = 0;
	backbone_info.ref_mismatch = 0;
	backbone_info_org.ref_matched = 0;
	backbone_info_org.ref_mismatch = 0;
	struct query_info query_info, query_info_org;
	query_info.Patch = Patch;
	query_info.Fill = Fill;
	query_info.Patch_K = Patch_K;
	query_info.Patch_D = Patch_D;
	query_info.Patch_G = Patch_G;
	query_info.report_b = report_begin;
	query_info.report_e = report_end;
	
	if (report_end > report_begin)
	{
		ofstream o_report_align("subalign.txt");
	}

	while (getline(mapped_in, str))
	{
		stringstream ss(str);
		
		ss >> query_info.qName >> query_info.qLength >> query_info.qStart >> query_info.qEnd >> query_info.qStrand >> query_info.tName >> query_info.tLength >> query_info.tStart >> query_info.tEnd >> query_info.tStrand >> query_info.score >> query_info.numMatch >> query_info.numMismatch >> query_info.numIns >> query_info.numDel >> query_info.mapQV >> query_info.qAlignedSeq >> query_info.matchPattern >> query_info.tAlignedSeq;
		
		query_info.n_exist = 0;
		query_info.n_new = 0;
		
		n_reads++;
		query_info.read_idx = n_reads;
		query_info_org = query_info;
		if ( gap == 1)
		{
			Add_Path_To_Backbone(&backbone_info_org, &query_info_org, K_size);
		}

		if (Debug == 1 || gap != 1)
		{
			
			Add_Path_To_Backbone_Sparse(&backbone_info, &query_info, K_size);
			
		}

	
		
	}

	//cout << "done." << endl;;
	if (gap == 1)
	{
		if (Debug)
		{
			cout << "Nodes: " << backbone_info_org.n_nodes << "." << endl;
			cout << "Edges: " << backbone_info_org.n_edges << "." << endl;

		}
		uint64_t cum_sum = 0;
		
		map<int, int> cov_cnt;
		
		backbone_info_org.cov_vec.resize(backbone_info_org.node_vec.size());
		int radius = 200;
		if (radius > backbone.size())
		{
			radius = backbone.size() - 1;
		}
		for (int i = 0; i < radius; ++i)
		{
			cov_cnt[backbone_info_org.node_vec[i]->cov]++;
		}

		for (int i = 0; i<backbone_info_org.node_vec.size(); ++i)
		{
			if (i >= radius)
			{
				cov_cnt[backbone_info_org.node_vec[i - radius]->cov]--;
				if (cov_cnt[backbone_info_org.node_vec[i - radius]->cov] <= 0)
				{
					cov_cnt.erase(backbone_info_org.node_vec[i - radius]->cov);
				}
			}
			if (i + radius < backbone_info_org.node_vec.size())
			{
				cov_cnt[backbone_info_org.node_vec[i+radius]->cov]++;
			}
			
			
			if (cov_cnt.size() > 0)
			{
				backbone_info_org.cov_vec[i] = cov_cnt.rbegin()->first;
			}
			else
			{
				backbone_info_org.cov_vec[i] = 0;
			}
		}

		
		FindBestPath(&backbone_info_org);
		string filename = "subgraph.dot";
		if (subgraph_end > subgraph_begin)
		{
			OutputSubGraph(&backbone_info_org, subgraph_begin, subgraph_end, filename);

		}

	}
	
	if (Debug == 1 || gap != 1)
	{
		if (Debug)
		{
			cout << "Nodes: " << backbone_info.n_nodes << "." << endl;

			cout << "Edges: " << backbone_info.n_edges << "." << endl;
		}


		map<int, int> cov_cnt;

		backbone_info.cov_vec.resize(backbone_info.sparse_node_vec.size());
		int radius = 200;
		for (int i = 0; i < radius; ++i)
		{
			if (backbone_info.sparse_node_vec[i] != NULL)
			{
				cov_cnt[backbone_info.sparse_node_vec[i]->cov]++;
			}
			
		}

		for (int i = 0; i<backbone_info.sparse_node_vec.size(); ++i)
		{
			if (i >= radius)
			{
				if (backbone_info.sparse_node_vec[i - radius] != NULL)
				{
					cov_cnt[backbone_info.sparse_node_vec[i - radius]->cov]--;
					if (cov_cnt[backbone_info.sparse_node_vec[i - radius]->cov] <= 0)
					{
						cov_cnt.erase(backbone_info.sparse_node_vec[i - radius]->cov);
					}
				}
				
			}
			if ((i + radius < backbone_info.sparse_node_vec.size()) && (backbone_info.sparse_node_vec[i + radius] != NULL))
			{
				cov_cnt[backbone_info.sparse_node_vec[i + radius]->cov]++;
			}


			if (cov_cnt.size() > 0)
			{
				backbone_info.cov_vec[i] = cov_cnt.rbegin()->first;
			}
			else
			{
				backbone_info.cov_vec[i] = 0;
			}
		}



		FindBestPathSparse(&backbone_info);
		string filename = "sparse_subgraph.dot";
		if (subgraph_end > subgraph_begin)
		{
			OutputSparseSubGraph(&backbone_info, subgraph_begin, subgraph_end, filename);

		}
	}
	
	

	string consensus;
	if (gap == 1)
	{

		uint64_t max_score = 0;
		size_t position = 0;

		for (int i = 0; i < backbone_info_org.node_vec.size(); ++i)
		{
			if (backbone_info_org.node_vec[i]->score>max_score)
			{
				max_score = backbone_info_org.node_vec[i]->score;
				position = i;
			}
			backbone_info_org.node_vec[i]->in_backbone = 0;
		}

		if (max_score == 0)
		{
			string OutputFilename2 = OutputFilename + ".consensus.fasta";
			ofstream o_cns(OutputFilename2.c_str());
			o_cns << tag << endl;
			o_cns << backbone << endl;
			cout << "Empty ouput. Backbone copied." << endl;
			return -1;
		}

		consensus_node * current_node = backbone_info_org.node_vec[position];

		char KmerStr[100];
		if (current_node != NULL)
		{
			uint64_t kmer_uint64 = (current_node->kmer);
			bitsarr2str(&kmer_uint64, K_size, KmerStr, 1);
			consensus = KmerStr;
			reverse(consensus.begin(), consensus.end());
			current_node = current_node->last_node;
			while (current_node != NULL)
			{
				kmer_uint64 = (current_node->kmer);
				bitsarr2str(&kmer_uint64, K_size, KmerStr, 1);

				consensus.push_back(KmerStr[0]);
				current_node = current_node->last_node;
			}

			reverse(consensus.begin(), consensus.end());

		}

		if (Debug)
		{
			string OutputFilename2 = "DEBUG_" + OutputFilename + ".consensus.fasta";
			ofstream o_cns(OutputFilename2.c_str());
			o_cns << ">" << OutputFilename << endl;
			o_cns << consensus << endl;

		}
		

	
		
		int cns_pos = consensus.size();
		current_node = backbone_info_org.node_vec[position];
		backbone_info_org.node_vec.clear();
		backbone_info_org.node_vec.push_back(current_node);
		if (current_node != NULL)
		{
			current_node = current_node->last_node;
			while (current_node != NULL)
			{
				backbone_info_org.node_vec.push_back(current_node);
				cns_pos--;
				current_node = current_node->last_node;
			}
		}
		reverse(backbone_info_org.node_vec.begin(), backbone_info_org.node_vec.end());

		for (int ii = 0; ii<backbone_info_org.node_vec.size(); ++ii)
		{
			backbone_info_org.node_vec[ii]->in_backbone = 1;
			backbone_info_org.node_vec[ii]->cns_coord = ii + 1;
		}
		string filename = "subgraph_cns.dot";
		if (cns_end >cns_begin)
		{
			OutputSubGraph(&backbone_info_org, cns_begin, cns_end, filename);

		}
		if (Debug)
		{
			for (int i = 0; i < backbone_info_org.node_vec.size(); ++i)
			{
				if (i % 10 == 0)
				{
					o_profile << i << " " << backbone_info_org.node_vec[i]->score;

					if (backbone_info_org.node_vec[i]->coord>0)
					{
						o_profile << " bb_coord: " << backbone_info_org.node_vec[i]->coord;
					}
					o_profile << " cns_coord: " << backbone_info_org.node_vec[i]->cns_coord << endl;

				}

			}
		}
		
		
	}

	if (Debug == 1 || gap != 1)
	{
		
		uint64_t max_score = 0;
		size_t position = 0;
		for (int i = 0; i < backbone_info.sparse_node_vec.size(); ++i)
		{
			if (backbone_info.sparse_node_vec[i]!=NULL&&backbone_info.sparse_node_vec[i]->score>max_score)
			{
				max_score = backbone_info.sparse_node_vec[i]->score;
				position = i;
				
			}
			if (backbone_info.sparse_node_vec[i] != NULL)
			{
				backbone_info.sparse_node_vec[i]->in_backbone = 0;
			}
		}

		sparse_consensus_node * current_node = backbone_info.sparse_node_vec[position];
		backbone_info.sparse_node_vec.clear();
		backbone_info.sparse_node_vec.push_back(current_node);
		if (max_score == 0)
		{
			string OutputFilename2 = OutputFilename + ".consensus.fasta";
			ofstream o_cns(OutputFilename2.c_str());
			o_cns << tag << endl;
			o_cns << backbone << endl;
			cout << "Empty ouput." << endl;
			return -1;
		}
		char KmerStr[100];
		if (current_node != NULL)
		{
			uint64_t kmer_uint64 = (current_node->kmer);
			bitsarr2str(&kmer_uint64, K_size, KmerStr, 1);
			consensus = KmerStr;
			reverse(consensus.begin(), consensus.end());
			sparse_consensus_edge_node * edge_node = current_node->left;
			sparse_consensus_edge_node * best_edge_node = NULL;
			int max_edge_cov = 0;
			while (edge_node != NULL)
			{
				if (edge_node->node_ptr == current_node->last_node)
				{
					if (edge_node->edge_cov >= max_edge_cov)
					{
						max_edge_cov = edge_node->edge_cov;
						best_edge_node = edge_node;
						//break;
					}
					
				}
				edge_node = edge_node->nxt_edge;
			}
			uint64_t edge_bits = best_edge_node->edge;
			
			char edge_str[100];
			bitsarr2str(&edge_bits, best_edge_node->len, edge_str, 1);
			for (int i = best_edge_node->len - 1; i >= 0; --i)
			{
				consensus.push_back(edge_str[i]);
				if (i > 0)
				{
					backbone_info.sparse_node_vec.push_back(NULL);
				}
			}
			
			current_node = current_node->last_node;

			while (current_node->last_node != NULL)
			{	
				backbone_info.sparse_node_vec.push_back(current_node);
				edge_node = current_node->left;
				best_edge_node = NULL;
				max_edge_cov = 0;
				while (edge_node != NULL)
				{
					if (edge_node->node_ptr == current_node->last_node)
					{
						if (edge_node->edge_cov > max_edge_cov||max_edge_cov==0)
						{
											
							max_edge_cov = edge_node->edge_cov;
							best_edge_node = edge_node;
							//break;
						}
					}
					edge_node = edge_node->nxt_edge;
				}

				uint64_t edge_bits = best_edge_node->edge;
				
				char edge_str[100];
				bitsarr2str(&edge_bits, best_edge_node->len, edge_str, 1);

				uint64_t kmer_bits = current_node->last_node->kmer;

				char kmer_str[100];
				bitsarr2str(&kmer_bits, K_size, kmer_str, 1);
				if (best_edge_node->len == K_size&&strcmp(kmer_str, edge_str) != 0)
				{
					cout << "Error." << endl;
					return -1;
					cout << kmer_str << endl;
					cout << edge_str << endl;
				}
				for (int i = best_edge_node->len - 1; i >= 0; --i)
				{
					consensus.push_back(edge_str[i]);
					if (i > 0)
					{
						backbone_info.sparse_node_vec.push_back(NULL);
					}
				}

				current_node = current_node->last_node;

			}

			reverse(consensus.begin(), consensus.end());
			reverse(backbone_info.sparse_node_vec.begin(), backbone_info.sparse_node_vec.end());

			for (int ii = 0; ii<backbone_info.sparse_node_vec.size(); ++ii)
			{
				if (backbone_info.sparse_node_vec[ii] != NULL)
				{
					backbone_info.sparse_node_vec[ii]->in_backbone = 1;
					backbone_info.sparse_node_vec[ii]->cns_coord = ii + 1;
				}
			}
			string filename = "sparse_subgraph_cns.dot";
			if (cns_end >cns_begin)
			{
				OutputSparseSubGraph(&backbone_info, cns_begin, cns_end, filename);

			}

			if (Debug)
			{
				for (int i = 0; i < backbone_info.sparse_node_vec.size(); ++i)
				{
					if (backbone_info.sparse_node_vec[i] != NULL)
					{
						if (i % 10 == 0)
						{
							o_profile << i << " " << backbone_info.sparse_node_vec[i]->score;
							if (backbone_info.sparse_node_vec[i]->coord>0)
							{
								o_profile << " bb_coord: " << backbone_info.sparse_node_vec[i]->coord;
							}
							o_profile << " cns_coord: " << backbone_info.sparse_node_vec[i]->cns_coord << endl;

						}

					}


				}
			}
			


		}

		
	}

	string OutputFilename2 = OutputFilename + ".consensus.fasta";
	ofstream o_cns(OutputFilename2.c_str());

	int t = 0;
	for (int i = 0; i < OutputFilename.size(); ++i)
	{
		if (OutputFilename[i] == '\\' ||OutputFilename[i]=='/')
		{
			t = i+1;
		}
	}

	o_cns << tag << endl;
	if (consensus.size()>0)
	{
		o_cns << consensus << endl;

	}
	else
	{
		o_cns << backbone << endl;

	}
	cout << "Finished." << endl;
	
	//ClearInfo(&backbone_info);
	
	/*
	o_profile << "#match, mismatch, deletion" << endl;
	for (int i = 0; i < backbone_align_profile.match_vec.size(); ++i)
	{
		o_profile << backbone_align_profile.match_vec[i] << ", ";
		o_profile << backbone_align_profile.mismatch_vec[i] << ", ";
		o_profile << backbone_align_profile.deletion_vec[i] << ", ";
		o_profile << backbone_align_profile.insertion_vec[i] << ", ";
		o_profile << endl;
	}
	*/

	return 0;
}