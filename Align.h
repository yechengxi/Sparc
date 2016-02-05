
#ifndef __ALIGN_H
#define __ALIGN_H


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



struct aln_t
{
	const static int MaxLen = 301;
	int score[MaxLen][MaxLen];
	uint8_t b[MaxLen][MaxLen];
	int insertions, deletions, substitutions;
};
void str2int(char *s, int * out)
{
	int str_sz = strlen(s);
	for (int i = 0; i<str_sz; ++i)
	{
		char temp = s[i];
		switch (temp)
		{
		case 'A':
			out[i] = 1;
			break;
		case 'a':
			out[i] = 1;
			break;
		case 'C':
			out[i] = 2;
			break;
		case 'c':
			out[i] = 2;
			break;
		case 'G':
			out[i] = 3;
			break;
		case 'g':
			out[i] = 3;
			break;
		case 'T':
			out[i] = 4;
			break;
		case 't':
			out[i] = 4;
			break;
		default:
			cout << "error" << endl;
			return;
		}
	}

}

int GlobalAlign(struct aln_t *aln_t, char *A, char *B, int match, int mismatch, int gap_cost, int band_width)
{
	int A_sz = strlen(A), B_sz = strlen(B);
	int Ai[500], Bi[500];
	memset(aln_t->score, 0, sizeof(int)*(aln_t->MaxLen*aln_t->MaxLen));
	memset(aln_t->b, 0, sizeof(uint8_t)*(aln_t->MaxLen*aln_t->MaxLen));

	str2int(A, Ai);
	str2int(B, Bi);
	for (int i = 1; i <= A_sz; ++i)
	{

		int idA = Ai[i - 1];
		int upper_bound, lower_bound;
		if (B_sz>A_sz)
		{
			upper_bound = min(B_sz - A_sz + i + band_width, B_sz);
			lower_bound = max(1, i - band_width);
		}
		else
		{
			lower_bound = i - band_width - (A_sz - B_sz);
			if (lower_bound < 1)
			{
				lower_bound = 1;
			}
			upper_bound = i + band_width;
			if (upper_bound > B_sz)
			{
				upper_bound = B_sz;
			}
		}
		for (int j = lower_bound; j <= upper_bound; ++j)
		{

			int idB = Bi[j - 1];
			int Sij = -1000000, Sij_max = -1000000, from = 0;
			if (j<=upper_bound)
			{

				Sij = aln_t->score[i - 1][j] + gap_cost;
				Sij_max = Sij;
				from = 1;

			}

			if (j>=lower_bound)
			{
				Sij = aln_t->score[i][j - 1] + gap_cost;
				if (Sij>Sij_max)
				{
					Sij_max = Sij;
					from = 2;
				}
			}
			if (idA == idB)
			{
				Sij = aln_t->score[i - 1][j - 1] + match;
			}
			else
			{
				Sij = aln_t->score[i - 1][j - 1] + mismatch;
			}
			if (Sij>Sij_max)
			{
				Sij_max = Sij;
				from = 3;
			}

			aln_t->score[i][j] = Sij_max;
			aln_t->b[i][j] = from;

		}


	}

	/*
	ofstream o_sc("score.txt");
	ofstream o_path("path.txt");

	for(int i=1;i<=A_sz;++i)
	{
	for(int j=1;j<=B_sz;++j)
	{
	o_sc<<aln_t->score[i][j]<<" ";
	}
	o_sc<<endl;
	}
	for(int i=1;i<=A_sz;++i)
	{
	for(int j=1;j<=B_sz;++j)
	{
	int path=aln_t->b[i][j];
	o_path<<(path)<<" ";
	}
	o_path<<endl;
	}
	*/
	
	int max_score = -1000000;
	for (int i = B_sz; i >= 0; --i)
	{
		if ((aln_t->score[A_sz][i])>max_score)
		{
			max_score = aln_t->score[A_sz][i];
		}
	}
	return max_score;//aln_t->score[A_sz][B_sz];



}

void printAlign(struct aln_t *aln_t, char *A, char *B, string &A_aln, string &B_aln)
{
	A_aln.clear();
	B_aln.clear();
	int i, j, max_score = -100000000;
	int A_sz = strlen(A), B_sz = strlen(B);
	aln_t->deletions = 0;
	aln_t->insertions = 0;
	aln_t->substitutions = 0;
	A_aln.resize(2 * max(A_sz, B_sz) + 10);
	B_aln.resize(2 * max(A_sz, B_sz) + 10);

	/*
	for(int k=1;k<=B_sz;++k)
	{
	if(aln_t->score[k][A_sz]>max_score)
	{i=k;}
	}
	*/
	i = A_sz;
	j = B_sz;

	int position = 0;


	max_score = aln_t->score[A_sz][B_sz];

	while (1)
	{
		//	cout<<i<<" "<<j<<endl;//36//43

		if (i == 0)
		{

			for (int k = j; k >= 1; --k)
			{
				
				A_aln[position] = '-';
				B_aln[position] = B[k - 1];
				position++;
			}

			break;
		}


		if (j == 0)
		{
			//the head of B is missing here
			for (int k = i; k >= 1; --k)
			{
				aln_t->insertions++;//A has more bases 
				A_aln[position] = A[k - 1];
				B_aln[position] = '-';
				position++;
			}
			break;
		}
		if (aln_t->b[i][j] == 3 || aln_t->b[i][j] == 0)
		{

			A_aln[position] = A[i - 1];
			B_aln[position] = B[j - 1];
			position++;
			if (A[i - 1] != B[j - 1])
			{
				aln_t->substitutions++;
			}
			i--; j--;
			continue;
		}
		if (aln_t->b[i][j] == 2)
		{
			
			aln_t->deletions++;//A is missing a base			
			A_aln[position] = '-';
			B_aln[position] = B[j - 1];
			position++;
			
			j--;
			continue;
		}
		if (aln_t->b[i][j] == 1)
		{
			aln_t->insertions++;//A has one more base
			A_aln[position] = A[i - 1];
			B_aln[position] = '-';
			position++;

			i--;
			continue;
		}
	}
	A_aln.resize(position);
	B_aln.resize(position);

	reverse(A_aln.begin(), A_aln.end());
	reverse(B_aln.begin(), B_aln.end());
}


#endif
