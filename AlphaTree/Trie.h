#pragma once
#include "allocator.h"
// 
// //tmptmptmptmpt tmp
// #include <fstream>
// using namespace std;


template<class Imgidx, class Trieidx>
class Trie
{
	Imgidx minidx;
	Trieidx **trie;
	//Imgidx *levelsize;
	Imgidx triesize, mask_field, mask_msb;
	int8 shamt, numlevels;
	//delayed non-leaf node push
		
// 	//tmptmptmp
// 	ofstream f;

public:
	Trie(Imgidx triesize)
	{
		Imgidx size;
		shamt = 2;
		for (int8 nbyte = sizeof(Trieidx); nbyte; nbyte >>= 1)
			shamt++;
		mask_field = (1 << shamt) - 1;
		mask_msb = 1 << (shamt - 1);
		numlevels = 1;
		for (size = triesize >> shamt; size; size >>= shamt)
			numlevels++;
		
		trie = (Trieidx**)Malloc(sizeof(Trieidx*) * numlevels);
		//levelsize = (Imgidx*)Malloc(sizeof(Imgidx) * numlevels);
		size = triesize;
		for (int8 i = 0; i < numlevels; i++)
		{
			size >>= shamt;
			trie[i] = (Trieidx*)Malloc(sizeof(Trieidx) * (size + 1));
			for (Imgidx j = 0; j < (size + 1); j++)
				trie[i][j] = 0;
			//levelsize[i] = lvlsz;
		}
		minidx = triesize;

		//tmp!
// 		f.open("C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/trie0rrr.dat",std::ofstream::out);
// 		f << triesize << endl;
	}
	~Trie()
	{
		for (int8 i = 0; i < numlevels; i++)
			Free(trie[i]);
		Free(trie);
		
// 		f.close();//tmptmp
	}

	inline Imgidx top() { return minidx; }
	inline Imgidx min_rank() { return minidx >> 1; }
	inline void push(Imgidx in, int8 incidence)
	{
		Imgidx n, s_in, shamt1;
		Trieidx *p;

		//tmp
/*		f << '0' << '\n' << in << endl;*/

		n = (in << 1) + incidence;
		s_in = n >> shamt;

		if (n < minidx)
			minidx = n;
		p = &(trie[0][s_in]);
		*p = *p | ((Trieidx)1 << (n & mask_field));
		n = s_in;
		s_in >>= shamt;
		for (int8 i = 1; i < numlevels; i++)
		{
			p = &(trie[i][s_in]);
			shamt1 = n & mask_field;
 			//if (((*p) >> shamt1) & 1)
 				//break;
			*p = *p | ((Trieidx)1 << shamt1);
			n = s_in;
			s_in >>= shamt;
		}
	}
	inline void pop()
	{
		Imgidx s_idx = minidx >> shamt, shamt1;
		Trieidx *p, tmp;
		int8 lvl;

// 		//tmp
// 		f << '1' << '\n' << (minidx>>1) << endl;

		shamt1 = minidx & mask_field;
		p = &(trie[0][s_idx]);
		*p = *p & (~((Trieidx)1 << shamt1++));
		for (lvl = 0; !*p;)
		{
			minidx = s_idx;
			s_idx >>= shamt;
			shamt1 = minidx & mask_field;
			p = &(trie[++lvl][s_idx]);
			*p = *p & (~((Trieidx)1 << shamt1));
		}
		for (tmp = *p >> shamt1; !(tmp & 1); tmp >>= 1)
			shamt1++;
		minidx = (minidx & ~mask_field) | shamt1;
		while (lvl)
		{
			tmp = trie[--lvl][minidx];
			p = &tmp;
			s_idx = minidx;
			minidx <<= shamt;
			for (shamt1 = 0; !(tmp & 1); shamt1++)
				tmp >>= 1;
			minidx |= shamt1;
		}
	}
};
