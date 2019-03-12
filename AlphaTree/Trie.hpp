#pragma once
#include "allocator.h"

template<class Imgidx, class Trieidx>
class Trie
{
	Imgidx minidx;
	Trieidx **trie;
	//Imgidx *levelsize;
	Imgidx imgsize, mask_field, mask_msb;
	int8 shamt, numlevels;
	//delayed non-leaf node push

public:
	Trie(Imgidx imgsize)
	{
		Imgidx size, lvlsz;
		shamt = 2;
		for (int8 nbyte = sizeof(Trieidx); nbyte; nbyte >>= 1)
			shamt++;
		mask_field = (1 << shamt) - 1;
		mask_msb = 1 << (shamt - 1);
		this->imgsize = imgsize;
		numlevels = 1;
		for (size = imgsize >> shamt; size; size >>= shamt)
			numlevels++;
		
		trie = (Trieidx**)Malloc(sizeof(Trieidx*) * numlevels);
		//levelsize = (Imgidx*)Malloc(sizeof(Imgidx) * numlevels);
		size = imgsize;
		for (int8 i = 0; i < numlevels; i++)
		{
			size >>= shamt;
			trie[i] = (Trieidx*)Malloc(sizeof(Trieidx) * (size + 1));
			for (Imgidx j = 0; j < (size + 1); j++)
				trie[i][j] = 0;
			//levelsize[i] = lvlsz;
		}
		minidx = imgsize;
		this->push(imgsize);
	}
	~Trie()
	{
		for (int8 i = 0; i < numlevels; i++)
			Free(trie[i]);
		Free(trie);
	}

	void push(Imgidx in)
	{
		Imgidx s_in = in >> shamt, shamt1;
		Trieidx *p;

		if (in < minidx)
			minidx = in;
		p = &(trie[0][s_in]);
		*p = *p | ((Trieidx)1 << (in & mask_field));
		in = s_in;
		s_in >>= shamt;
		for (int8 i = 1; i < numlevels; i++)
		{
			p = &(trie[i][s_in]);
			shamt1 = in & mask_field;
			if (((*p) >> shamt1) & 1)
				break;
			*p = *p | ((Trieidx)1 << shamt1);
			in = s_in;
			s_in >>= shamt;
		}
	}
	void pop()
	{
		Imgidx s_idx = minidx >> shamt, shamt1;
		Trieidx *p, tmp;
		int8 lvl;

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

