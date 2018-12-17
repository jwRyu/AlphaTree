#pragma once

#include "defines.h"
#include "allocator.h"

template <int NN>
class NeighbourList;


// 4-neighbour List
template <>
class NeighbourList<4>
{
	uint8 *count; //2bit counter
	uint8 *list;  //list of closest neighbours. 8bit for each pixel. 2bit for each neighbour. The closest is at the LSB position
	uint32 imgsize;
	/*
	inline uint8 get_count(uint32& idx) const
	{
		return count[idx >> 2] >> ((idx & 3) << 1);
	}
	*/
	inline uint8 countdown(uint32& idx)
	{
		uint32 idx1 = idx >> 1, idx2 = (idx & 3) << 1;
		uint8 cnt = count[idx1] >> idx2;
		count[idx1] = (count[idx1] & (~(3 << idx2))) | (((--cnt) & 3) << idx2);
		return cnt == 0;
	}

	inline uint8 countup(uint32& idx)
	{
		uint32 idx1 = idx >> 2, idx2 = (idx & 3) << 1;
		uint8 cnt = count[idx1] >> idx2;
		uint8 ret = cnt & 3;
		count[idx1] = (count[idx1] & (~(3 << idx2))) | (((cnt + 1) & 3) << idx2);
		return ret;
	}

public:
	NeighbourList(uint32 img_size): imgsize(img_size)
	{
		list = (uint8*)Calloc(1, img_size);
		count = (uint8*)Calloc(1, (img_size +3)>>2);
	}
	~NeighbourList()
	{
		Free(list);
	}

	inline void insert(uint32& idx, uint8& neighbour)
	{
		uint8 entry = list[idx], cnt = countup(idx);
		list[idx] = entry | (neighbour << (cnt<<1));
	}
};

