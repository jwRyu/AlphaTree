#pragma once

#include "defines.h"
#include "allocator.h"

template <typename T>
class HQueue
{
	T *queue;
	T *bottom, *cur;
public:
	uint64 qsize;
	T min_level, max_level;
	HQueue(uint64 qsize, T *dhist, uint32 dhistsize)
	{
		queue = (T*)Malloc((size_t)qsize * sizeof(T));
		bottom = (T*)Malloc((size_t)(dhistsize + 1) * sizeof(T));
		cur = (T*)Malloc((size_t)(dhistsize + 1) * sizeof(T));

		qsize = qsize;
		min_level = max_level = dhistsize;

		T sum_hist = 0;
		for (uint32 i = 0; i < dhistsize; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		bottom[dhistsize] = 0;
		cur[dhistsize] = 1;
	}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void hqueue_push(T pidx, T level)
	{
		min_level = min(level, min_level);
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
	}

	inline T hqueue_pop()
	{
		return queue[--cur[min_level]];
	}

	inline void hqueue_find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};
/*

struct neighbouridx {
	//  -  3  -
	//  2  p  1
	//  -  0  -
	uint8 neighbour;
	uint32 pidx;
};
template <>
class HQueue <neighbouridx>
{
	neighbouridx *queue;
	uint64 *bottom, *cur;
	uint64 qsize;
	uint64 min_level, max_level;
public:
	HQueue(uint64 qsize, uint64 *dhist, uint32 dhistsize, uint8 neighbours)
	{
		uint64 nn = neighbours >> 1;
		int shamt;
		queue = (neighbouridx*)Malloc((size_t)qsize * nn * sizeof(neighbouridx));
		bottom = (uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(uint64));
		cur = (uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(uint64));

		this->qsize = qsize;
		min_level = max_level = dhistsize;

		for (shamt = -1; nn; nn >>= 1)
			shamt++;

		uint64 sum_hist = 0;
		for (uint64 i = 0; i < dhistsize; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i] << shamt;
		}
		bottom[dhistsize] = 0;
		cur[dhistsize] = 1;
}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void hqueue_push(uint64 pidx, uint64 level)
	{
		min_level = min(level, min_level);
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
	}

	inline T hqueue_pop()
	{
		return queue[--cur[min_level]];
	}

	inline void hqueue_find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};
inline void hqueue_push(HQueue<neighidx>* hqueue, uint32 idx, uint8 neighbor, uint32 level)
{
	hqueue->min_level = min(level, hqueue->min_level);
#if DEBUG
	assert(level < hqueue->max_level);
	assert(hqueue->cur[level] < hqueue->qsize);
#endif
	hqueue->queue[hqueue->cur[level]].pidx = idx;
	hqueue->queue[hqueue->cur[level]++].neighbour = neighbor;
}
*/

/*
void hqueue_new(HQueue<uint32>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize);
void hqueue_new(HQueue<neighidx>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize, uint8 neighbours);
*/


