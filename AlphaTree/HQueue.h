#pragma once

#include "defines.h"
#include "allocator.h"
#include "iostream"

template <class Imgidx>
class HQueue
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
public:
	int64 qsize;
	int64 min_level;
	HQueue(uint64 qsize_in, Imgidx *dhist, int32 numlevels)
	{
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));


		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void push(Imgidx pidx, int64 level)
	{
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
		if (level < min_level)
		{
			min_level = level;
			//return 1;
		}
//		else
			//return 0;
//		min_level = min(level, min_level);
	}

	inline Imgidx pop()
	{
		return queue[--cur[min_level]];
	}

	inline Imgidx top()
	{
		return queue[cur[min_level] - 1];
	}

	int64 get_minlev()
	{
		return min_level;
	}

	inline void find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};



template <class Imgidx>
class HQueue_hdr
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
	uint64 *seeker;
public:
	int64 qsize, seekersize;
	int64 min_level;
	HQueue_hdr(uint64 qsize_in, Imgidx *dhist, int32 numlevels)
	{
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		seekersize = (numlevels + 1 + 63) >> 6;
		seeker = (uint64 *)Malloc((size_t)(seekersize) * sizeof(uint64));

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		for (int64 i = 0; i < seekersize; i++)
			seeker[i] = 0;
		seeker[numlevels >> 6] |= (uint64)1 << (numlevels & 63);
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue_hdr()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
		Free(seeker);
	}

	inline int push(Imgidx pidx, int64 level)
	{
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		int64 qidx = cur[level]++;
		queue[qidx] = pidx;
		seeker[level >> 6] |= (uint64)1 << (level & 63);
		if (level <= min_level)
		{
			min_level = level;
			return 1;
		}
		return 0;
		//		min_level = min(level, min_level);
	}

	inline Imgidx pop()
	{
		Imgidx popidx = --cur[min_level];

		if(bottom[min_level] == cur[min_level])
			seeker[min_level >> 6] &= ~((uint64)1 << (min_level & 63));
		return queue[popidx];
	}

	inline Imgidx top()
	{
		return queue[cur[min_level] - 1];
	}

	inline int64 get_minlev()
	{
		return min_level;
	}

	inline void find_min_level()
	{
		Imgidx qidx, widx;
		uint64 w;

		for (qidx = min_level >> 6; !seeker[qidx]; qidx++)
			;

		w = seeker[qidx];
		
		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(uint64)1))
		{
			w >>= 1;
			widx++;
		}

		min_level = ((qidx << 6) + widx);
// 
// 		while (bottom[min_level] == cur[min_level])
// 			min_level++;
// 
// 
// 		if ( != min_level)
// 			qidx = qidx;
	}
};


template <class Imgidx>
class HQueue_hdr2
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
	uint64 *seeker,*seeker2;
public:
	int64 qsize;
	int64 min_level;
	HQueue_hdr2(uint64 qsize_in, Imgidx *dhist, int32 numlevels)
	{
		int64 seekersize, seeker2size;
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		
		seekersize = (numlevels + 1 + 63) >> 6;
		seeker2size = (seekersize + 63) >> 6;
		
		seeker = (uint64 *)Malloc((size_t)(seekersize) * sizeof(uint64));
		seeker2 = (uint64 *)Malloc((size_t)(seeker2size) * sizeof(uint64));

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		for (int64 i = 0; i < seekersize; i++)
			seeker[i] = 0;
		seeker[numlevels >> 6] |= (uint64)1 << (numlevels & 63);
		for (int64 i = 0; i < seeker2size; i++)
			seeker2[i] = 0;
		seeker2[numlevels >> 12] |= (uint64)1 << ((numlevels >> 6) & 63);
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue_hdr2()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
		Free(seeker);
		Free(seeker2);
	}

	inline void push(Imgidx pidx, int64 level)
	{
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif	
		int64 qidx = cur[level]++;
		queue[qidx] = pidx;
		seeker[level >> 6] |= (uint64)1 << (level & 63);
		seeker2[level >> 12] |= (uint64)1 << ((level >> 6) & 63);
		if (level < min_level)
		{
			min_level = level;
			//return 1;
		}
		//		else
					//return 0;
		//		min_level = min(level, min_level);
	}

	inline Imgidx pop()
	{
		Imgidx popidx = --cur[min_level];

		if (bottom[min_level] == cur[min_level])
		{
			seeker[min_level >> 6] &= ~((uint64)1 << (min_level & 63));
			if (!seeker[min_level >> 6])
				seeker2[min_level >> 12] &= ~((uint64)1 << ((min_level >> 6) & 63));
		}
		return queue[popidx];
	}

	inline Imgidx top()
	{
		return queue[cur[min_level] - 1];
	}

	int64 get_minlev()
	{
		return min_level;
	}

	inline void find_min_level()
	{
		Imgidx qidx, widx;
		uint64 w;

		for (qidx = min_level >> 12; !seeker2[qidx]; qidx++)
			;

		w = seeker2[qidx];
		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(uint64)1))
		{
			w >>= 1;
			widx++;
		}

		qidx = ((qidx << 6) + widx);

		w = seeker[qidx];
		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(uint64)1))
		{
			w >>= 1;
			widx++;
		}

		min_level = ((qidx << 6) + widx);
		// 
		// 		while (bottom[min_level] == cur[min_level])
		// 			min_level++;
		// 
		// 
		// 		if ( != min_level)
		// 			qidx = qidx;
	}
};



/*
template <class Imgidx>
class HQueue_ubr
{
	int8 *queue;
//	Imgidx *bottom, *cur;
public:
	uint64 qsize;
	Imgidx min_level;
	HQueue_ubr(uint64 qsize)
	{
		//tmp
		queue = (int8*)Malloc((size_t)(qsize + 1) * sizeof(int8));
		min_level = (Imgidx)(qsize);

		for (uint64 i = 0; i < qsize; i++)
			queue[i] = 0;
		//queue[qsize] = 1;
//		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx), 1);
	//	cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx), 1);


		// 		qsize = qsize_in;
		// 		min_level = numlevels - 1;
		// 
		// 		Imgidx sum_hist = 0;
		// 		for (int32 i = 0; i < numlevels; i++)
		// 		{
		// 			bottom[i] = cur[i] = sum_hist;
		// 			sum_hist += dhist[i];
		// 		}
		// 		bottom[numlevels] = 0;
		// 		cur[numlevels] = 1;
	}
	~HQueue_ubr()
	{
		Free(queue);
		//Free(bottom, 1);
		//Free(cur, 1);
	}

	inline void push(Imgidx pidx)
	{
		min_level = min(pidx, min_level);
		queue[pidx] = 1;
	}

	inline Imgidx top()
	{
		return min_level;
	}

	inline void pop()
	{
		queue[min_level] = 0;
		//std::cout << "popping from" << min_level << std::endl;
		while (!queue[min_level])
			min_level++;
	}

	inline void find_min_level()
	{
		
	}
};
*/

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
