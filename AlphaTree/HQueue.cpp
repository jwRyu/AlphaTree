#include "defines.h"
#include "HQueue.hpp"
#include "allocator.h"

HQueue* hqueue_new(uint64 qsize, uint32 *dhist, uint32 dhistsize)
{
	uint32 i;
	HQueue* hqueue = (HQueue*)Malloc(sizeof(HQueue));
	hqueue->queue = (uint32*)Malloc((size_t)qsize * sizeof(uint32));
	hqueue->bottom = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));
	hqueue->cur = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));

	hqueue->qsize = qsize;
	hqueue->min_level = hqueue->max_level = dhistsize;

	int sum_hist = 0;
	for (i = 0; i < dhistsize; i++)
	{
		hqueue->bottom[i] = hqueue->cur[i] = sum_hist;
		sum_hist += dhist[i];
	}
	hqueue->bottom[dhistsize] = 0;
	hqueue->cur[dhistsize] = 1;

	return hqueue;
}

void hqueue_free(HQueue* hqueue)
{
	Free(hqueue->queue);
	Free(hqueue->bottom);
	Free(hqueue->cur);
	Free(hqueue);
}
