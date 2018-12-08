#pragma once

#include "defines.h"

typedef struct HQueue
{
	uint32 *queue, *bottom, *cur;
	uint64 qsize;
	uint32 min_level, max_level;
}HQueue;

HQueue* hqueue_new(uint64 qsize, uint32 *dhist, uint32 dhistsize);
void hqueue_free(HQueue* hqueue);
inline void hqueue_push(HQueue* hqueue, uint32 newidx, uint32 level)
{
	hqueue->min_level = min(level, hqueue->min_level);
#if DEBUG
	assert(level < hqueue->max_level);
	assert(hqueue->cur[level] < hqueue->qsize);
#endif
	hqueue->queue[hqueue->cur[level]++] = newidx;
}

inline uint32 hqueue_pop(HQueue* hqueue)
{
	return hqueue->queue[--hqueue->cur[hqueue->min_level]];
}

inline void hqueue_find_min_level(HQueue* hqueue)
{
	while (hqueue->bottom[hqueue->min_level] == hqueue->cur[hqueue->min_level])
		hqueue->min_level++;
}