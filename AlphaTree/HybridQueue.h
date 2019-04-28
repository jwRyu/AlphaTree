#pragma once
#include "allocator.h"
#include "defines.h"

template<class Imgidx>
struct MinList
{
	Imgidx idx;
	Imgidx *next;
};

template<class Imgidx, class Qidx>
class HybridQueue
{
	MinList *list, *head;
	Qidx *queue;

	Imgidx size, qsize;
	int8 shamt;
	void initHQ(Imgidx size, size_t listsize)
	{
		this->size = size;
		shamt = 2;
		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
			shamt++;
		Imgidx roundup = (1 << shamt) - 1;
		qsize = (size + roundup) >> shamt;
		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		
		list = (MinList*)Malloc(listsize * sizeof(MinList));
		head = NULL;

	}
public:
	HybridQueue(Imgidx size)
	{
		this->size = size;
		shamt = 2;
		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
			shamt++;
		Imgidx roundup = (1 << shamt) - 1;
		qsize = (size + roundup) >> shamt;
		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));

		//list = 
	}
	HybridQueue(Imgidx size, size_t listsize)
	{
		
	}
	~HybridQueue()
	{
		Free(queue);
	}
};

