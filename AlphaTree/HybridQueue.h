#pragma once
#include "allocator.h"
#include "defines.h"

//#include <iostream>

#define LISTSIZE_DEFAULT 32

template<class Imgidx>
struct MinList
{
	Imgidx idx;
	int8 incidence;
	MinList *next;
};

template<class Imgidx>//, class Qidx>
class HybridQueue
{
	MinList<Imgidx> *list, *list_end, *head, *tail;
	int8 *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;	
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;

//	int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
/*		cnt = 0;//tmp*/
		Imgidx i;
		this->maxSize_queue = size;
/*		shamt = 2;*/
// 		nbit = sizeof(Qidx) * 8;
// 		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
// 			shamt++;
// 		mask_field = (1 << shamt) - 1;
// 		qsize = (size + mask_field) >> shamt;
// 
// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));	
		queue = (int8*)Malloc((size + 1) * sizeof(int8));
		list = (MinList<Imgidx>*)Malloc(listsize * sizeof(MinList<Imgidx>));
		list_end = list + listsize;
		maxSize_list = listsize;
		head = tail = 0;

		for (i = 0; i < size; i++)
			queue[i] = -1;
		queue[size] = 0;
		for (i = 0; i < listsize; i++)
			list[i].idx = -1;
		curSize_list = 0;
//		minidx_queue = size >> shamt;
		minidx_queue = size;
	}
public:
	HybridQueue(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}
	inline Imgidx top_idx() { return head->idx; }
	inline Imgidx top_incidence() { return head->incidence; }
	inline void push(Imgidx idx, int8 incidence)
	{
		MinList<Imgidx> *p, *q;
		int16 i;

// 		cnt++;//tmp
// 
// 		if (cnt == 786)//tmp
// 			idx = idx;
		if (curSize_list < maxSize_list) //spare room in the list
		{
			curSize_list++;
			
			for (q = &list[0]; q < list_end; q++)
			{
				if (q->idx == -1)
				{
					q->idx = idx;
					q->incidence = incidence;
					break;
				}
			}
			if (head == 0)
			{
				head = tail = q;
				q->next = 0;
			}
			else if (idx < head->idx)
			{
				q->next = head;
				head = q;
			}
			else
			{
				for (p = head; p->next && p->next->idx < idx; p = p->next)
					;
				q->next = p->next;
				p->next = q;
				if (tail->idx < q->idx)
					tail = q;
			}
// 			if (checklist())//tmp
// 				head = head;
		}
		else if (idx < tail->idx) // push to the full list
		{
			push_queue(tail->idx, tail->incidence);
			
			if (idx < head->idx)
			{
				tail->next = head;
				head = tail;
				for (p = list; p->next != tail; p++)
					;
				p->next = 0;
				tail->idx = idx;
				tail->incidence = incidence;
				tail = p;
			}
			else
			{
				//std::cout << cnt++ ;
				for (p = head; p->next && p->next->idx < idx; p = p->next)
					;


				tail->next = p->next;
				p->next = tail;
				for (q = tail; q->next != tail; q = q->next)
					;

				//std::cout << "end" << std::endl;//yogi debug
				q->next = 0;
				tail->idx = idx;
				tail->incidence = incidence;
				tail = q;
			}

// 			if (checklist())//tmp
// 				head = head;
		}
		else
			push_queue(idx, incidence); // push to the queue
	}
	inline void push_queue(Imgidx idx, int8 incidence)
	{
// 		int8 bitfield_idx;
// 		idx = (idx << 1) | (Imgidx)incidence;
// 		bitfield_idx = idx & mask_field;
// 		idx = idx >> shamt;
// 		if (idx < minidx_queue)
// 			minidx_queue = idx;
// 		queue[idx] = queue[idx] | (1 << bitfield_idx);
		if (idx < minidx_queue)
			minidx_queue = idx;
		queue[idx] = incidence;
	}
	void pop()
	{
		int8 qi;
		Imgidx idx;
		MinList<Imgidx> *p;
// 		cnt++;//tmp
// 		if (cnt == 776)//tmp
// 			cnt = cnt;
		if (head == tail)
		{
			head->idx = -1;
			head = tail = &list[0];
			head->idx = minidx_queue;
			head->incidence = queue[minidx_queue];
			head->next = 0;
			tail = head;

			queue[minidx_queue] = -1;
			while (queue[++minidx_queue] == -1)
				;
// 
// 			while (1)
// 			{
// 				p = &list[curSize_list];
// 				p->idx = minidx_queue;
// 				p->incidence = queue[minidx_queue];
// 				
// 				queue[minidx_queue] = -1;
// 				while (queue[++minidx_queue] == -1)
// 					;
// 				if (minidx_queue >= maxSize_queue)
// 					minidx_queue = minidx_queue;
// 				//algorithm might crash without early termination
// 				//if (curSize_list == maxSize_list - 1 || minidx_queue == maxSize_queue) 
// 				if (curSize_list == maxSize_list - 1)
// 				{
// 					p->next = 0;
// 					break;
// 				}
// 				p->next = p + 1;
// 				curSize_list++;
// 			}
// 			while (curSize_list < maxSize_list)
// 			{
// 				qi = queue[minidx_queue];
// 				idx = minidx_queue << shamt;
// 				for (i = 0; qi && (curSize_list < maxSize_list); i++)
// 				{
// 					if (qi & 1)
// 					{
// 						list[curSize_list].idx = (idx | i) >> 1;
// 						list[curSize_list].incidence = i & 1;
// 						list[curSize_list].next = &list[curSize_list + 1];
// 						curSize_list++;
// 					}
// 					qi >>= 1;
// 				}
// 				
// 				for (queue[minidx_queue] = qi << i; queue[minidx_queue++] == 0;)
// 					;
// 				//end?
// 			}
// 			if (checklist())//tmp
// 				head = head;
		}
		else
		{
			curSize_list--;
// 			if (head->next->idx < 0)//tmp
// 				head = head;
			head->idx = -1;
			head = head->next;
// 			if (checklist())//tmp
// 				head = head;
		}
	}

// 	int8 checklist()//tmp
// 	{
// 		MinList<Imgidx> *p;
// 		if (head)
// 		{
// 			for (p = head; p; p = p->next)
// 			{
// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
// 					return 1;
// 			}
// 			if (tail->next)
// 				return 1;
// 		}
// 		int n = 0;
// 		for (int i = 0; i < maxSize_list; i++)
// 		{
// 			if (list[i].idx != -1)
// 				n++;
// 		}
// 		if (n != curSize_list)
// 			return 1;
// 		return 0;
// 	}

	~HybridQueue()
	{
		Free(list);
		Free(queue);
	}
};

