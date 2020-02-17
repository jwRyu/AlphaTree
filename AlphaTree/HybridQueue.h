#pragma once
#include "allocator.h"
#include "defines.h"
#include "Trie.h"
#include "HierarQueue.h"

#include <iostream>
using namespace std;

#define LISTSIZE_DEFAULT 64
#define HEAPSIZE_DEFAULT 128

#define TRACK_QUEUEING	0

//tmp
#if TRACK_QUEUEING
#include <fstream>
using namespace std;
#endif


template<class Imgidx, class Trieidx>//, class Qidx>
class HybridQueue_Trie
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	Trie<Imgidx, Trieidx> *trie;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
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
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		trie = new Trie<Imgidx, int64>(size);
		list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
		list[0] = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_Trie(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_Trie(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		int16 i;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (idx < trie->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list]);

				for (i = curSize_list - 1; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		trie->push(idx);
	}
	inline void pop()
	{
		int8 i;
		Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0] = trie->top();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		trie->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
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

	~HybridQueue_Trie()
	{
		delete trie;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};



template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue_Rank
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	HQueue_l1idx_rank<Imgidx> *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
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
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		queue = new HQueue_l1idx_rank<Imgidx>(size);
		list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
		list[0] = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue_Rank(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue_Rank(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

// 	double get_jumpnum() { return queue->jumpnum; }
// 	double get_jumpdist() { return queue->jumpdist; }

	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		int16 i;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (idx < queue->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list]);

				for (i = curSize_list - 1; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		queue->push(idx);
	}
	inline void pop()
	{
		int8 i;
		Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0] = queue->top();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
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

	~HybridQueue_HQueue_Rank()
	{
		delete queue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};



template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue_Rank1
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	HQueue_l1idx_rank<Imgidx> *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;

	int16 l0, mask;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
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
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		queue = new HQueue_l1idx_rank<Imgidx>(size);

		if (listsize > 256)
			listsize = 256;
		else if (listsize > 128)
			listsize = 256;
		else if (listsize > 64)
			listsize = 128;
		else if (listsize > 32)
			listsize = 64;
		else if (listsize > 16)
			listsize = 32;
		else if (listsize > 8)
			listsize = 16;
		else if (listsize > 4)
			listsize = 8;
		else
			listsize = 4;

		mask = listsize - 1;
		list = (Imgidx*)Malloc((listsize) * sizeof(Imgidx));
		maxSize_list = listsize;
		curSize_list = 0;
		l0 = listsize >> 1;
		list[l0] = size;
		
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue_Rank1(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue_Rank1(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	// 	double get_jumpnum() { return queue->jumpnum; }
	// 	double get_jumpdist() { return queue->jumpdist; }

	inline Imgidx get_minlev() { return list[l0]; }
	inline Imgidx top() { return list[l0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		int16 i, j, lm;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		lm = (l0 - 1) & mask;
		if (idx < queue->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				if (idx < list[l0])
				{
					list[lm] = idx;
					l0 = lm;
				}
				else
				{
					i = (l0 + curSize_list) & mask;
					j = (l0 + curSize_list - 1) & mask;
					while (idx < list[j])
					{
						list[i] = list[j];
						i = j;
						j = (j - 1) & mask;
					}
					list[i] = idx;
				}
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[lm]);

				if (idx < list[l0])
				{
					list[lm] = idx;
					l0 = lm;
				}
				else
				{
					i = lm;
					j = (lm - 1) & mask;
					while (idx < list[j])
					{
						list[i] = list[j];
						i = j;
						j = (j - 1) & mask;
					}
					list[i] = idx;
				}
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		queue->push(idx);
	}
	inline void pop()
	{
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 1)
		{
			list[l0] = queue->top();

			pop_queue();
		}
		else
		{
			l0 = (l0 + 1) & mask;
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
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

	~HybridQueue_HQueue_Rank1()
	{
		delete queue;
		Free(list);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};


// template<class Imgidx>//, class Qidx>
// class HybridQueue_HH //Heap + Hierarchical
// {
// 	//MinList1<Imgidx> *list, *list_end, *head, *tail;
// 	HeapQueue_rank *heapqueue;				//prior queue
// 	HQueue_l1idx_rank<Imgidx> *hierarqueue; //secondary queue
// 	Imgidx minidx_hierarqueue, maxidx_heapqueue;
// //	int16 curSize_list, maxSize_list; heapqueue->cursize
// 	Imgidx maxSize_queue, mask_field;
// 	int8 shamt, nbit;
// 
// 
// #if TRACK_QUEUEING
// 	Imgidx *in_size;
// 
// 	ofstream f;
// #endif
// 	//	int32 cnt;
// 	void initHQ(Imgidx size, size_t heapsize)
// 	{
// 		/*		cnt = 0;//tmp*/
// 		Imgidx i;
// 		this->maxSize_queue = size;
// 		/*		shamt = 2;*/
// 		// 		nbit = sizeof(Qidx) * 8;
// 		// 		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
// 		// 			shamt++;
// 		// 		mask_field = (1 << shamt) - 1;
// 		// 		qsize = (size + mask_field) >> shamt;
// 		// 
// 		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));	
// 		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
// 		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
// 		hierarqueue = new HQueue_l1idx_rank<Imgidx>(size);
// 		heapqueue = new HeapQueue_rank<Imgidx>(heapsize);
// 		maxidx_heapqueue = 0;
// // 		list = (Imgidx*)Malloc((heapsize) * sizeof(Imgidx));
// //		list[0] = 0;
// //		list++;
// //		maxSize_list = heapsize - 1;
// //		curSize_list = -1;
// 		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
// 		//list_end = list + listsize;
// 		//maxSize_list = listsize;
// 		//head = tail = 0;
// 
// 
// 		//		for (i = 0; i < size; i++)
// 		//			queue[i] = -1;
// 		//		queue[size] = 0;
// 		//for (i = 0; i < listsize; i++)
// 			//list[i].idx = -1;
// 		//curSize_list = 0;
// 		//		minidx_queue = size >> shamt;
// 		minidx_hierarqueue = size;
// 
// 
// 		//tmp
// #if TRACK_QUEUEING
// 		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
// 		f << -1 << '\n' << size << endl;
// #endif
// 	}
// public:
// 	HybridQueue_HQueue_Rank(Imgidx size)
// 	{
// 		initHQ(size, HEAPSIZE_DEFAULT);
// 	}
// 	HybridQueue_HQueue_Rank(Imgidx size, size_t heapsize)
// 	{
// 		initHQ(size, heapsize);
// 	}
// 
// 	// 	double get_jumpnum() { return queue->jumpnum; }
// 	// 	double get_jumpdist() { return queue->jumpdist; }
// 
// 	inline Imgidx get_minlev() { return heapqueue->get_minlev(); }
// 	inline Imgidx top() { return heapqueue->top(); }
// 	inline void push(Imgidx idx)
// 	{
// 		//MinList1<Imgidx> *p, *q;
// 
// #if TRACK_QUEUEING
// 		//tmp
// 		f << '0' << '\n' << idx << endl;
// #endif
// 
// 		// 		cnt++;//tmp
// 		// 
// 		// 		if (cnt == 786)//tmp
// 		// 			idx = idx;
// 		if (idx < hierarqueue->top())
// 		{
// 			if (heapqueue->cursize < heapqueue->maxsize) //spare room in the list
// 			{
// 				maxidx_heapqueue = maxidx_heapqueue > idx ? maxidx_heapqueue: idx;
// 				heapqueue->push(idx);
// 			}
// 			else if (idx < )// push to the full list
// 			{
// 				for (Imgidx i = heapqueue->cursize; heapqueue->get_arr(i) != maxidx_heapqueue; i--)
// 					;
// 				
// 			}
// 			else
// 				push_queue(idx); // push to the queue
// 		}
// 		else
// 			push_queue(idx); // push to the queue
// 	}
// 	inline void push_queue(Imgidx idx)
// 	{
// 		hierarqueue->push(idx);
// 	}
// 	inline void pop()
// 	{
// 		int8 i;
// 		Imgidx idx;
// 		// 		cnt++;//tmp
// 		// 		if (cnt == 776)//tmp
// 		// 			cnt = cnt;
// 
// 
// 		//tmp
// 
// #if TRACK_QUEUEING
// 		f << '1' << '\n' << list[0] << endl;
// #endif
// 
// 		if (curSize_list == 0)
// 		{
// 			list[0] = hierarqueue->top();
// 
// 			pop_queue();
// 		}
// 		else
// 		{
// 			for (i = 0; i < curSize_list; i++)
// 				list[i] = list[i + 1];
// 			curSize_list--;
// 		}
// 
// #if TRACK_QUEUEING
// 		f << list[0] << endl;
// #endif
// 	}
// 	inline void pop_queue()
// 	{
// 		hierarqueue->pop();
// 		// 		queue[minidx_queue] = -1;
// 		// 		while (queue[++minidx_queue] == -1)
// 		// 			;
// 	}
// 
// 	// 	int8 checklist()//tmp
// 	// 	{
// 	// 		MinList<Imgidx> *p;
// 	// 		if (head)
// 	// 		{
// 	// 			for (p = head; p; p = p->next)
// 	// 			{
// 	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
// 	// 					return 1;
// 	// 			}
// 	// 			if (tail->next)
// 	// 				return 1;
// 	// 		}
// 	// 		int n = 0;
// 	// 		for (int i = 0; i < maxSize_list; i++)
// 	// 		{
// 	// 			if (list[i].idx != -1)
// 	// 				n++;
// 	// 		}
// 	// 		if (n != curSize_list)
// 	// 			return 1;
// 	// 		return 0;
// 	// 	}
// 
// 	~HybridQueue_HQueue_Rank()
// 	{
// 		delete hierarqueue;
// 		Free(list - 1);
// 		//Free(queue);
// 
// #if TRACK_QUEUEING
// 		//tmp
// 		f.close();
// #endif
// 	}
// };



template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	int64 *levels;
	HQueue_l1idx<Imgidx> *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 minlevnotfixed;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	int32 cnt;
	void initHQ(uint64 qsize_in, Imgidx *dhist, int32 numlevels, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		// 
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));	
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		queue = new HQueue_l1idx<Imgidx>(qsize_in, dhist, numlevels);
		list = (Imgidx*)Malloc((listsize) * sizeof(Imgidx));
		levels = (int64*)Malloc((listsize + 1) * sizeof(int64));
		levels[0] = 0;
		levels++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = qsize_in;

		minlevnotfixed = 0;
		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue(uint64 qsize_in, Imgidx *dhist, int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue(uint64 qsize_in, Imgidx *dhist, int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}
	inline Imgidx top(){return list[0];}
	inline int64 get_minlev() { return levels[0]; }
	inline void find_minlev()
	{
//		if (minlevnotfixed)
			queue->find_minlev();

	}
	inline void push(Imgidx idx, int64 level)
	{
		//MinList1<Imgidx> *p, *q;
		int16 i;

//		cout << "- pushing " << (int)idx << " at level " << (int)level << endl;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;

		if (curSize_list == -1) //should be run only the first time
		{
			list[0] = idx;
			levels[0] = level;
			curSize_list = 0;
			push_queue(idx, level);
			return;
		}

		if (level < queue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; level < levels[i]; i--)
				{
					list[i + 1] = list[i];
					levels[i + 1] = levels[i];
				}
				list[i + 1] = idx;
				levels[i + 1] = level;
				curSize_list++;
			}
			else if (level < levels[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list], levels[curSize_list]);

				for (i = curSize_list - 1; level < levels[i]; i--)
				{
					list[i + 1] = list[i];
					levels[i + 1] = levels[i];
				}
				list[i + 1] = idx;
				levels[i + 1] = level;
			}
			else
				push_queue(idx, level); // push to the queue
		}
		else
			push_queue(idx, level); // push to the queue
	}
	inline void push_queue(Imgidx idx, int64 level)
	{
		if (queue->push(idx, level))
			minlevnotfixed = 0;
	}
	inline Imgidx pop()
	{
		Imgidx ret;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;

//		cout << "- popping " << (int)list[0] << " at level " << (int)levels[0] << endl;

		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		ret = list[0];

		if (curSize_list == 0)
		{
			list[0] = queue->top();
			levels[0] = queue->get_minlev();

			pop_queue();
		}
		else
		{
			for (int8 i = 0; i < curSize_list; i++)
			{
				list[i] = list[i + 1];
				levels[i] = levels[i + 1];
			}
			curSize_list--;
		}
		return ret;

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		minlevnotfixed = 1;
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
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

	~HybridQueue_HQueue()
	{
		delete queue;
		Free(list);
		Free(levels - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

