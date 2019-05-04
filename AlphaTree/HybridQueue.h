#pragma once
#include "allocator.h"
#include "defines.h"
#include "Trie.h"

//#include <iostream>

#define LISTSIZE_DEFAULT 16

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
	Trie<Imgidx, int64> *trie;
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
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		trie = new Trie<Imgidx, int64>(size << 1);
		list = (MinList<Imgidx>*)Malloc(listsize * sizeof(MinList<Imgidx>));
		list_end = list + listsize;
		maxSize_list = listsize;
		head = tail = 0;

//		for (i = 0; i < size; i++)
//			queue[i] = -1;
//		queue[size] = 0;
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
//		if (idx < minidx_queue)
//			minidx_queue = idx;
//		queue[idx] = incidence;
		trie->push(idx, incidence);
	}
	inline void pop()
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
			//head->idx = minidx_queue;
			//head->incidence = queue[minidx_queue];
			head->idx = trie->min_rank();
			head->incidence = trie->min_incidence();
			head->next = 0;
			tail = head;
			
			pop_queue();
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

	~HybridQueue()
	{
		delete trie;
		Free(list);
		//Free(queue);
	}
};

template<class Imgidx>//, class Qidx>
class HybridQueue_1
{
	MinList<Imgidx> *list, *head, *tail, *first, *last;
	Trie<Imgidx, int64> *trie;
	int8 *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;

	//int32 cnt;//tmp
	void initHQ(Imgidx size, size_t listsize)
	{
		//cnt = 0;//tmp
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
		trie = new Trie<Imgidx, int64>(size << 1);
		list = (MinList<Imgidx>*)Malloc((listsize + 2) * sizeof(MinList<Imgidx>));
		tail = &list[listsize];
		head = &list[listsize + 1];
		head->next = tail;
		tail->idx = size;
		first = last = tail;
		maxSize_list = listsize;

		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		for (i = 0; i < listsize; i++)
			list[i].idx = -1;
		curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;
	}
public:
	HybridQueue_1(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_1(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}
	inline Imgidx top_idx() { return first->idx; }
	inline Imgidx top_incidence() { return first->incidence; }
	inline void push(Imgidx idx, int8 incidence)
	{
		MinList<Imgidx> *p, *q;
		int16 i;

		//cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (curSize_list < maxSize_list) //spare room in the list
		{
			curSize_list++;

			for (q = list; q < tail; q++)
			{
				if (q->idx == -1)
				{
					q->idx = idx;
					q->incidence = incidence;
					break;
				}
			}

			for (p = head; p->next->idx < idx; p = p->next)
				;
			q->next = p->next;
			p->next = q;
			if (p == head)
				first = q;
			if (q->next == tail)
				last = q;
		}
		else if (idx < last->idx) // push to the full list
		{
			push_queue(last->idx, last->incidence);
			last->idx = idx;
			last->incidence = incidence;
			for (p = list; p->next != last; p++)
				;
			p->next = tail;

			if (idx < first->idx)
			{
				last->next = first;
				head->next = last;
				first = last;
				last = p;
			}
			else
			{
				//std::cout << cnt++ ;
				for (q = first; q->next->idx < idx; q = q->next)
					;

				last->next = q->next;
				q->next = last;
				if (p->next == tail)
					last = p;
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
//		if (idx < minidx_queue)
//			minidx_queue = idx;
//		queue[idx] = incidence;
		trie->push(idx, incidence);
	}
	inline void pop()
	{
		int8 qi;
		Imgidx idx;
		MinList<Imgidx> *p;
		//cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
		if (curSize_list == 1)
		{
			//first->idx = -1;
			//head = tail = &list[0];
			//head->idx = minidx_queue;
			//head->incidence = queue[minidx_queue];
			first->idx = trie->min_rank();
			first->incidence = trie->min_incidence();
			//head->next = 0;
			//tail = head;

			pop_queue();
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
			first->idx = -1;
			first = first->next;
			head->next = first;
			// 			if (checklist())//tmp
			// 				head = head;
		}
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
	// 		if (curSize_list == 1)
	// 		{
	// 			if (head->next != first)
	// 				return 1;
	// 			if (head->next != last)
	// 				return 1;
	// 			if (head->next->next != tail)
	// 				return 1;
	// 		}
	// 		else
	// 		{
	// 			int i;
	// 			p = first;
	// 			for (i = 0; i < curSize_list; i++)
	// 			{
	// 				if (i == curSize_list - 1 && p != last)
	// 					return 1;
	// 				if (p == tail)
	// 					return 1;
	// 				if (p->idx > p->next->idx)
	// 					return 1;
	// 				p = p->next;
	// 			}
	// 		}
	// 		return 0;
	// 	}

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

	~HybridQueue_1()
	{
		delete trie;
		Free(list);
		//Free(queue);
	}
};



//////////////////////////////////////////////////////////////////////////
//HQ NEW?
// #pragma once
// #include "allocator.h"
// #include "defines.h"
// 
// //#include <iostream>
// 
// //#define LISTSIZE_DEFAULT 32
// 
// // template<class Imgidx>
// // struct MinList
// // {
// // 	Imgidx idx;
// // 	int8 incidence;
// // 	MinList *next;
// // };
// 
// template<class Imgidx>
// struct MinList
// {
// 	MinList<Imgidx>* next;
// 	int8 incidence;
// };
// 
// template<class Imgidx>//, class Qidx>
// class HybridQueue
// {
// 	//MinList<Imgidx> *list, *list_end, *head, *tail;
// 	MinList<Imgidx> *list, *head;
// 	int8 *queue;
// 	uint8 *qlocalhist;
// 	Imgidx qhsize, qsize, hsize, listsize;
// 	Imgidx minidx_queue;
// 	//int16 curSize_list;
// //	Imgidx maxSize_queue, mask_field;
// 	int8 shamt_qhist;
// 
// 	//int32 cnt;
// 	void initHQ(Imgidx size)
// 	{
// 		//	cnt = 0;//tmp
// 		Imgidx i;
// 		//this->maxSize_queue = size;
// /*		shamt = 2;*/
// // 		nbit = sizeof(Qidx) * 8;
// // 		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
// // 			shamt++;
// // 		mask_field = (1 << shamt) - 1;
// // 		qsize = (size + mask_field) >> shamt;
// // 
// // 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));	
// 		//qsize = (((size + 1) >> 1) + 63) >> 6;
// 		shamt_qhist = 8;
// 		listsize = (size + 1) >> 1;
// 		listsize = (size + 1);//tmp
// 		//qsize = size - listsize + 1;
// 		this->qhsize = size;
// 		qsize = size + 1;// -listsize + 1;
// 		hsize = (size + 1 + 255) >> shamt_qhist;
// 
// 		queue = (int8*)Malloc(qsize * sizeof(int8));
// 		qlocalhist = (uint8*)Malloc(hsize * sizeof(uint8));
// 		for (i = 0; i < qsize; i++)
// 			queue[i] = -1;
// 
// 		//	queue -= listsize;
// 
// 		queue[size] = 0;
// 		for (i = 0; i < hsize; i++)
// 			qlocalhist[i] = 0;
// 		qlocalhist[hsize - 1] = 1;
// 		list = (MinList<Imgidx>*)Malloc(listsize * sizeof(MinList<Imgidx>));//(MinList<Imgidx>*)Malloc(listsize * sizeof(MinList<Imgidx>));
// 		head = 0;
// 		// 		list_end = list + listsize;
// 		// 		maxSize_list = listsize;
// 		// 		head = tail = 0;
// 
// 				//for (i = 0; i < listsize; i++)
// 					//list[i].idx = -1;
// 				//curSize_list = 0;
// 		//		minidx_queue = size >> shamt;
// 		minidx_queue = size;
// 	}
// public:
// 	HybridQueue(Imgidx size)
// 	{
// 		initHQ(size);
// 	}
// 	// 	HybridQueue(Imgidx size, size_t listsize)
// 	// 	{
// 	// 		initHQ(size);
// 	// 	}
// 	inline void top(Imgidx& idx, int8& incidence)
// 	{
// 		if (head)
// 		{
// 			idx = (Imgidx)(head - list);
// 			incidence = head->incidence;
// 		}
// 		else
// 		{
// 			idx = minidx_queue;
// 			incidence = queue[minidx_queue];
// 		}
// 	}
// 	// 	inline Imgidx top_idx() { return head->idx; }
// 	// 	inline Imgidx top_incidence() { return head->incidence; }
// 	inline void push(Imgidx idx, int8 incidence)
// 	{
// 		MinList<Imgidx> *p, *q;
// 		/*		cnt++;*/
// 		if (idx < (qhsize >> 1))
// 		{
// 			q = &list[idx];
// 			q->incidence = incidence;
// 
// 			if (head)
// 			{
// 				if (q < head)
// 				{
// 					q->next = head;
// 					head = q;
// 				}
// 				else
// 				{
// 					for (p = head; p->next && p->next < q; p = p->next)
// 						;
// 					q->next = p->next;
// 					p->next = q;
// 				}
// 			}
// 			else
// 			{
// 				head = q;
// 				q->next = 0;
// 			}
// 			// 			if (checklist())//tmp
// 			// 				head = head;
// 		}
// 		// 		else if (idx < tail->idx) // push to the full list
// 		// 		{
// 		// 			push_queue(tail->idx, tail->incidence);
// 		// 			
// 		// 			if (idx < head->idx)
// 		// 			{
// 		// 				tail->next = head;
// 		// 				head = tail;
// 		// 				for (p = list; p->next != tail; p++)
// 		// 					;
// 		// 				p->next = 0;
// 		// 				tail->idx = idx;
// 		// 				tail->incidence = incidence;
// 		// 				tail = p;
// 		// 			}
// 		// 			else
// 		// 			{
// 		// 				//std::cout << cnt++ ;
// 		// 				for (p = head; p->next && p->next->idx < idx; p = p->next)
// 		// 					;
// 		// 
// 		// 
// 		// 				tail->next = p->next;
// 		// 				p->next = tail;
// 		// 				for (q = tail; q->next != tail; q = q->next)
// 		// 					;
// 		// 
// 		// 				//std::cout << "end" << std::endl;//yogi debug
// 		// 				q->next = 0;
// 		// 				tail->idx = idx;
// 		// 				tail->incidence = incidence;
// 		// 				tail = q;
// 		// 			}
// 		// 
// 		// // 			if (checklist())//tmp
// 		// // 				head = head;
// 		// 		}
// 		else
// 			push_queue(idx, incidence); // push to the queue
// 
// // 		if (checklist())//tmp
// // 			head = head;
// 	}
// 	inline void push_queue(Imgidx idx, int8 incidence)
// 	{
// 		// 		int8 bitfield_idx;
// 		// 		idx = (idx << 1) | (Imgidx)incidence;
// 		// 		bitfield_idx = idx & mask_field;
// 		// 		idx = idx >> shamt;
// 		// 		if (idx < minidx_queue)
// 		// 			minidx_queue = idx;
// 		// 		queue[idx] = queue[idx] | (1 << bitfield_idx);
// 
// 
// 		if (idx < minidx_queue)
// 			minidx_queue = idx;
// 		queue[idx] = incidence;
// 		//qlocalhist[idx >> shamt_qhist]++;
// 	}
// 	inline void pop_queue()
// 	{
// 		Imgidx idx = minidx_queue >> shamt_qhist;
// 
// 		// 		if (minidx_queue == 876) //tmptmptmptmptmp
// 		// 			idx = idx;
// 
// 			//	qlocalhist[idx]--;
// 		queue[minidx_queue] = -1;
// 		//if (qlocalhist[idx])
// 		if (1)
// 		{
// 			minidx_queue = minidx_queue + 1;
// 			while (minidx_queue < qhsize && queue[minidx_queue] == -1)
// 				minidx_queue++;
// 			if (minidx_queue == qhsize)
// 				head = head;
// 		}
// 		else
// 		{
// 			while (!qlocalhist[++idx])
// 				;
// 			minidx_queue = idx << shamt_qhist;
// 			while (queue[minidx_queue] < 0)
// 				minidx_queue++;
// 		}
// 	}
// 	inline void pop()
// 	{
// 		// 		cnt++;//tmp
// 		//  		if (cnt == 3233)//tmp
// 		//  			cnt = cnt;
// 		if (head)
// 			head = head->next;
// 		else
// 			pop_queue();
// 		// 		if (checklist())//tmp
// 		// 			head = head;
// 	}
// 
// 	int8 checklist()//tmp
// 	{
// 		MinList<Imgidx> *p;
// 
// 		if (head)
// 		{
// 			for (p = head; p; p = p->next)
// 			{
// 				if (p->incidence < 0 || p->incidence > 1 || p->next && p > p->next)
// 					return 1;
// 			}
// 		}
// 		if (queue[minidx_queue] < 0 || queue[minidx_queue] > 1)
// 			return 1;
// 		if (minidx_queue < listsize || minidx_queue >= listsize + qsize)
// 			return 1;
// 
// 		uint8 hsum = 0;
// 		for (Imgidx hidx = 0; hidx < hsize; hidx++)
// 		{
// 			Imgidx qidx = hidx << 8;
// 			hsum = 0;
// 			for (int i = 0; i < 256; i++)
// 			{
// 				if (qidx <= qhsize && queue[qidx++] >= 0)
// 					hsum++;
// 			}
// 			if (qlocalhist[hidx] != hsum)
// 				return 1;
// 		}
// 
// 		return 0;
// 		// 		if (head)
// 		// 		{
// 		// 			for (p = head; p; p = p->next)
// 		// 			{
// 		// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
// 		// 					return 1;
// 		// 			}
// 		// 			if (tail->next)
// 		// 				return 1;
// 		// 		}
// 		// 		int n = 0;
// 		// 		for (int i = 0; i < maxSize_list; i++)
// 		// 		{
// 		// 			if (list[i].idx != -1)
// 		// 				n++;
// 		// 		}
// 		// 		if (n != curSize_list)
// 		// 			return 1;
// 		// 		return 0;
// 	}
// 
// 	~HybridQueue()
// 	{
// 		Free(list);
// 		//Free(queue + listsize);
// 		Free(queue);
// 		Free(qlocalhist);
// 	}
// };
// 
// 
// 
