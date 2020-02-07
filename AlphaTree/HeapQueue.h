
#define PQ_DEBUG 0
#if PQ_DEBUG
#include <iostream>
#include <fstream>
#endif

template<class Imgidx, class Pixel>
class HQentry
{
public:
	HQentry() {}
	~HQentry() {}
	Imgidx pidx;
	Pixel alpha;
	inline void operator=(const HQentry& item)
	{
		this->pidx = item.pidx;
		this->alpha = item.alpha;
	}
};


//Heap-based priority queue
template<class Imgidx, class Pixel>
class HeapQueue
{
	Imgidx cursize;
	Imgidx maxsize;
	HQentry<Imgidx, Pixel> *arr;
	Pixel pop_level;
	Pixel max_level;
	HQentry<Imgidx, Pixel> pushlist[7];//Connectivity - 1
	int8 pushlistidx;
#if PQ_DEBUG
	std::fstream fs;
#endif
public:
//	Pixel min_level;
	HeapQueue(Imgidx maxsize) : cursize(0), maxsize(maxsize)
	{
		arr = new HQentry<Imgidx, Pixel>[maxsize + 1];
		pushlistidx = 0;
#if PQ_DEBUG
		this->fs.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", std::fstream::out);
#endif
		max_level = (Pixel)-1;
	}

	~HeapQueue()
	{
		delete[] arr;

#if PQ_DEBUG
		fs.close();
#endif
	}

	inline int minidx_pushlist()
	{
		int minidx = 0;
		for (int i = 0; i < pushlistidx; i++)
		{
			if (pushlist[i].alpha < pushlist[minidx].alpha)
				minidx = i;
		}
		return minidx;
	}

	inline void find_minlev()
	{
		int minidx = minidx_pushlist();
		if (pushlist[minidx].alpha <= arr[1].alpha)
		{
			arr[1].alpha = pushlist[minidx].alpha;
			arr[1].pidx = pushlist[minidx].pidx;
		}
		else
		{
			pop_run();
			push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
		}
		
		for (int i = 1; i < pushlistidx; i++)
		{
			pushlist[minidx].alpha = max_level;
			minidx = minidx_pushlist();
			push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
		}
//		min_level = arr[1].alpha;
	}

	inline Pixel get_minlev() { return arr[1].alpha; }

	inline void validate()
	{
		Imgidx i, j;

		for (i = cursize; i > 1; i--)
		{
			j = i >> 1;
			if (arr[i].alpha < arr[j].alpha)
			{
				std::cout << "invalidate queue" << std::endl;
				std::cin >> i;
			}
		}
	}

	inline Imgidx top() { return arr[1].pidx; }

	inline Imgidx pop()
	{
		//pop_level = arr[1].alpha;
		pushlistidx = 0;
		return arr[1].pidx;
	}
	inline Imgidx pop_run()
	{
		Imgidx outval = arr[1].pidx;
		Imgidx current = 1, next, next0, next1, curidx;
		Pixel curalpha;

		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;
		

		curidx = arr[cursize].pidx;
		curalpha = arr[cursize].alpha;
//		if (curidx < 0)
	//		curidx = curidx;
		cursize--;

#if PQ_DEBUG
		fs << "1" << std::endl << (int)arr[1].alpha << std::endl;
#endif
		while (1)
		{
			next0 = current << 1;
			next1 = next0 + 1;
			if (next0 > cursize)
				break;
			if (next1 <= cursize && arr[next1].alpha < arr[next0].alpha)
				next = next1;
			else
				next = next0;

			if (curalpha < arr[next].alpha)
				break;

			arr[current] = arr[next];
			current = next;
		}
		arr[current].alpha = curalpha;
		arr[current].pidx = curidx;
// 
// 		if (cursize)
// 			min_level = arr[1].alpha;
// 		else
// 			min_level = (Pixel)-1;

		//validate();

		return outval;
	}

	inline void push(Imgidx pidx, Pixel alpha)
	{
		pushlist[pushlistidx].pidx = pidx;
		pushlist[pushlistidx++].alpha = alpha;
	}

// 	void push(Imgidx pidx)
// 	{
// 		push_run(pidx, (Pixel)pidx);
// 	}

	inline void push_run(Imgidx pidx, Pixel alpha)
	{
		Imgidx current, next;

#if PQ_DEBUG
		fs << "0" << std::endl << (int)alpha << std::endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		current = cursize;

	//	if (pidx < 0)
	//		pidx = pidx;

		next = current >> 1;
		while (next && (arr[next].alpha > alpha))
		{
			arr[current] = arr[next];
			current = next;
			next = next >> 1;
		}

		arr[current].pidx = pidx;
		arr[current].alpha = alpha;

// 		if (current == 1)
// 			min_level = alpha;
		//validate();
	}
};


//Heap-based priority queue
template<class Imgidx>
class HeapQueue_rank
{
	Imgidx cursize;
	Imgidx maxsize;
	Imgidx *arr;
#if PQ_DEBUG
	std::fstream fs;
#endif
public:
	Imgidx min_level;
	HeapQueue_rank(Imgidx maxsize) : cursize(0), maxsize(maxsize)
	{
		arr = new Imgidx[maxsize + 1];
#if PQ_DEBUG
		this->fs.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", std::fstream::out);
#endif

	}
	~HeapQueue_rank()
	{
		delete[] arr;

#if PQ_DEBUG
		fs.close();
#endif
	}

	inline Imgidx get_minlev() { return arr[1]; }
	inline void find_minlev() {};
	inline void validate()
	{
		Imgidx i, j;

		for (i = cursize; i > 1; i--)
		{
			j = i >> 1;
			if (arr[i].alpha < arr[j].alpha)
			{
				std::cout << "invalidate queue" << std::endl;
				std::cin >> i;
			}
		}
	}
	inline Imgidx top() { return arr[1]; }
	inline Imgidx pop()
	{
		Imgidx outval;
		Imgidx current = 1, next, next0, next1, curidx;
		Imgidx curalpha;

		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;


		outval = arr[current];
		curidx = arr[cursize--];
		//curalpha = arr[cursize].alpha;
		//		if (curidx < 0)
			//		curidx = curidx;

#if PQ_DEBUG
		fs << "1" << std::endl << (int)arr[1].alpha << std::endl;
#endif
		while (1)
		{
			next0 = current << 1;
			next1 = next0 + 1;
			if (next0 > cursize)
				break;
			if (next1 <= cursize && arr[next1] < arr[next0])
				next = next1;
			else
				next = next0;

			if (curidx < arr[next])
				break;

			arr[current] = arr[next];
			current = next;
		}
		arr[current] = curidx;

		if (cursize)
			min_level = arr[1];
		else
			min_level = (Imgidx)-1;

		//validate();

		return outval;
	}
	inline void push(Imgidx pidx)
	{
		Imgidx current, next;

#if PQ_DEBUG
		fs << "0" << std::endl << (int)alpha << std::endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		current = cursize;

		//	if (pidx < 0)
		//		pidx = pidx;

		next = current >> 1;
		while (next && (arr[next] > pidx))
		{
			arr[current] = arr[next];
			current = next;
			next = next >> 1;
		}

		arr[current] = pidx;

		if (current == 1)
			min_level = pidx;
		//validate();
	}
};