
#define HEAPQ_STAT 1
#if HEAPQ_STAT
#include <iostream>
#include <fstream>
using namespace std;
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
#if HEAPQ_STAT
	ofstream f;
	Imgidx popcnt;
#endif
public:
//	Pixel min_level;
	HeapQueue(Imgidx maxsize) : cursize(0), maxsize(maxsize)
	{
		arr = new HQentry<Imgidx, Pixel>[maxsize + 1];
		pushlistidx = 0;
#if HEAPQ_STAT
		this->f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", ofstream::out);
		popcnt = 0;
#endif
		max_level = (Pixel)-1;
	}

	~HeapQueue()
	{
		delete[] arr;

#if HEAPQ_STAT
		f.close();
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

#if HEAPQ_STAT
		f << "1" << std::endl << (int)arr[1].alpha << std::endl;
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

#if HEAPQ_STAT
		f << "0" << std::endl << (int)alpha << std::endl;
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
	
	Imgidx *arr;
#if HEAPQ_STAT
	ofstream f;
	Imgidx *depth;
	Imgidx cur_depth;
	Imgidx popcnt;
#endif
public:
	Imgidx cursize;
	Imgidx maxsize;
	Imgidx min_level;
	HeapQueue_rank(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in)
	{
		arr = new Imgidx[maxsize];
		arr--;
#if HEAPQ_STAT
		popcnt = 0;
		f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", std::ofstream::app);
		f << -1 << '\n' << maxsize << endl;
		depth = (Imgidx*)Calloc(maxsize * sizeof(Imgidx));
#endif

	}
	~HeapQueue_rank()
	{
		delete[] (arr + 1);

#if HEAPQ_STAT
		f.close();
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
#if HEAPQ_STAT
		Imgidx curdepth = depth[cursize + 1];
		cur_depth = cur_depth > depth[1] ? cur_depth : depth[1];
		f << '1' << std::endl << arr[1] << std::endl << cur_depth++ << std::endl;
		if (popcnt++ % 3000 == 0)
			cout << "pop " << popcnt << '/' << maxsize << endl;
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
#if HEAPQ_STAT
			depth[current] = depth[next];
#endif
			current = next;
		}
		arr[current] = curidx;
#if HEAPQ_STAT
		depth[current] = curdepth;
		f << arr[1] - min_level << endl;
#endif

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

#if HEAPQ_STAT
		f << '0' << '\n' << pidx << endl;
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
#if HEAPQ_STAT
			depth[current] = depth[next];
#endif
			current = next;
			next = next >> 1;
		}

		arr[current] = pidx;
#if HEAPQ_STAT
		depth[current] = 0;
#endif

		if (current == 1)
		{
			min_level = pidx;
#if HEAPQ_STAT
			if (cursize > 2)
			{
				if (arr[2] < arr[3])
					depth[2] = depth[2] > cur_depth? depth[2] : cur_depth;
				else
					depth[3] = depth[3] > cur_depth ? depth[3] : cur_depth;
			}
			else if (cursize > 1)
				depth[2] = depth[2] > cur_depth ? depth[2] : cur_depth;
			cur_depth = 0;
#endif
		}
		//validate();
	}

	//for hybrid queueing
// 	inline Imgidx get_arr(Imgidx arridx) { return arr[arridx]; }
// 	inline void set_arr(Imgidx arridx, Imgidx val) { arr[arridx] = val; }
};