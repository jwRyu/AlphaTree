
#define PQ_DEBUG 0
#if PQ_DEBUG
#include <iostream>
#include <fstream>
#endif

template<class Imgidx, class Pixel>
class PQarray
{
public:
	PQarray() {}
	~PQarray() {}
	Imgidx pidx;
	Pixel alpha;
	inline void operator=(const PQarray& item)
	{
		this->pidx = item.pidx;
		this->alpha = item.alpha;
	}
};


//Heap-based priority queue
template<class Imgidx, class Pixel>
class PriorityQueue
{
	Imgidx cursize;
	Imgidx maxsize;
	PQarray<Imgidx, Pixel> *arr;
#if PQ_DEBUG
	std::fstream fs;
#endif
public:
	Pixel min_level;
	PriorityQueue(Imgidx maxsize) : cursize(0), maxsize(maxsize)
	{
		arr = new PQarray<Imgidx, Pixel>[maxsize + 1];
#if PQ_DEBUG
		this->fs.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", std::fstream::out);
#endif

	}

	~PriorityQueue()
	{
		delete[] arr;

#if PQ_DEBUG
		fs.close();
#endif
	}

	inline void find_min_level() {};

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


	Imgidx pop()
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

		if (cursize)
			min_level = arr[1].alpha;
		else
			min_level = (Pixel)-1;

		//validate();

		return outval;
	}

	void push(Imgidx pidx, Pixel alpha)
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

		if (current == 1)
			min_level = alpha;
		//validate();
	}
};