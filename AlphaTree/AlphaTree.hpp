#pragma once
#include<iostream>
#include "defines.h"
#include "allocator.h"
#include "HQueue.hpp"


#define DELAYED_NODE_ALLOC		0
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define ANODE_CANDIDATE		0xfffffffe

/*
#define dimg_idx_v(pidx) ((pidx)<<1)
#define dimg_idx_h(pidx) ((pidx)<<1)+1

#define dimg_idx_0(pidx) ((pidx)<<2)
#define dimg_idx_1(pidx) ((pidx)<<2)+1
#define dimg_idx_2(pidx) ((pidx)<<2)+2
#define dimg_idx_3(pidx) ((pidx)<<2)+3

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))
*/
#define A		1.10184
#define SIGMA	-2.6912
#define B		-0.0608
#define M		0.1

//Memory allocation reallocation schemes
/*
#define TSE 0
#define MAXIMUM 1
#define LINEAR 2
#define EXP 3
int mem_scheme = -1;
double size_init[4] = { -1, 1, .2, .15 };
double size_mul[4] = { 1, 1, 1, 2 };
double size_add[4] = { .05, 0, 0.15, 0 };
*/

template<class Tindex, class Tpixel>
struct AlphaNode
{
	Tindex area;
	Tpixel level;  /* alpha of flat zone */
	double sumPix;
	Tpixel minPix;
	Tpixel maxPix;
	Tindex parentidx;
};

template<class Tindex, class Tpixel>
class AlphaTree
{
	void init_isAvailable(uint8* isAvailable)
	{
		int32 i, j, k;
		int32 imgsize = width * height;

		if (connectivity == 4)
		{
			for (i = 0; i < ((imgsize + 1) >> 1); i++)
				isAvailable[i] = 255;

			set_field(isAvailable, 0, 3);
			set_field(isAvailable, width - 1, 5);
			set_field(isAvailable, width * (height - 1), 10);
			set_field(isAvailable, width * height - 1, 12);


			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				set_field(isAvailable, i, 7);
				set_field(isAvailable, j, 14);
				j++;
			}

			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				set_field(isAvailable, j, 11);
				set_field(isAvailable, k, 13);
				j += width;
				k += width;
			}
		}
		else
		{
			for (i = 0; i < imgsize; i++)
				isAvailable[i] = 255;

			isAvailable[0] = 14;
			isAvailable[width - 1] = 131;
			isAvailable[width*(height - 1)] = 56;
			isAvailable[width * height - 1] = 224;


			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				isAvailable[i] = 143;
				isAvailable[j] = 248;
				j++;
			}

			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				isAvailable[j] = 62;
				isAvailable[k] = 227;
				j += width;
				k += width;
			}
		}
	}
	inline uint8 is_available(uint8 isAvailable, uint8 iNeighbour)
	{
		//return	(((isAvailable[idx >> 1] >> ((idx & 1) << 2)) & 0x0f) >> iNeighbour) & 1;
		return	(isAvailable >> iNeighbour) & 1;
	}
	inline void set_field(uint8* arr, Tindex idx, uint8 in)
	{
		uint8 shamt = (idx & 1) << 2;
		arr[idx >> 1] &= (in << shamt) | ((uint8)(0x0f) << (4 - shamt));
	}
	inline uint8 get_field(uint8* arr, Tindex idx)
	{
		return (arr[idx >> 1] >> ((idx & 1) << 2)) & 0x0f;
	}
	inline void push_neighbour(HQueue<Tindex> *hqueue, Tindex* levelroot, Tpixel* dimg, Tindex idx, Tindex dimgidx)
	{
		Tpixel dissim = dimg[dimgidx];
		hqueue->hqueue_push(idx, dissim);
		if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
			levelroot[dissim] = ANODE_CANDIDATE;
#else
			levelroot[dissim] = NewAlphaNode(dissim);
#endif
	}
	void compute_dimg(Tpixel* dimg, Tindex* dhist, Tpixel* img)
	{
		Tindex dimgidx, imgidx, i, j;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width] - (int)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + 1] - (int)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width] - (int)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + 1] - (int)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  -
			//   -  p  3
			//   0  1  2
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				dimgidx++; //skip 0
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width] - (int)img[imgidx]));//1
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width + 1] - (int)img[imgidx]));//2
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//3
				dhist[dimg[dimgidx++]]++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width - 1] - (int)img[imgidx]));//0
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width] - (int)img[imgidx]));//1
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width + 1] - (int)img[imgidx]));//2
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//3
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width - 1] - (int)img[imgidx]));//0
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + width] - (int)img[imgidx]));//1
				dhist[dimg[dimgidx]]++;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				dimg[dimgidx] = (Tpixel)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//3
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
	}
	void Flood(Pixel* img)
	{
		Tindex imgsize, dimgsize, nredges, x0;
		Tpixel current_level, next_level;
		Tindex numlevels;
		HQueue<Tindex>* hqueue;
		Tindex *dhist;
		Tpixel *dimg;
		Tindex iChild, *levelroot;
		uint8 *isVisited, *isAvailable, isAv;
		Tindex *pParentAry;
		Tindex p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
		
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << (8 * sizeof(uint8));

		dhist = (Tindex*)Malloc((size_t)numlevels * sizeof(Tindex));
		dimg = (uint8*)Malloc((size_t)dimgsize * sizeof(uint8));
		levelroot = (Tindex*)Malloc((Tindex)(numlevels + 1) * sizeof(Tindex));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		if (connectivity == 4)
			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
			isAvailable = (uint8*)Malloc((size_t)(imgsize));
		for (p = 0; p < numlevels; p++)
			levelroot[p] = NULL_LEVELROOT;
		memset(dhist, 0, (size_t)numlevels * sizeof(Tindex));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		compute_dimg(dimg, dhist, img);
		dhist[255]++;
		hqueue = new HQueue<Tindex>(nredges + 1, dhist);

		//tree size estimation (TSE)
		nrmsd = 0;
		for (p = 0; p < numlevels; p++)
			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		maxSize = min(imgsize, (Tindex)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		//maxSize = imgsize;
		Free(dhist);

		parentAry = (Tindex*)Malloc((size_t)imgsize * sizeof(Tindex));
		node = (AlphaNode<int32, uint8>*)Malloc((size_t)maxSize * sizeof(AlphaNode<int32, uint8>));
		pParentAry = parentAry;

		levelroot[255] = NewAlphaNode(255);
		node[0].parentidx = 0;

		current_level = 255;
		x0 = imgsize >> 1;
		hqueue->hqueue_push(x0, current_level);

		iChild = 0;

		while (1)
		{
			while (hqueue->min_level <= current_level)
			{
				p = hqueue->hqueue_pop();
				if (is_visited(isVisited, p))
				{
					hqueue->hqueue_find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				hqueue->hqueue_find_min_level();
#endif
				if (connectivity == 4)
				{
					isAv = get_field(isAvailable, p);
					q = p << 1;
					
					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))
						push_neighbour(hqueue, levelroot, dimg, p + width, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))
						push_neighbour(hqueue, levelroot, dimg, p + 1, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))
						push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))
						push_neighbour(hqueue, levelroot, dimg, p - width, q - (width << 1));
						
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))
						push_neighbour(hqueue, levelroot, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))
						push_neighbour(hqueue, levelroot, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))
						push_neighbour(hqueue, levelroot, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))
						push_neighbour(hqueue, levelroot, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))
						push_neighbour(hqueue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))
						push_neighbour(hqueue, levelroot, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))
						push_neighbour(hqueue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))
						push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
				}


				if (current_level > hqueue->min_level)
					current_level = hqueue->min_level;
#if HQUEUE_COST_AMORTIZE
				else
					hqueue->hqueue_find_min_level();
#endif

#if DELAYED_ANODE_ALLOC
				if (levelroot[current_level] == ANODE_CANDIDATE)
					levelroot[current_level] = NewAlphaNode((uint8)current_level);
#endif
				connectPix2Node(p, img[p], levelroot[current_level]);

			}
			//		if(curSize > 22051838 && (curSize))
				//		printf("curSize: %d\n",curSize);
					//Redundant node removal
			if (node[iChild].parentidx == levelroot[current_level] &&
				node[levelroot[current_level]].area == node[iChild].area)
			{
				levelroot[current_level] = iChild;
#if DELAYED_ANODE_ALLOC
				curSize--;
				//memset((uint8*)(node + curSize), 0, sizeof(AlphaNode));
#endif
			}

			next_level = current_level + 1;
			while (~next_level && levelroot[next_level] == NULL_LEVELROOT)
				next_level++;
#if DELAYED_ANODE_ALLOC
			if (levelroot[next_level] == ANODE_CANDIDATE)
				levelroot[next_level] = NewAlphaNode((uint8)next_level);
#endif
			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			if (node[levelroot[next_level]].area == imgsize)
			{
				node[levelroot[next_level]].parentidx = 0;
				break;
			}

			iChild = levelroot[current_level];
			levelroot[current_level] = NULL_LEVELROOT;
			current_level = next_level;
		}


		delete hqueue;
		Free(dimg);
		Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}
public:
	Tindex maxSize;
	Tindex curSize;
	Tindex height, width, channel, connectivity;
	AlphaNode<int32, uint8>* node;
	Tindex* parentAry;
	double nrmsd;

	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree() { if (node) { Free(node); Free(parentAry); } }
	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }
	
	inline void connectPix2Node(Tindex pidx, Tpixel pix_val, Tindex iNode)
	{
		AlphaNode<int32, uint8> *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->area++;
		pNode->maxPix = max(pNode->maxPix, pix_val);
		pNode->minPix = min(pNode->minPix, pix_val);
		pNode->sumPix += pix_val;
	}

	inline void connectNode2Node(AlphaNode<int32, uint8>* pPar, Tindex iPar, AlphaNode<int32, uint8>* pNode)
	{
		pNode->parentidx = iPar;
		pPar->area += pNode->area;
		pPar->maxPix = max(pNode->maxPix, pPar->maxPix);
		pPar->minPix = min(pNode->minPix, pPar->minPix);
		pPar->sumPix += pNode->sumPix;
	}

	inline Tindex NewAlphaNode(uint8 level)
	{
		AlphaNode<int32, uint8> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout<<"Reallocating...\n";
			maxSize = min(height * width, maxSize + (Tindex)(height * width * 0.1));

			node = (AlphaNode<int32, uint8>*)Realloc(node, maxSize * sizeof(AlphaNode<int32, uint8>));
			pNew = node + curSize;
		}
		pNew->level = level;
		pNew->minPix = (uint8)-1;
		pNew->maxPix = 0;
		pNew->sumPix = 0.0;
		pNew->parentidx = 0;
		pNew->area = 0;

		return curSize++;
	}

	inline uint8 is_visited(uint8* isVisited, Tindex p)
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}

	inline void visit(uint8* isVisited, Tindex p)
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

	void BuildAlphaTree(Tpixel *img, Tindex height, Tindex width, Tindex channel, Tindex connectivity)
	{
		this->height = height;
		this->width = width;
		this->channel = channel;
		this->connectivity = connectivity;
		curSize = 0;
		if (connectivity != 4 && connectivity != 8)
		{
			std::cout << "connectivity should be 4 or 8\n" << std::endl;
		}
		else
			Flood(img);
	}
};