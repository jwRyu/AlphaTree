#pragma once
#include<iostream>
#include "defines.h"
#include "allocator.h"
#include "HQueue.hpp"
#include "Trie.hpp"


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

template<class Imgidx, class Pixel>
struct AlphaNode
{
	Imgidx area;
	Pixel level;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	Imgidx parentidx;
};

template<class Imgidx, class Pixel>
class AlphaTree
{
	inline Pixel abs_diff(Pixel p, Pixel q)
	{
		if (p > q)
			return p - q;
		else
			return q - p;
	}
	void compute_dimg(Pixel* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, i, j;
		
		imgidx = dimgidx = 0;
		if (sizeof(Pixel) < 8)
		{
			if (connectivity == 4)
			{
				for (i = 0; i < height - 1; i++)
				{
					for (j = 0; j < width - 1; j++)
					{						
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
						dhist[dimg[dimgidx++]]++;
						imgidx++;
					}
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					dimgidx++;
					imgidx++;
				}
				for (j = 0; j < width - 1; j++)
				{
					dimgidx++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
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
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					dhist[dimg[dimgidx++]]++;
					imgidx++;
					for (j = 1; j < width - 1; j++)
					{
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
						dhist[dimg[dimgidx++]]++;
						imgidx++;
					}
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					dhist[dimg[dimgidx]]++;
					dimgidx += 3;//skip 2,3
					imgidx++;
				}

				//bottom
				for (j = 0; j < width - 1; j++)
				{
					dimgidx += 3; //skip 0,1,2
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
			}
		}
		else
		{
			if (connectivity == 4)
			{
				for (i = 0; i < height - 1; i++)
				{
					for (j = 0; j < width - 1; j++)
					{
						dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]);
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = abs_diff(img[imgidx + 1], img[imgidx]);
						dhist[dimg[dimgidx++]]++;
						imgidx++;
					}
					dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]);
					dhist[dimg[dimgidx++]]++;
					dimgidx++;
					imgidx++;
				}
				for (j = 0; j < width - 1; j++)
				{
					dimgidx++;
					dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]);
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
					dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]); 
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = abs_diff(img[imgidx + width + 1], img[imgidx]);  
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = abs_diff(img[imgidx + 1], img[imgidx]);
					dhist[dimg[dimgidx++]]++;
					imgidx++;
					for (j = 1; j < width - 1; j++)
					{
						dimg[dimgidx] = abs_diff(img[imgidx + width - 1], img[imgidx]); 
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]); 
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = abs_diff(img[imgidx + width + 1], img[imgidx]);
						dhist[dimg[dimgidx++]]++;
						dimg[dimgidx] = abs_diff(img[imgidx + 1], img[imgidx]);
						dhist[dimg[dimgidx++]]++;
						imgidx++;
					}
					dimg[dimgidx] = abs_diff(img[imgidx + width - 1], img[imgidx]);
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = abs_diff(img[imgidx + width], img[imgidx]);
					dhist[dimg[dimgidx]]++;
					dimgidx += 3;//skip 2,3
					imgidx++;
				}

				//bottom
				for (j = 0; j < width - 1; j++)
				{
					dimgidx += 3; //skip 0,1,2
					dimg[dimgidx] = abs_diff(img[imgidx + 1], img[imgidx]);
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
			}
		}

	}
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
	inline void set_field(uint8* arr, Imgidx idx, uint8 in)
	{
		uint8 shamt = (idx & 1) << 2;
		arr[idx >> 1] &= (in << shamt) | ((uint8)(0x0f) << (4 - shamt));
	}
	inline uint8 get_field(uint8* arr, Imgidx idx)
	{
		return (arr[idx >> 1] >> ((idx & 1) << 2)) & 0x0f;
	}
	inline void push_neighbour(HQueue<Imgidx> *hqueue, Imgidx* levelroot, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		hqueue->push(idx, dissim);
		if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
			levelroot[dissim] = ANODE_CANDIDATE;
#else
			levelroot[dissim] = NewAlphaNode(dissim);
#endif
	}
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode)
	{
		AlphaNode<Imgidx, Pixel> *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->area++;
		pNode->maxPix = max(pNode->maxPix, pix_val);
		pNode->minPix = min(pNode->minPix, pix_val);
		pNode->sumPix += pix_val;
	}
	inline void connectNode2Node(AlphaNode<Imgidx, Pixel>* pPar, Imgidx iPar, AlphaNode<Imgidx, Pixel>* pNode)
	{
		pNode->parentidx = iPar;
		pPar->area += pNode->area;
		pPar->maxPix = max(pNode->maxPix, pPar->maxPix);
		pPar->minPix = min(pNode->minPix, pPar->minPix);
		pPar->sumPix += pNode->sumPix;
	}
	inline Imgidx NewAlphaNode(uint8 level)
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = min(height * width, maxSize + (Imgidx)(height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
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
	inline uint8 is_visited(uint8* isVisited, Imgidx p)
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}
	inline void visit(uint8* isVisited, Imgidx p)
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}
	void Flood_LDR(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		Pixel current_level, next_level;
		Imgidx numlevels;
		HQueue<Imgidx>* hqueue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx iChild, *levelroot;
		uint8 *isVisited, *isAvailable, isAv;
		Imgidx *pParentAry;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << (8 * sizeof(uint8));

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		levelroot = (Imgidx*)Malloc((Imgidx)(numlevels + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		if (connectivity == 4)
			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
			isAvailable = (uint8*)Malloc((size_t)(imgsize));
		for (p = 0; p < numlevels; p++)
			levelroot[p] = NULL_LEVELROOT;
		memset(dhist, 0, (size_t)numlevels * sizeof(Imgidx));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		compute_dimg(dimg, dhist, img);
		dhist[255]++;
		hqueue = new HQueue<Imgidx>(nredges + 1, dhist);

		//tree size estimation (TSE)
		nrmsd = 0;
		for (p = 0; p < numlevels; p++)
			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		//maxSize = imgsize;
		Free(dhist);

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(Imgidx));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		pParentAry = parentAry;

		levelroot[255] = NewAlphaNode(255);
		node[0].parentidx = 0;

		current_level = 255;
		x0 = imgsize >> 1;
		hqueue->push(x0, current_level);

		iChild = 0;

		while (1)
		{
			while (hqueue->min_level <= current_level)
			{
				p = hqueue->pop();
				if (is_visited(isVisited, p))
				{
					hqueue->find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				hqueue->find_min_level();
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
					hqueue->find_min_level();
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
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Imgidx, Pixel>* node;
	Imgidx* parentAry;
	double nrmsd;

	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree() { if (node) { Free(node); Free(parentAry); } }
	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }
	
	void BuildAlphaTree(Pixel *img, Imgidx height, Imgidx width, Imgidx channel, Imgidx connectivity)
	{
		this->height = height;
		this->width = width;
		this->channel = channel;
		this->connectivity = connectivity;
		curSize = 0;
		if (connectivity != 4 && connectivity != 8)
		{
			std::cout << "connectivity should be 4 or 8\n" << std::endl;
			return;
		}


		Flood_LDR(img);
	}
};