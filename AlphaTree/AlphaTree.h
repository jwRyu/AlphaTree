#pragma once
#include<iostream>
#include "defines.h"
#include "allocator.h"
#include "HQueue.h"
#include "Trie.h"



#define DELAYED_NODE_ALLOC		1
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define NODE_CANDIDATE		0xfffffffe

#define dimg_idx_v(pidx) ((pidx)<<1)
#define dimg_idx_h(pidx) ((pidx)<<1)+1

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

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

#define A		1.3901
#define SIGMA	-2.1989
#define B		-0.1906
#define M		0.05

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
class AlphaNode
{
public:
	Imgidx area;
	Pixel level;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	Imgidx parentidx;

	inline void set(Imgidx area, Pixel level, double sumPix, Pixel minPix, Pixel maxPix)
	{
		this->area = area;
		this->level = level;
		this->sumPix = sumPix;
		this->minPix = minPix;
		this->maxPix = maxPix;
	}
	inline void add(AlphaNode* q)
	{
		this->area += q->area;
		this->sumPix += q->sumPix;
		this->minPix = min(this->minPix, q->minPix);
		this->maxPix = max(this->maxPix, q->maxPix);
	}
	inline void add(Pixel pix_val)
	{
		this->area++;
		this->sumPix += (double)pix_val;
		this->minPix = min(this->minPix, pix_val);
		this->maxPix = max(this->maxPix, pix_val);
	}
	inline void copy(AlphaNode* q)
	{
		this->area = q->area;
		this->sumPix = q->sumPix;
		this->minPix = q->minPix;
		this->maxPix = q->maxPix;
	}
};

template<class Imgidx, class Pixel>
class RankItem
{
public:
	Pixel alpha;
	Imgidx dimgidx;

	inline void operator=(const RankItem& q)
	{
		this->alpha = q.alpha;
		this->dimgidx = q.dimgidx;
	}
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
	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>* rankinfo, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;

		if (connectivity == 4)
			nredges = (width - 1) * height + width * (height - 1);
		else
			nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);

		//rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		tmp = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc((sizeof(Pixel) >> 1) * 65536 * sizeof(Imgidx));

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
					rankinfo[contidx++].dimgidx = dimgidx++; 
					imgidx++;
				}
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 2;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
				rankinfo[contidx++].dimgidx = dimgidx++;
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
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
			}
		}

		for (i = 0; i < (sizeof(Pixel) >> 1) * 65536; i++)
			hist[i] = 0;
		for (i = 0; i < nredges; i++)
		{
			hidx = rankinfo[i].alpha;
			hist[hidx & mask]++;

			shamt = (Pixel)16;
			h_offset = (Pixel)65536;
			for (nbyte = 2; nbyte < sizeof(Pixel); nbyte += 2)
			{
				hidx = hidx >> shamt;
				hist[h_offset + (hidx & mask)]++;
				h_offset += (Pixel)65536;
			}
		}

		h = hist;
		shamt = 0;
		for (nbyte = 0; nbyte < sizeof(Pixel); nbyte += 2)
		{
			hsum = 0;
			for (i = 0; i < hist_size; i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = (rankinfo[i].alpha >> shamt) & mask;
				j = --h[hidx];
				tmp[j] = rankinfo[i];
			}
			r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
			shamt += 16;
			h += hist_size;
		}

		for (i = 0; i < nredges; i++)
		{
			//rank_to_alpha[i] = rankinfo[i].alpha;
			rank[rankinfo[i].dimgidx] = i;
		}

		//Free(rankinfo);
		Free(tmp);
		Free(hist);
	}
	void compute_dimg(Pixel* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + stride_w] - (int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + stride_w] - (int64)img[imgidx]));
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
	void init_isAvailable(uint8* isAvailable)
	{
		int32 i, j, k;
		int32 imgsize = width * height;

		if (connectivity == 4)
		{

			//		    Neighbour Index
			// 			-      3      -
			// 			2    pixel    1
			// 			-      0      -
			//          
			//			Neighbour indices to bit field
			//			7 6 5 4 3 2 1 0
			//         MSB			 LSB
			//			0: Neighbour pixel not available (corner of Image or partition)
			//			1: available
			for (i = 0; i < ((imgsize + 1) >> 1); i++)
				isAvailable[i] = 0xff;

			set_field(isAvailable, 0, 0x3);
			set_field(isAvailable, width - 1, 0x5);
			set_field(isAvailable, width * (height - 1), 0xa);
			set_field(isAvailable, width * height - 1, 0xc);


			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				set_field(isAvailable, i, 0x7);
				set_field(isAvailable, j, 0xe);
				j++;
			}

			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				set_field(isAvailable, j, 0xb);
				set_field(isAvailable, k, 0xd);
				j += width;
				k += width;
			}
		}
		else
		{
//		    Neighbour Index
// 			6      5      4
// 			7    pixel    3
// 			0      1      2
//          
//			Neighbour indices to bit field
//			7 6 5 4 3 2 1 0
//         MSB			 LSB
//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
//			1: available

			//initialize to all available
			for (i = 0; i < imgsize; i++)
				isAvailable[i] = 0xff;
			
			//four corners
			isAvailable[0] = 0x0e;
			isAvailable[width - 1] = 0x83;
			isAvailable[width*(height - 1)] = 0x38;
			isAvailable[width * height - 1] = 0xe0;

			//top and bottom row
			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				isAvailable[i] = 0x8f;
				isAvailable[j] = 0xf8;
				j++;
			}

			//leftest and rightest column
			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				isAvailable[j] = 0x3e;
				isAvailable[k] = 0xe3;
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
#if DELAYED_NODE_ALLOC
			levelroot[dissim] = NODE_CANDIDATE;
#else
			levelroot[dissim] = NewAlphaNode(dissim);
#endif
	}	
	inline void connectPix2Node(Imgidx* parentAry, Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		if (pNode->area) //possibly unnecessary branch..
			pNode->add(pix_val);			
		else
			pNode->set(1, level, (double)pix_val, pix_val, pix_val);
	}

#if DELAYED_NODE_ALLOC
	inline void connectPix2Node(Imgidx* parentAry, Imgidx pidx, Pixel pix_val, Imgidx *levelroot, int32 level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		Imgidx iNode = levelroot[level];
		if (iNode == NODE_CANDIDATE)
		{
			iNode = NewAlphaNode();
			levelroot[level] = iNode;
			parentAry[pidx] = iNode;
			pNode = node + iNode;

			pNode->set(1, level, (double)pix_val, pix_val, pix_val);
		}
		else
		{
			pNode = node + iNode;
			parentAry[pidx] = iNode;
			pNode->add(pix_val);
		}
	}
#else
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode)
	{
		AlphaNode<Imgidx, Pixel> *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->add(pix_val);
	}
#endif
	inline void connectNode2Node(Imgidx iChild, Imgidx iPar, Pixel level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		pChild = node + iChild;
		pPar = node + iPar;
		pChild->parentidx = iPar;
		if (pPar->area)
			pPar->add(pChild);
		else
		{
			pPar->level = level;
			pPar->copy(pChild);
		}
	}
#if DELAYED_NODE_ALLOC
	inline void connectNode2Node(Imgidx* levelroot, Imgidx iChild, int32 level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		Imgidx iPar = levelroot[level];
		if (iPar == NODE_CANDIDATE)
		{
			iPar = NewAlphaNode();
			levelroot[level] = iPar;
			pPar = node + iPar;
			pChild = node + iChild;
			pChild->parentidx = iPar;
			pPar->level = level;
			pPar->copy(pChild);
		}
		else
		{
			pPar = node + iPar;
			pChild = node + iChild;
			pChild->parentidx = iPar;
			pPar->add(pChild);
		}
	}
#else
	inline void connectNode2Node(AlphaNode<Imgidx, Pixel>* pPar, Imgidx iPar, AlphaNode<Imgidx, Pixel>* pNode)
	{
		pNode->parentidx = iPar;
		pPar->add(pNode);
	}
#endif
#if DELAYED_NODE_ALLOC
	inline Imgidx NewAlphaNode()
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = min(2 * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			pNew = node + curSize;
		}
		return curSize++;
	}
#else
	inline Imgidx NewAlphaNode(Pixel level) //Fix it later - no need to initialize
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
		//		pNew->parentidx = 0;
		pNew->area = 0;

		return curSize++;
	}
#endif
	inline uint8 is_visited(uint8* isVisited, Imgidx p)
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}
	inline void visit(uint8* isVisited, Imgidx p)
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

	void Flood_HQueue(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int32 numlevels, max_level, current_level, next_level;
		HQueue<Imgidx>* hqueue;
		Imgidx *dhist;
		Pixel *dimg;
		Pixel dissim;
		Imgidx iChild, *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << (8 * sizeof(uint8));

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		levelroot = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		if (connectivity == 4)
			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
			isAvailable = (uint8*)Malloc((size_t)(imgsize));
		for (p = 0; p < numlevels; p++)
			levelroot[p] = NULL_LEVELROOT;
		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		max_level = (uint8)(numlevels - 1);
		
		compute_dimg(dimg, dhist, img);

		dhist[max_level]++;
		hqueue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels);
		curSize = 0;

		/////////////////////////////////////////
		//tree size estimation (TSE)
		nrmsd = 0;
		for (p = 0; p < numlevels; p++)
			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
		/////////////////////////////////////////
	//	printf("NRMSD: %f\tEst. NTS: %f\tEst. Tree size: %d\n", nrmsd, ((A * exp(SIGMA * nrmsd) + B) + M), tree->maxSize);
		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
		Free(dhist);

		parentAry = (int32*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		
#if DELAYED_NODE_ALLOC
		levelroot[max_level + 1] = NODE_CANDIDATE;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = imgsize >> 1;
		hqueue->push(x0, current_level);

		iChild = levelroot[max_level + 1];
		while (current_level <= max_level)
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

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(hqueue, levelroot, dimg, p + width, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))			push_neighbour(hqueue, levelroot, dimg, p + 1, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))			push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(hqueue, levelroot, dimg, p - width, q - (width << 1));
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(hqueue, levelroot, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))			push_neighbour(hqueue, levelroot, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(hqueue, levelroot, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))				push_neighbour(hqueue, levelroot, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(hqueue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))			push_neighbour(hqueue, levelroot, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(hqueue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))				push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
				}

				if (current_level > hqueue->min_level)
					current_level = hqueue->min_level;
#if HQUEUE_COST_AMORTIZE
				else
					hqueue->find_min_level();
#endif

#if DELAYED_NODE_ALLOC
				connectPix2Node(parentAry, p, img[p], levelroot, current_level);
#else
				connectPix2Node(parentAry, p, img[p], node + levelroot[current_level], levelroot[current_level]);
#endif

			}

			if (node[iChild].parentidx == levelroot[current_level] &&
				node[levelroot[current_level]].area == node[iChild].area)
			{
				levelroot[current_level] = iChild;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			next_level = current_level + 1;
			while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
				next_level++;

#if DELAYED_NODE_ALLOC
			connectNode2Node(levelroot, levelroot[current_level], next_level);
#else
			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
#endif


			iChild = levelroot[current_level];
			levelroot[current_level] = NULL_LEVELROOT;
			current_level = next_level;

		}
		node[iChild].parentidx = iChild;
		rootidx = iChild;
		curSize--;

		for (p = 0; p < imgsize; p++)
		{
			if (node[parentAry[p]].level)//Singleton 0-CC
			{
				x0 = NewAlphaNode();
				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
				node[x0].parentidx = parentAry[p];
				parentAry[p] = x0;
			}
		}

		delete hqueue;
		Free(dimg);
		Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}

// 	void Flood_HQueue(Pixel* img)
// 	{
// 		Imgidx imgsize, dimgsize, nredges, x0;
// 		int16 current_level, next_level, numlevels, max_level;
// 		HQueue<Imgidx>* hqueue;
// 		Imgidx *dhist;
// 		Pixel *dimg;
// 		Imgidx iChild, *levelroot;
// 		uint8 *isVisited, *isAvailable, isAv;
// 		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
// 
// 		imgsize = width * height;
// 		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
// 		dimgsize = (connectivity >> 1) * width * height;
// 		numlevels = 1 << (8 * sizeof(uint8));
// 		max_level = numlevels - 1;
// 
// 		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
// 		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
// 		levelroot = (Imgidx*)Malloc((Imgidx)(numlevels + 1) * sizeof(Imgidx));
// 		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
// 		if (connectivity == 4)
// 			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
// 		else
// 			isAvailable = (uint8*)Malloc((size_t)(imgsize));
// 		for (p = 0; p < numlevels; p++)
// 			levelroot[p] = NULL_LEVELROOT;
// 		memset(dhist, 0, (size_t)numlevels * sizeof(Imgidx));
// 		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
// 		init_isAvailable(isAvailable);
// 
// 		compute_dimg(dimg, dhist, img);
// 		dhist[max_level]++;
// 		hqueue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels);
// 
// 		//tree size estimation (TSE)
// 		nrmsd = 0;
// 		for (p = 0; p < numlevels; p++)
// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
// 		maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
// 		//maxSize = imgsize;
// 		Free(dhist);
// 
// 		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(Imgidx));
// 		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
// 
// 
// #if DELAYED_NODE_ALLOC
// 		levelroot[numlevels] = ANODE_CANDIDATE;
// #else
// 		levelroot[numlevels - 1] = NewAlphaNode((uint8)max_level);
// 		node[levelroot[numlevels - 1]].parentidx = levelroot[numlevels - 1];
// #endif
// 
// 		current_level = max_level;
// 		x0 = imgsize >> 1;
// 		hqueue->push(x0, current_level);
// 		iChild = 0;
// 		while (current_level <= max_level)
// 		{
// 			while (hqueue->min_level <= current_level)
// 			{
// 				p = hqueue->pop();
// 				if (is_visited(isVisited, p))
// 				{
// 					hqueue->find_min_level();
// 					continue;
// 				}
// 				visit(isVisited, p);
// #if !HQUEUE_COST_AMORTIZE
// 				hqueue->find_min_level();
// #endif
// 				if (connectivity == 4)
// 				{
// 					isAv = get_field(isAvailable, p);
// 					q = p << 1;
// 
// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(hqueue, levelroot, dimg, p + width, q);
// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))			push_neighbour(hqueue, levelroot, dimg, p + 1, q + 1);
// 					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))			push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
// 					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(hqueue, levelroot, dimg, p - width, q - (width << 1));
// 
// 				}
// 				else
// 				{
// 					isAv = isAvailable[p];
// 					q = p << 2;
// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(hqueue, levelroot, dimg, p + wstride1, q);
// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))			push_neighbour(hqueue, levelroot, dimg, p + width, q + 1);
// 					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(hqueue, levelroot, dimg, p + wstride0, q + 2);
// 					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))				push_neighbour(hqueue, levelroot, dimg, p + 1, q + 3);
// 					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(hqueue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
// 					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))			push_neighbour(hqueue, levelroot, dimg, p - width, q - wstride_d + 1);
// 					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(hqueue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
// 					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))				push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
// 				}
// 
// 
// 				if (current_level > hqueue->min_level)
// 					current_level = hqueue->min_level;
// #if HQUEUE_COST_AMORTIZE
// 				else
// 					hqueue->find_min_level();
// #endif
// 
// #if DELAYED_NODE_ALLOC
// 				connectPix2Node(parentAry, p, img[p], levelroot, current_level);
// #else
// 				connectPix2Node(p, img[p], levelroot[current_level]);
// #endif
// 
// 			}
// 			//		if(curSize > 22051838 && (curSize))
// 				//		printf("curSize: %d\n",curSize);
// 					//Redundant node removal
// 
// 			if (node[iChild].parentidx == levelroot[current_level] &&
// 				node[levelroot[current_level]].area == node[iChild].area)
// 			{
// 				levelroot[current_level] = iChild;
// #if DELAYED_NODE_ALLOC
// 				curSize--;
// 				//memset((uint8*)(node + curSize), 0, sizeof(AlphaNode));
// #endif
// 			}
// 			next_level = current_level + 1;
// 			while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
// 				next_level++;
// #if DELAYED_NODE_ALLOC
// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
// #else
// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + currentlvlroot);
// #endif
// 
// 			iChild = levelroot[current_level];
// 			levelroot[current_level] = NULL_LEVELROOT;
// 			current_level = next_level;
// 		}
// 		node[iChild].parentidx = iChild;
// 		curSize--;
// 		rootidx = curSize - 1;
// 		
// 		for (p = 0; p < imgsize; p++)
// 		{
// 			if (node[parentAry[p]].level)//Singleton 0-CC
// 			{
// 				x0 = NewAlphaNode();
// 				node[x0].level = 0;
// 				node[x0].area = 1;
// 				node[x0].maxPix =
// 				node[x0].minPix = img[p];
// 				node[x0].sumPix = (double)img[p];
// 				node[x0].parentidx = parentAry[p];
// 				parentAry[p] = x0;
// 			}
// 		}		
// 		delete hqueue;
// 		Free(dimg);
// 		Free(levelroot);
// 		Free(isVisited);
// 		Free(isAvailable);
// 	}
	void Flood_Trie(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		Imgidx current_rank, next_rank;
		Trie<Imgidx, trieidx> *trie;
		RankItem<Imgidx, Pixel>* rankinfo;
		Imgidx dimgidx;
		Imgidx *rank, top_rank;
		int8 incidence, shamt, mask;
//		Pixel *rank2alpha;
		Imgidx iChild, *levelroot;
		uint8 *isVisited, *isAvailable, isAv;
		Imgidx trietop, p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
		Imgidx incidence_map[4];

		if (connectivity == 4)
		{
			incidence_map[0] = width; 
			incidence_map[1] = 1;		
			shamt = 1;
			mask = 1;
		}
		else
		{
			incidence_map[0] = width - 1;
			incidence_map[1] = width;     
			incidence_map[2] = width + 1; 
			incidence_map[3] = 1;		  
			shamt = 2;
			mask = 3;
		}

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		
		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		if (connectivity == 4)
			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
			isAvailable = (uint8*)Malloc((size_t)(imgsize));
		
		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; // 
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		compute_dimg(rank, rankinfo, img);
		trie = new Trie<Imgidx, trieidx>(nredges << 1);
		//trie->push(nredges, 0);


// 		//tree size estimation (TSE)
// 		
// 		nrmsd = 0;
// 		for (p = 0; p < numlevels; p++)
// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
// 		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
// 		//maxSize = imgsize;
// 		Free(dhist);
		maxSize = nredges + imgsize;

		parentAry = (int32*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		for (p = 0; p < maxSize; p++)
		{
			node[p].parentidx = -1;
			node[p].area = 0;
		}

		visit(isVisited, 0);
		if (connectivity == 4)
		{
			trie->push(rank[0], 0);
			trie->push(rank[1], 0);
		}
		else
		{
			trie->push(rank[1], 0);
			trie->push(rank[2], 0);
			trie->push(rank[3], 0);
		}

		current_rank = trie->top() >> 1;
		connectPix2Node(parentAry, 0, img[0], current_rank, rankinfo[current_rank].alpha);
		//x0 = imgsize >> 1;
		iChild = current_rank;

		while (1)//(current_rank <= nredges)
		{
			while (1)//((trie->top() >> 1) <= current_rank)
			{
				trietop = trie->top();		//remove tmp variables later if possible
				incidence = trietop & 1;	//0 is outgoing, 1 is incoming
				top_rank = trietop >> 1;	//remove tmp variables later if possible
				dimgidx = rankinfo[top_rank].dimgidx;
				p = (dimgidx >> shamt) + (incidence_map[dimgidx & mask] & (incidence - 1)); //current pixel idx

				if (is_visited(isVisited, p))//yogi
					break;
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				//find_min_level();
#endif
				if (connectivity == 4)
				{
					isAv = get_field(isAvailable, p);
					q = p << 1;
					
					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		
						trie->push(rank[q], 0); 
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))			
						trie->push(rank[q + 1], 0);
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))			
						trie->push(rank[q - 1], 1);
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		
						trie->push(rank[q - (width << 1)], 1);
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))
						trie->push(rank[q], 0);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))
						trie->push(rank[q + 1], 0);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))
						trie->push(rank[q + 2], 0);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))
						trie->push(rank[q + 3], 0);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))
						trie->push(rank[q - wstride_d + 4], 1);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))
						trie->push(rank[q - wstride_d + 1], 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))
						trie->push(rank[q - wstride_d - 2], 1);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))
						trie->push(rank[q - 1], 1);
				}
		
				if (current_rank > trie->min_rank())
				{
					current_rank = trie->min_rank(); 
					connectPix2Node(parentAry, p, img[p], current_rank, rankinfo[current_rank].alpha);
				}
				else
				{
					connectPix2Node(parentAry, p, img[p], current_rank, rankinfo[current_rank].alpha);
					break;
				}				
			}

			//Redundant node removal
// 			if (node[iChild].parentidx == levelroot[current_rank] &&
// 				node[levelroot[current_rank]].area == node[iChild].area)
// 			{
// 				levelroot[current_rank] = iChild;
// #if DELAYED_NODE_ALLOC
// 				curSize--;
// 				//memset((uint8*)(node + curSize), 0, sizeof(AlphaNode));
// #endif
// 			}

			trie->pop();
			next_rank = trie->top() >> 1;

			//Redundant node removal
			if (node[iChild].parentidx == current_rank &&
				node[iChild].area == node[current_rank].area)
				current_rank = iChild;
			
			connectNode2Node(current_rank, next_rank, rankinfo[next_rank].alpha);
// #if DELAYED_NODE_ALLOC
// 			connectNode2Node(current_rank, next_rank, rankinfo[current_rank].alpha);
// #else
// 			connectNode2Node(node + levelroot[next_rank], levelroot[next_rank], node + levelroot[current_rank]);
//#endif
			if (node[next_rank].area == imgsize)
			{
				if (node[current_rank].area == imgsize)
					next_rank = current_rank;
				node[next_rank].parentidx = next_rank;
				break;
			}

			iChild = current_rank;
			current_rank = next_rank;
		}

		curSize = 0;
		for (p = 0; p < nredges; p++)
		{
			if (node[p].parentidx > 0)
				curSize++;
		}

		q = nredges;
		for (p = 0; p < imgsize; p++)
		{
			if (node[parentAry[p]].level)//Singleton 0-CC
			{
				(&node[q])->set(1, 0, (double)img[p], img[p], img[p]);
				node[q++].parentidx = parentAry[p];
				parentAry[p] = q;
				curSize++;
			}
		}


//		delete hqueue;
		delete trie;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
//		Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}
	
public:
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Imgidx, Pixel>* node;
	Imgidx rootidx;
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
		if(sizeof(Pixel) == 1)
			Flood_HQueue(img);
		else
			Flood_Trie(img);
	}
};