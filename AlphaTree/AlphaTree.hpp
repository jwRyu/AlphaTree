#pragma once

#include "defines.h"
#include "allocator.h"
#include "NeighbourList.hpp"

#define NULL_LEVELROOT		0xffffffff
#define ANODE_CANDIDATE		0xfffffffe

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

#define A		1.10184
#define SIGMA	-2.6912
#define B		-0.0608
#define M		0.03
#define SIZE_ADD	0.05

class AlphaTree;
typedef uint32 (AlphaTree::*idxConvFn)(uint32, uint8);

typedef struct AlphaNode
{
	uint32 area;
	uint8 alpha;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	uint32 parentidx;
} AlphaNode;

class AlphaTree
{
	uint32 maxSize;
	uint32 curSize;
	uint32 height, width, channel;
	uint8 connectivity;
	AlphaNode* node;
	uint32* parentAry;
	int neighbours[8];
	idxConvFn img_2_dimg_idx;
	double nrmsd;

	void compute_dimg(uint8 * dimg, NeighbourList<4>& nlist, uint32 * dhist, Pixel * img);
	//void init_isVisited(uint8 *isVisited);
	void Flood(Pixel* img);
	inline void connectPix2Node(uint32* parentAry, uint32 pidx, Pixel pix_val, AlphaNode* pNode, uint32 iNode) const
	{	
		parentAry[pidx] = iNode;
		pNode->area++;
		pNode->maxPix = max(pNode->maxPix, pix_val);
		pNode->minPix = min(pNode->minPix, pix_val);
		pNode->sumPix += pix_val;
	}

	inline void connectNode2Node(AlphaNode* pPar, uint32 iPar, AlphaNode* pNode) const
	{
		pNode->parentidx = iPar;
		pPar->area += pNode->area;
		pPar->maxPix = max(pNode->maxPix, pPar->maxPix);
		pPar->minPix = min(pNode->minPix, pPar->minPix);
		pPar->sumPix += pNode->sumPix;
	}

	inline uint32 NewAlphaNode(uint8 level)
	{
		AlphaNode *pNew = this->node + this->curSize;

		if (this->curSize == this->maxSize)
		{
			printf("Reallocating...\n");
			this->maxSize = min(this->height * this->width, (uint32)(this->maxSize + (uint32)(this->height * this->width * SIZE_ADD)));

			this->node = (AlphaNode*)Realloc(this->node, this->maxSize * sizeof(AlphaNode));
			pNew = this->node + this->curSize;
		}
		pNew->alpha = level;
		pNew->minPix = (uint8)-1;
		pNew->minPix = 0;
		pNew->sumPix = 0.0;
		pNew->parentidx = 0;
		pNew->area = 0;

		return this->curSize++;
	}
	
	inline uint8 is_visited(uint8* isVisited, uint32 p) const
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}

	inline void visit(uint8* isVisited, uint32 p) const
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}


	inline uint32 img_2_dimg_idx_4N(uint32 pidx, uint8 neighbour)
	{
		uint8 b0 = (uint8)(pidx & 1), b1 = (uint8)((pidx >> 1) & 1);

		return ((pidx + b1 * neighbours[neighbour]) << 1) + (b0^b1);
	}

	inline uint32 img_2_dimg_idx_8N(uint32 pidx, uint8 neighbour)
	{
//		uint8 b0 = (uint8)(pidx & 1), b1 = (uint8)((pidx >> 1) & 1);

	//	return ((pidx + b1 * neighbours[neighbour]) << 1) + (b0^b1);
		return 0;
	}
	

public:
	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree() { Free(node); Free(parentAry); }

	void AlphaFilter(Pixel* outimg, double lambda, uint32 area);
	void BuildAlphaTree(Pixel *img, uint32 height, uint32 width, uint32 channel, uint8 connectivity);
	inline void clear() { Free(node); Free(parentAry); node = 0; parentAry = 0; }
};