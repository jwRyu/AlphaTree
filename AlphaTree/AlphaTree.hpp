#pragma once

#include "defines.h"
#include "allocator.h"

#define PADDED_VISIT_ARRAY	1

#define NULL_LEVELROOT		0xffffffff
#define ANODE_CANDIDATE		0xfffffffe

//#define dimg_idx_v(pidx) ((pidx)<<1)
//#define dimg_idx_h(pidx) ((pidx)<<1)+1

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
typedef uint32 (AlphaTree::*idxConv)(uint32&, uint8&) const;

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
	int64 maxSize;
	int64 curSize;
	int64 height, width, channel;
	uint8 connectivity;
	AlphaNode* node;
	uint32* parentAry; //replace later
	int64 neighbours[8];
	idxConv imgidx_to_dimgidx;
	double nrmsd;

	void compute_dimg(uint8 * dimg, uint32 * dhist, Pixel * img);
	void init_isVisited(uint8 *isVisited, size_t arrsize);
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

	inline uint32 imgidx_to_dimgidx_4N(uint32& idx, uint8& neighbour) const
	{
		uint8 p0 = neighbour & 1;
		uint8 p1 = (neighbour >> 1) & 1;
		return ((idx - (uint32)p1 - (uint32)p1 * (uint32)p0 * ((uint32)width - 1)) << 1) + ((uint32)p0 ^ (uint32)p1);
	}
	
	inline int64 imgidx_to_dimgidx_4N(int64& idx, uint8& neighbour) const
	{
		uint8 p0 = neighbour & 1;
		uint8 p1 = (neighbour >> 1) & 1;
		return ((idx - (int64)p1 - (int64)p1 * (int64)p0 * ((int64)width - 1)) << 1) + ((int64)p0^(int64)p1);
	}

	inline uint32 imgidx_to_dimgidx_8N(uint32& idx, uint8& neighbour) const
	{
//		uint8 p0 = neighbour & 1;
//		uint8 p1 = (neighbour >> 1) & 1;
//		return ((idx - p1 - p1 * p0 * (width - 1)) << 1) + (p0^p1);
		return 0;
	}

	inline uint8 is_visited(uint8* isVisited, uint32 p) const
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}
	
	inline uint8 is_visited(uint8* isVisited, int64 p) const
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}

	inline void visit(uint8* isVisited, uint32 p) const
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

	inline void visit(uint8* isVisited, int64 p) const
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

public:
	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree() { Free(node); Free(parentAry); }

	void AlphaFilter(Pixel* outimg, double lambda, uint32 area);
	void BuildAlphaTree(Pixel *img, uint32 height, uint32 width, uint32 channel, uint8 connectivity);
	inline void clear() { Free(node); Free(parentAry); node = 0; parentAry = 0; }
};
  
template <typename T>
inline T img2Vidx(T idx, int64 width)
{
	return (T)((int64)idx + (int64)idx / width + width + 2);
}

