#include <cstdlib>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <cfloat>
#include <cstring>
#include "defines.h"
#include "AlphaTree.hpp"
#include "allocator.h"
#include "HQueue.hpp"

#define img2iVidx(idx, width) (((idx)/(width)) + (width) + 3)

void AlphaTree::AlphaFilter(Pixel* outimg, double lambda, uint32 area)
{
	uint32 i;
	Pixel *outval = new Pixel[curSize];
	uint8 *filtered = new uint8[curSize];
	
	for (i = 0; i < curSize; i++)
	{
		if (node[i].alpha >= lambda && node[i].area > area)
		{
			filtered[i] = 1;
			outval[i] = (Pixel)(node[i].sumPix / (double)node[i].area);
//			outval[i] = (Pixel)((double)node[i].area / (double)(width * height) * (255.0));
		}
		else
			filtered[i] = 0;
	}
	if (lambda >= node[0].alpha)
		return;

	for (i = 0; i < height * width; i++)
	{
		uint32 j = parentAry[i];
		while (!filtered[j])
			j = node[j].parentidx;
		outimg[i] = outval[j];
	}
	delete outval;
	delete filtered;
}

void AlphaTree::compute_dimg(uint8* dimg, uint32* dhist, Pixel* img)
{
	uint32 dimgidx, imgidx, i, j;
	HQueue<neighidx> *hqueue;

	if (connectivity == 4)
	{
		//   -  3  -
		//   2  p  1
		//   -  0  -
		imgidx = dimgidx = 0;
		for (i = 0; i < height - 1; i++)
		{
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width] - (int)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width] - (int)img[imgidx]));
			dhist[dimg[dimgidx++]]++;
			dimgidx++;
			imgidx++;
		}
		for (j = 0; j < width - 1; j++)
		{
			dimgidx++;
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));
			dhist[dimg[dimgidx++]]++;
			imgidx++;
		}
		img += width * height;

		hqueue_new(&hqueue, height * (width - 1) + (height - 1) * width, dhist, 256, 4);

		imgidx = dimgidx = 0;
		for (i = 0; i < height - 1; i++)
		{
			for (j = 0; j < width - 1; j++)
			{
				hqueue_push(hqueue, imgidx, 0, dimg[dimgidx]);
				hqueue_push(hqueue, imgidx + width, 2, dimg[dimgidx++]);
				hqueue_push(hqueue, imgidx, 1, dimg[dimgidx]);
				hqueue_push(hqueue, imgidx + 1, 3, dimg[dimgidx++]);
				imgidx++;
			}
			hqueue_push(hqueue, imgidx, 0, dimg[dimgidx]);
			hqueue_push(hqueue, imgidx + width, 2, dimg[dimgidx++]);
			dimgidx++;
			imgidx++;
		}
		for (j = 0; j < width - 1; j++)
		{
			dimgidx++;
			hqueue_push(hqueue, imgidx, 1, dimg[dimgidx]);
			hqueue_push(hqueue, imgidx + 1, 3, dimg[dimgidx++]);
			imgidx++;
		}
		img += width * height;
	}
	else if (connectivity == 8)
	{
		//   -  -  -
		//   -  p  4
		//   1  2  3
		imgidx = dimgidx = 0;
		for (i = 0; i < height - 1; i++)
		{
			dimgidx++; //skip 1
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width] - (int)img[imgidx]));//2
			dhist[dimg[dimgidx++]]++;
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width + 1] - (int)img[imgidx]));//3
			dhist[dimg[dimgidx++]]++;
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//4
			dhist[dimg[dimgidx++]]++;
			imgidx++;
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width - 1] - (int)img[imgidx]));//1
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width] - (int)img[imgidx]));//2
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width + 1] - (int)img[imgidx]));//3
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//4
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width - 1] - (int)img[imgidx]));//1
			dhist[dimg[dimgidx++]]++;
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + width] - (int)img[imgidx]));//2
			dhist[dimg[dimgidx++]]++;
			dimgidx += 2;//skip 3,4
			imgidx++;
		}
		for (j = 0; j < width - 1; j++)
		{
			dimgidx += 3; //skip 1,2,3
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));//4
			dhist[dimg[dimgidx++]]++;
			imgidx++;
		}
		img += width * height;

	}
}

void AlphaTree::init_isVisited(uint8 * isVisited)
{
	uint64 *p = (uint64*)isVisited;
	uint32 i, j;

	memset(isVisited, 0, ((width + 2)*(height + 2) + 7) >> 3);
	//First row
	for (i = 0; i < (width + 2) >> 6; i++)
	{
		p[i] = 0xffffffffffffffff;
	}
	for (i = (width + 2) & ~(0x3f); i < width + 2; i++)
		visit(isVisited, i);

	//Middle rows
	j = width + 2;
	for (i = 0; i < height; i++)
	{
		visit(isVisited, j);
		visit(isVisited, j + width + 1);
		j += width + 2;
	}

	//Last row
	i = (width + 2) * (height + 1);
	while (i & 0x3f)
		visit(isVisited, i++);
	for (; i < (width + 2) * (height + 2) - 64; i++)
	{
		*(uint64*)(&isVisited[i>>3]) = 0xffffffffffffffff;
		i += 64;
	}
	for (; i < (width + 2)*(height + 2); i++)
		visit(isVisited, i++);
}

void AlphaTree::Flood(Pixel* img)
{
	uint32 imgsize, dimgsize, nredges, max_level, current_level, next_level, x0, p, q, dissim;
	uint32 numlevels;
	HQueue<uint32>* hqueue;
	uint32 *dhist;
	uint8 *dimg;
	uint32 iChild, *levelroot;
	uint8 *isVisited;
	uint32 *pParentAry;
	uint8 iNeighbour;

	height = 60;
	width = 120;

	imgsize = width * height;
	nredges = width * (height - 1) + (width - 1) * height + (connectivity == 8) ? (width - 1) * (height - 1) * 2 : 0;
	dimgsize = (connectivity >> 1) * width * height; //To make indexing easier
	numlevels = 1 << (8 * sizeof(uint8));

	//Temp memory allocation
	dhist = (uint32*)Malloc((size_t)numlevels * sizeof(uint32));
	dimg = (uint8*)Malloc((size_t)dimgsize * sizeof(uint8));
	levelroot = (uint32*)Malloc((uint32)(numlevels + 1) * sizeof(uint32));
	isVisited = (uint8*)Malloc((size_t)((width + 2) * (height + 2) + 7) >> 3);
	for (p = 0; p < numlevels; p++)
		levelroot[p] = NULL_LEVELROOT;
	memset(dhist, 0, (size_t)numlevels * sizeof(uint32));
	init_isVisited(isVisited);
	
	max_level = (uint8)(numlevels - 1);

	//Calc. dissim, Histogram
	compute_dimg(dimg, dhist, img);
	dhist[max_level]++; //make a room for the root
	hqueue_new(&hqueue, nredges + 1, dhist, numlevels);
	
	//tree size estimation (TSE)
	nrmsd = 0;
	for (p = 0; p < numlevels; p++)
		nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	this->maxSize = min(imgsize, (uint32)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));

	Free(dhist);

	//Allocate parent array and nodes
	this->parentAry = (uint32*)Malloc((size_t)imgsize * sizeof(uint32));
	this->node = (AlphaNode*)Malloc((size_t)this->maxSize * sizeof(AlphaNode));
	pParentAry = this->parentAry;

	//put root node in levelroot
	levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	this->node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];

	current_level = max_level;
	x0 = imgsize >> 1; // First pixel (anything between 0 ~ imgsize-1)
	hqueue_push(hqueue, x0, current_level);

	//	free(dhist);

	iChild = levelroot[max_level + 1];
	while (current_level <= max_level)
	{
		while (hqueue->min_level <= current_level)
		{
			p = hqueue_pop(hqueue);
			if (is_visited(isVisited, img2iVidx(p, width)))
			{
				hqueue_find_min_level(hqueue);
				continue;
			}
			visit(isVisited, img2iVidx(p, width));

			for (iNeighbour = 0; iNeighbour < connectivity; iNeighbour++)
			{
				q = p + neighbours[iNeighbour];
				if (!is_visited(isVisited, img2iVidx(q, width)))
				{
					dissim = (uint32)dimg[dimg_idx_h(p - 1)];
					hqueue_push(hqueue, p - 1, dissim);
					if (levelroot[dissim] == NULL_LEVELROOT)
						levelroot[dissim] = ANODE_CANDIDATE;
				}
			}

			if (LEFT_AVAIL(p, width) && !is_visited(isVisited, p - 1))
			{
				dissim = (uint32)dimg[dimg_idx_h(p - 1)];
				hqueue_push(hqueue, p - 1, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
					levelroot[dissim] = ANODE_CANDIDATE;
			}
			if (RIGHT_AVAIL(p, width) && !is_visited(isVisited, p + 1))
			{
				dissim = (uint32)dimg[dimg_idx_h(p)];
				hqueue_push(hqueue, p + 1, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
					levelroot[dissim] = ANODE_CANDIDATE;
			}
			if (UP_AVAIL(p, width) && !is_visited(isVisited, p - width))
			{
				dissim = (uint32)dimg[dimg_idx_v(p - width)];
				hqueue_push(hqueue, p - width, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
					levelroot[dissim] = ANODE_CANDIDATE;
			}
			if (DOWN_AVAIL(p, width, imgsize) && !is_visited(isVisited, p + width))
			{
				dissim = (uint32)dimg[dimg_idx_v(p)];
				hqueue_push(hqueue, p + width, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
					levelroot[dissim] = ANODE_CANDIDATE;
			}

			if (current_level > hqueue->min_level)
				current_level = hqueue->min_level;
			hqueue_find_min_level(hqueue);

			if (levelroot[current_level] == ANODE_CANDIDATE)
				levelroot[current_level] = NewAlphaNode((uint8)current_level);
			connectPix2Node(pParentAry, p, img[p], this->node + levelroot[current_level], levelroot[current_level]);

		}
		//		if(this->curSize > 22051838 && (this->curSize))
			//		printf("curSize: %d\n",this->curSize);
				//Redundant node removal
		if (this->node[iChild].parentidx == levelroot[current_level] &&
			this->node[levelroot[current_level]].area == this->node[iChild].area)
		{
			levelroot[current_level] = iChild;
			this->curSize--;

			memset((uint8*)(this->node + this->curSize), 0, sizeof(AlphaNode));
		}

		next_level = current_level + 1;
		while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
			next_level++;
		if (levelroot[next_level] == ANODE_CANDIDATE)
			levelroot[next_level] = NewAlphaNode((uint8)next_level);
		connectNode2Node(this->node + levelroot[next_level], levelroot[next_level], this->node + levelroot[current_level]);

		iChild = levelroot[current_level];
		levelroot[current_level] = NULL_LEVELROOT;
		current_level = next_level;

	}

	hqueue_free(hqueue);
	Free(dimg);
	Free(levelroot);
	Free(isVisited);
}


void AlphaTree::BuildAlphaTree(Pixel *img, uint32 height, uint32 width, uint32 channel, uint8 connectivity)
{
	this->height = height;
	this->width = width;
	this->channel = channel;
	this->connectivity = connectivity;
	if (connectivity == 4)
	{
		neighbours[0] = -(int)width;
		neighbours[1] = +1;
		neighbours[2] = -1;
		neighbours[3] = width;

	}
	else if (connectivity == 8)
	{
		neighbours[0] = -(int)width - 1;
		neighbours[1] = -(int)width;
		neighbours[2] = -(int)width + 1;
		neighbours[3] = +1;
		neighbours[4] = -1;
		neighbours[5] = width - 1;
		neighbours[6] = width;
		neighbours[7] = width + 1;
	}

	Flood(img);
}