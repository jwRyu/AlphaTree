#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "HQueue.hpp"
#include "allocator.h"
using namespace std;

#define OUTPUT_FNAME "C:/Users/jwryu/RUG/2018/AlphaTree/test.dat"

#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Grey"
#define INPUTIMAGE_DIR_COLOUR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Colour" //colour images are used after rgb2grey conversion
#define REPEAT 20
#define RUN_TSE_ONLY 0

#define DEBUG 0

#define DELAYED_ANODE_ALLOC		0
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define ANODE_CANDIDATE		0xfffffffe

#define dimg_idx_v(pidx) ((pidx)<<1)
#define dimg_idx_h(pidx) ((pidx)<<1)+1

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

#define A		1.10184
#define SIGMA	-2.6912
#define B		-0.0608
#define M		0.03

//Memory allocation reallocation schemes
#define TSE 0
#define MAXIMUM 1
#define LINEAR 2
#define EXP 3
int mem_scheme = -1;
double size_init[4] = { -1, 1, .2, .15 };
double size_mul[4] = { 1, 1, 1, 2 };
double size_add[4] = { .05, 0, 0.15, 0 };

double nrmsd;

#if DEBUG
void* buf;
uint64 bufsize;
void save_buf(void* src, uint64 size)
{
	memcpy(buf, src, size);
	bufsize = size;
}
uint8 isChanged(void *src)
{
	uint64 i;
	for (i = 0; i < bufsize; i++)
	{
		if (((uint8*)buf)[i] != ((uint8*)src)[i])
			return 1;
	}
	return 0;
}
#endif


typedef struct AlphaNode
{
	uint32 area;
	uint8 level;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	uint32 parentidx;
} AlphaNode;

class AlphaTree
{
public:
	uint32 maxSize;
	uint32 curSize;
	uint32 height, width, channel;
	AlphaNode* node;
	uint32* parentAry;

	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree() { Free(node); Free(parentAry); }
	inline void clear() { Free(node); Free(parentAry); curSize = 0; }

	inline void connectPix2Node(uint32 pidx, Pixel pix_val, uint32 iNode)
	{
		AlphaNode *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->area++;
		pNode->maxPix = max(pNode->maxPix, pix_val);
		pNode->minPix = min(pNode->minPix, pix_val);
		pNode->sumPix += pix_val;
	}

	inline void connectNode2Node(AlphaNode* pPar, uint32 iPar, AlphaNode* pNode)
	{
		pNode->parentidx = iPar;
		pPar->area += pNode->area;
		pPar->maxPix = max(pNode->maxPix, pPar->maxPix);
		pPar->minPix = min(pNode->minPix, pPar->minPix);
		pPar->sumPix += pNode->sumPix;
	}

	inline uint32 NewAlphaNode(uint8 level)
	{
		AlphaNode *pNew = node + curSize;

		if (curSize == maxSize)
		{
			printf("Reallocating...\n");
			maxSize = min(height * width, (uint32)(size_mul[mem_scheme] * maxSize) + (uint32)(height * width * size_add[mem_scheme]));

			node = (AlphaNode*)Realloc(node, maxSize * sizeof(AlphaNode));
			pNew = node + curSize;
		}
		pNew->level = level;
		pNew->minPix = (uint8)-1;
		pNew->minPix = 0;
		pNew->sumPix = 0.0;
		pNew->parentidx = 0;
		pNew->area = 0;

		return curSize++;
	}

	inline uint8 is_visited(uint8* isVisited, uint32 p)
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}

	inline void visit(uint8* isVisited, uint32 p)
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

	void compute_dimg(uint8* dimg, uint32* dhist, Pixel* img);

	void Flood(Pixel* img);
	void BuildAlphaTree(Pixel *img, uint32 height, uint32 width, uint32 channel);
};

void AlphaTree::compute_dimg(uint8* dimg, uint32* dhist, Pixel* img)
{
	uint32 dimgidx, imgidx, stride_w = width, i, j;

	imgidx = dimgidx = 0;
	for (i = 0; i < height - 1; i++)
	{
		for (j = 0; j < width - 1; j++)
		{
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + stride_w] - (int)img[imgidx]));
			dhist[dimg[dimgidx++]]++;
			dimg[dimgidx] = (uint8)(abs((int)img[imgidx + 1] - (int)img[imgidx]));
			dhist[dimg[dimgidx++]]++;
			imgidx++;
		}
		dimg[dimgidx] = (uint8)(abs((int)img[imgidx + stride_w] - (int)img[imgidx]));
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
}



void AlphaTree::Flood(Pixel* img)
{
	uint32 imgsize, dimgsize, nredges, max_level, current_level, next_level, x0, p, dissim;
	uint32 numlevels;
	HQueue<uint32>* hqueue;
	uint32 *dhist;
	uint8 *dimg;
	uint32 iChild, *levelroot;
	uint8 *isVisited;
	uint32 *pParentAry;

	imgsize = width * height;
	nredges = width * (height - 1) + (width - 1) * height;
	dimgsize = 2 * width * height; //To make indexing easier
	numlevels = 1 << (8 * sizeof(uint8));

	dhist = (uint32*)Malloc((size_t)numlevels * sizeof(uint32));
	dimg = (uint8*)Malloc((size_t)dimgsize * sizeof(uint8));
	levelroot = (uint32*)Malloc((uint32)(numlevels + 1) * sizeof(uint32));
	isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
	for (p = 0; p < numlevels; p++)
		levelroot[p] = NULL_LEVELROOT;
	memset(dhist, 0, (size_t)numlevels * sizeof(uint32));
	memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));

	max_level = (uint8)(numlevels - 1);

	compute_dimg(dimg, dhist, img);
	dhist[max_level]++;
	hqueue = new HQueue<uint32>(nredges + 1, dhist, numlevels);
	
	//tree size estimation (TSE)
	nrmsd = 0;
	for (p = 0; p < numlevels; p++)
		nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	if (mem_scheme == TSE)
		maxSize = min(imgsize, (uint32)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
	else
		maxSize = (uint32)(imgsize * size_init[mem_scheme]);

	Free(dhist);

	parentAry = (uint32*)Malloc((size_t)imgsize * sizeof(uint32));
	node = (AlphaNode*)Malloc((size_t)maxSize * sizeof(AlphaNode));
	pParentAry = parentAry;

	levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];

	current_level = max_level;
	x0 = imgsize >> 1;
	hqueue->hqueue_push(x0, current_level);

	iChild = levelroot[max_level + 1];
	while (current_level <= max_level)
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

			if (LEFT_AVAIL(p, width) && !is_visited(isVisited, p - 1))
			{
				dissim = (uint32)dimg[dimg_idx_h(p - 1)];
				hqueue->hqueue_push(p - 1, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
					levelroot[dissim] = ANODE_CANDIDATE;
#else
					levelroot[dissim] = NewAlphaNode((uint8)dissim);
#endif
			}
			if (RIGHT_AVAIL(p, width) && !is_visited(isVisited, p + 1))
			{
				dissim = (uint32)dimg[dimg_idx_h(p)];
				hqueue->hqueue_push(p + 1, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
					levelroot[dissim] = ANODE_CANDIDATE;
#else
					levelroot[dissim] = NewAlphaNode((uint8)dissim);
#endif
			}
			if (UP_AVAIL(p, width) && !is_visited(isVisited, p - width))
			{
				dissim = (uint32)dimg[dimg_idx_v(p - width)];
				hqueue->hqueue_push(p - width, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
					levelroot[dissim] = ANODE_CANDIDATE;
#else
					levelroot[dissim] = NewAlphaNode((uint8)dissim);
#endif
			}
			if (DOWN_AVAIL(p, width, imgsize) && !is_visited(isVisited, p + width))
			{
				dissim = (uint32)dimg[dimg_idx_v(p)];
				hqueue->hqueue_push(p + width, dissim);
				if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_ANODE_ALLOC
					levelroot[dissim] = ANODE_CANDIDATE;
#else
					levelroot[dissim] = NewAlphaNode((uint8)dissim);
#endif
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
			curSize--;

			memset((uint8*)(node + curSize), 0, sizeof(AlphaNode));
		}

		next_level = current_level + 1;
		while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
			next_level++;
#if DELAYED_ANODE_ALLOC
		if (levelroot[next_level] == ANODE_CANDIDATE)
			levelroot[next_level] = NewAlphaNode((uint8)next_level);
#endif
		connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);

		iChild = levelroot[current_level];
		levelroot[current_level] = NULL_LEVELROOT;
		current_level = next_level;

	}
	delete hqueue;
	Free(dimg);
	Free(levelroot);
	Free(isVisited);
}


void AlphaTree::BuildAlphaTree(Pixel *img, uint32 height, uint32 width, uint32 channel)
{
	this->height = height;
	this->width = width;
	this->channel = channel;
	curSize = 0;
	Flood(img);
}
/*
void DeleteAlphaTree(AlphaTree* tree)
{
	Free(tree->parentAry);
	Free(tree->node);
	Free(tree);
}*/

int main(int argc, char **argv)
{
	AlphaTree *tree;
	uint32 width, height, channel;
	uint32 cnt = 0;
	ofstream f;
	ifstream fcheck;
	char in;
	uint32 i, contidx;
	std::string path;
	double time_elapsed = 0;
	double pixels_processed = 0;


	contidx = 0;
	//	f.open("C:/Users/jwryu/RUG/2018/AlphaTree/AlphaTree_grey_Exp.dat", std::ofstream::app);
	fcheck.open(OUTPUT_FNAME);
	if (fcheck.good())
	{
		cout << "Output file \"" << OUTPUT_FNAME << "\" already exists. Overwrite? (y/n/a)";
		//cin >> in;
		in = 'y';
		if (in == 'a')
		{
			f.open(OUTPUT_FNAME, std::ofstream::app);
			cout << "Start from : ";
			cin >> contidx;
		}
		else if (in == 'y')
			f.open(OUTPUT_FNAME);
		else
			exit(-1);
	}
	else
		f.open(OUTPUT_FNAME);

	cnt = 0;
	for (mem_scheme = 0; mem_scheme < 4; mem_scheme++) // memory scheme loop (TSE, Max, Linear, Exp)
	{
#if RUN_TSE_ONLY
		if (mem_scheme > 0)
			break;
#endif
		for (i = 0; i < 2; i++) // grey, colour loop
		{
			if (i == 0)
				path = INPUTIMAGE_DIR;
			else
				path = INPUTIMAGE_DIR_COLOUR;

			for (auto & p : std::experimental::filesystem::directory_iterator(path))
			{
				if (++cnt < contidx)
				{
					cout << cnt << ": " << p << endl;
					continue;
				}
				cv::String str1(p.path().string().c_str());
				cv::Mat cvimg;
				if (i == 0)
					cvimg = imread(str1, cv::IMREAD_GRAYSCALE);
				else
				{
					cvimg = imread(str1, cv::IMREAD_COLOR);
					cv::cvtColor(cvimg, cvimg, CV_BGR2GRAY);
				}

				/*
				cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE);// Create a window for display.
				cv::imshow("Display window", cvimg);                   // Show our image inside it.
				cv::waitKey(0);
				getc(stdin);
				*/

				height = cvimg.rows;
				width = cvimg.cols;
				channel = cvimg.channels();

				cout << cnt << ": " << str1 << ' ' << height << 'x' << width << endl;

				if (channel != 1)
				{
					cout << "input should be a greyscale image" << endl;
					getc(stdin);
					exit(-1);
				}

				double runtime, minruntime;
				for (int testrep = 0; testrep < REPEAT; testrep++)
				{
					memuse = max_memuse = 0;
					auto wcts = std::chrono::system_clock::now();

					tree = (AlphaTree*)Malloc(sizeof(AlphaTree));
					//		start = clock();
					tree->BuildAlphaTree(cvimg.data, height, width, channel);

					std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
					runtime = wctduration.count();
					minruntime = testrep == 0 ? runtime : min(runtime, minruntime);

					if (testrep < (REPEAT - 1))
						tree->clear();
				}
				f << p.path().string().c_str() << '\t' << height << '\t' << width << '\t' << max_memuse << '\t' << nrmsd << '\t' << tree->maxSize << '\t' << tree->curSize << '\t' << minruntime << mem_scheme << i << endl;

				pixels_processed += width * height;
				time_elapsed += minruntime;
				cout << "Time Elapsed: " << minruntime << " mean processing speed(Mpix/s): " << pixels_processed / (time_elapsed * 1000000) << endl;
				cvimg.release();
				str1.clear();
				tree->clear();
				//return 0;
			}
		}
	}

	f.close();
	return 0;
}