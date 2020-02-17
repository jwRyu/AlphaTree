#pragma once
#pragma once
#include<iostream>
#include<math.h>
#include "defines.h"
#include "HierarQueue.h"
#include "Trie.h"
#include "HybridQueue.h"
#include "HeapQueue.h"

//using namespace std;
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


#define	IMGIDX_32BITS		0
#define	IMGIDX_64BITS		1
#define	PIXEL_8BIT			0
#define	PIXEL_16BIT			1
#define	PIXEL_32BIT			2
#define	PIXEL_64BIT			3
#define	PIXEL_FLOAT			4
#define	PIXEL_DOUBLE		5

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

//tmptmptmptmptmptmtmp
int disp = 0;

template<class Imgidx, class Pixel>
class AlphaNode
{
public:
	Imgidx area;
	Pixel alpha;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	Imgidx parentidx;

	Imgidx rootidx;
	Pixel filter_val;

	inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	{
		this->area = area_in;
		this->alpha = level;
		this->sumPix = sumPix_in;
		this->minPix = minPix_in;
		this->maxPix = maxPix_in;
		this->parentidx = -1;
	}

	inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in, int8 pix_type)
	{
		this->area = area_in;
		this->alpha = level;
		if (pix_type < 4)
			this->sumPix = (double)sumPix_in;
		else if (pix_type == PIXEL_DOUBLE)
			this->sumPix = *((double*)&sumPix_in);
		else
			this->sumPix = *((float*)&sumPix_in);
		this->minPix = minPix_in;
		this->maxPix = maxPix_in;
		this->parentidx = -1;
	}

	inline void add(AlphaNode* q)//, int8 pix_type)
	{
		if (disp) printf("Entering add_this\n");
		if (disp) this->print(this);
		this->area += q->area;
		this->sumPix += q->sumPix;

		this->minPix = min(this->minPix, q->minPix);
		this->maxPix = max(this->maxPix, q->maxPix);
		if (disp) printf("Exiting add_this\n");
		if (disp) this->print(this);
		if (disp) getchar();
	}
	inline void add(Pixel pix_val, int8 pix_type)
	{
		this->area++;
		if (pix_type < 4)
			this->sumPix += (double)pix_val;
		else if (pix_type == PIXEL_DOUBLE)
			this->sumPix += *((double*)&pix_val);
		else
			this->sumPix += *((float*)&pix_val);
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
	inline void connect_to_parent(AlphaNode* pPar, Imgidx iPar)
	{
		if (disp) printf("conn2par_enter\n");
		this->parentidx = iPar;
		pPar->add(this);
	}

	void print(AlphaNode* node)
	{
		printf("Node idx: %lld\nparent: %lld\narea: %lld\nsumpix: %lf\nminpix: %lld\nmaxpix: %lld\n\n"
			, (uint64)(this - node)
			, (uint64)this->parentidx
			, (uint64)this->area
			, (double)this->sumPix
			, (uint64)this->minPix
			, (uint64)this->maxPix);
	}
};

template<class Imgidx, class Pixel>
class RankItem
{
public:
	Pixel alpha;
	Imgidx dimgidx;
	Imgidx p, q;

	inline void operator=(const RankItem& item)
	{
		this->alpha = item.alpha;
		this->dimgidx = item.dimgidx;
		this->p = item.p;
		this->q = item.q;
	}
};

struct pixdata
{
	uint8 isVisited;
	uint8 isAvailable;
	uint8 edgerank;
	uint8 cnt;
};

template<class Imgidx, class Pixel>
class ATree
{
	inline Pixel abs_diff(Pixel p, Pixel q)
	{
		if (p > q)
			return p - q;
		else
			return q - p;
	}

	inline uint8 compute_incidedge_queue(Pixel d0, Pixel d1)
	{
		if (d0 <= d1)
			return 0x4;
		else
			return 0x1;
	}
	inline uint8 compute_incidedge_queue(Pixel d0, Pixel d1, Pixel d2)
	{
		if (d0 <= d1)
		{
			if (d0 < d2)
			{
				if (d1 < d2) //00 10 01 00
					return 0x24;
				else //00 01 10 00
					return 0x18;
			}
			else //00 01 00 10
				return 0x12;
		}
		else
		{
			if (d1 < d2)
			{
				if (d0 < d2) //00 10 00 01
					return 0x21;
				else //00 00 10 01
					return 0x09;
			}
			else //00 00 01 10
				return 0x06;
		}
	}
	inline uint8 compute_incidedge_queue(Pixel d0, Pixel d1, Pixel d2, Pixel d3)
	{ //when dx == dy, d3 < d2 < d0 < d1
		if (d0 <= d1)
		{
			if (d0 < d2)
			{
				if (d0 < d3)
				{
					if (d1 < d2)
					{
						if (d1 < d3)
						{
							if (d2 < d3) //11 10 01 00
								return 0xe4;
							else //10 11 01 00
								return 0xb4;
						}
						else // 10 01 11 00
							return 0x9c;
					}
					else
					{
						if (d3 <= d1)
						{
							if (d3 <= d2) //01 10 11 00
								return 0x6c;
							else //01 11 10 00
								return 0x78;
						}
						else //11 01 10 00!
							return 0xd8;
					}
				}
				else
				{
					if (d1 < d2)//10 01 00 11
						return 0x93;
					else//01 10 00 11
						return 0x63;
				}
			}
			else
			{
				if (d3 <= d2)     //01 00 10 11
					return 0x4b;
				else if (d3 <= d0)//01 00 11 10
					return 0x4e;
				else if (d3 <= d1)//01 11 00 10
					return 0x72;
				else             //11 01 00 10
					return 0xd2;
			}
		}
		else
		{
			if (d1 < d2)
			{
				if (d1 < d3)
				{
					if (d0 < d2)
					{
						if (d0 < d3)
						{
							if (d2 < d3) //11 10 00 01
								return 0xe1;
							else //10 11 00 01
								return 0xb1;
						}
						else // 10 00 11 01
							return 0x8d;
					}
					else
					{
						if (d3 <= d0)
						{
							if (d3 <= d2) //00 10 11 01
								return 0x2d;
							else //00 11 10 01
								return 0x39;
						}
						else //11 00 10 01
							return 0xc9;
					}
				}
				else
				{
					if (d0 < d2)//10 00 01 11
						return 0x87;
					else//00 10 01 11
						return 0x27;
				}
			}
			else
			{
				if (d3 <= d2)     //00 01 10 11
					return 0x1b;
				else if (d3 <= d1)//00 01 11 10
					return 0x1e;
				else if (d3 <= d0)//00 11 01 10
					return 0x36;
				else             //11 00 01 10
					return 0xc6;
			}
		}
	}

#define _COMPUTE_RANK_PIXEL \
	rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))\
			: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));\
	rankinfo[contidx].p = imgidx;\
	rankinfo[contidx].q = imgidx + width;\
	rankinfo[contidx++].dimgidx = dimgidx++;\
	rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))\
			: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));\
	rankinfo[contidx].p = imgidx;\
	rankinfo[contidx].q = imgidx + 1;\
	rankinfo[contidx++].dimgidx = dimgidx++;

#define _COMPUTE_RANK_PIXEL_LASTCOL \
	rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))\
		: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));\
	rankinfo[contidx].p = imgidx;\
	rankinfo[contidx].q = imgidx + width;\
	rankinfo[contidx++].dimgidx = dimgidx;\
	dimgidx += 2;

#define _COMPUTE_RANK_PIXEL_LASTROW \
	dimgidx++;\
	rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))\
		: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));\
	rankinfo[contidx].p = imgidx;\
	rankinfo[contidx].q = imgidx + 1;\
	rankinfo[contidx++].dimgidx = dimgidx++;

	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>*& rankinfo, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;
		int8 p64 = (pix_type == PIXEL_64BIT) || (pix_type == PIXEL_DOUBLE);
		Pixel *incidents = (Pixel*)Malloc(connectivity * width * height * sizeof(Pixel));

		if (connectivity == 4)
			nredges = (width - 1) * height + width * (height - 1);
		else
			nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);



		//rankinfo = (RankItem1<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem1<Imgidx, Pixel>));
		tmp = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc(max(hist_size * sizeof(Imgidx) * sizeof(Pixel), (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1));

		//int aa = max(hist_size, (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1);

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					_COMPUTE_RANK_PIXEL
						imgidx++;
				}
				_COMPUTE_RANK_PIXEL_LASTCOL
					imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				_COMPUTE_RANK_PIXEL_LASTROW
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
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width - 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width - 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width + 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width - 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width - 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
			}
		}

		for (i = 0; i < (Imgidx)(sizeof(Pixel) * 65536 >> 1); i++)
			hist[i] = 0;

		if (bit_depth < 16) // count sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankinfo[i].alpha;
				hist[hidx]++;
			}

			h = hist;

			hsum = 0;
			for (i = 0; i < (Imgidx)(hist_size); i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = rankinfo[i].alpha;
				j = --h[hidx];
				tmp[j] = rankinfo[i]; //slow!
			}
			r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
		}
		else //Radix sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankinfo[i].alpha;
				hist[hidx & mask]++;

				shamt = (Pixel)16;
				h_offset = (Pixel)65536;
				for (nbyte = 2; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
				{
					hidx = hidx >> shamt;
					hist[h_offset + (hidx & mask)]++;
					h_offset += (Pixel)65536;
				}
			}

			h = hist;
			shamt = 0;
			for (nbyte = 0; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
			{
				hsum = 0;
				for (i = 0; i < (Imgidx)(hist_size); i++)
				{
					hsum += h[i];
					h[i] = hsum;
				}
				for (i = nredges - 1; i >= 0; i--)
				{
					hidx = (rankinfo[i].alpha >> shamt) & mask;
					j = --h[hidx];
					tmp[j] = rankinfo[i]; //slow!
				}
				r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
				shamt += 16;
				h += hist_size;
			}
		}
		for (i = 0; i < nredges; i++)
		{
			//rank_to_alpha[i] = rankinfo[i].alpha;
			rank[rankinfo[i].dimgidx] = i;
		}



		//std::cout << "in_func2: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		//Free(rankinfo);
		Free(tmp);
		Free(hist);
		Free(incidents);
	}

	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>*& rankinfo, pixdata *pd, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges,wstride;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;
		int8 p64 = (pix_type == PIXEL_64BIT) || (pix_type == PIXEL_DOUBLE);
		Pixel *incidents = (Pixel*)Malloc(connectivity * width * height * sizeof(Pixel));
		Pixel maxpixval = (Pixel)-1;
		wstride = width * 2 - 1;
		if (connectivity == 4)
			nredges = (width - 1) * height + width * (height - 1);
		else
			nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);

		// 		uint8 tmptmp[4],ret;
// 		for (int i0 = 0; i0 < 4; i0++)
// 		{
// 			tmptmp[0] = i0;
// 			for (int i1 = 0; i1 < 4; i1++)
// 			{
// 				if (tmptmp[0] == i1)
// 					continue;
// 				tmptmp[1] = i1;
// 				for (int i2 = 0; i2 < 4; i2++)
// 				{
// 					if (tmptmp[0] == i2 || tmptmp[1] == i2)
// 						continue;
// 					tmptmp[2] = i2;
// 					ret = compute_incidedge_queue(tmptmp[0], tmptmp[1], tmptmp[2]);
// 					for (int shamt = 0; shamt < 4; shamt += 2)
// 					{
// 						if (tmptmp[(ret >> shamt) & 3] > tmptmp[(ret >> (shamt+2)) & 3])
// 							compute_incidedge_queue(tmptmp[0], tmptmp[1], tmptmp[2]);;
// 					}
// 					for (int i3 = 0; i3 < 4; i3++)
// 					{
// 						if (tmptmp[0] == i3 || tmptmp[1] == i3 || tmptmp[2] == i3)
// 							continue;
// 						tmptmp[3] = i3;
// 						ret = compute_incidedge_queue(tmptmp[0], tmptmp[1], tmptmp[2], tmptmp[3]);
// 						for (int shamt = 0; shamt < 6; shamt += 2)
// 						{
// 							if (tmptmp[(ret >> shamt) & 3] > tmptmp[(ret >> (shamt + 2)) & 3])
// 								compute_incidedge_queue(tmptmp[0], tmptmp[1], tmptmp[2], tmptmp[3]);
// 						}
// 					}
// 				}
// 			}
// 		}

		//rankinfo = (RankItem1<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem1<Imgidx, Pixel>));
		tmp = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc(max(hist_size * sizeof(Imgidx) * sizeof(Pixel), (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1));

		//int aa = max(hist_size, (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1);

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			//first row
			_COMPUTE_RANK_PIXEL //first pix
			pd[imgidx].cnt = 2;
			pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 2].alpha, rankinfo[contidx - 1].alpha);
			pd[imgidx].isVisited = 0;
			pd[imgidx].isAvailable = 0x3;
//			pd[imgidx + width].edgerank = 0xff;
			imgidx++;
//			pd[imgidx].edgerank = 0xfe;
			for (j = 1; j < width - 1; j++)
			{
				_COMPUTE_RANK_PIXEL
				pd[imgidx].cnt = 3;
				pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 2].alpha, rankinfo[contidx - 1].alpha, rankinfo[contidx - 3].alpha);
				pd[imgidx].isVisited = 0;
				pd[imgidx].isAvailable = 0x7;
				imgidx++;
			}
			_COMPUTE_RANK_PIXEL_LASTCOL
			pd[imgidx].cnt = 2;
			pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 1].alpha, maxpixval, rankinfo[contidx - 2].alpha);
			pd[imgidx].isVisited = 0;
			pd[imgidx].isAvailable = 0x5;
			imgidx++;

			for (i = 1; i < height - 1; i++)
			{
				_COMPUTE_RANK_PIXEL
				pd[imgidx].cnt = 3;
				pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 2].alpha, rankinfo[contidx - 1].alpha, maxpixval, rankinfo[contidx - 2 - wstride].alpha);
				pd[imgidx].isVisited = 0;
				pd[imgidx].isAvailable = 0xb;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					_COMPUTE_RANK_PIXEL
					pd[imgidx].cnt = 4;
					pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 2].alpha, rankinfo[contidx - 1].alpha, rankinfo[contidx - 3].alpha, rankinfo[contidx - 2 - wstride].alpha);
					pd[imgidx].isVisited = 0;
					pd[imgidx].isAvailable = 0xf;
					imgidx++;
				}
				_COMPUTE_RANK_PIXEL_LASTCOL
				pd[imgidx].cnt = 3;
				pd[imgidx].edgerank = compute_incidedge_queue(rankinfo[contidx - 1].alpha, maxpixval, rankinfo[contidx - 2].alpha, rankinfo[contidx - 1 - wstride].alpha);
				pd[imgidx].isVisited = 0;
				pd[imgidx].isAvailable = 0xd;
				imgidx++;
			}
			
			_COMPUTE_RANK_PIXEL_LASTROW
			pd[imgidx].cnt = 2;
			pd[imgidx].edgerank = compute_incidedge_queue(maxpixval, rankinfo[contidx - 1].alpha, maxpixval, rankinfo[contidx - 1 - wstride].alpha);
			pd[imgidx].isVisited = 0;
			pd[imgidx].isAvailable = 0xa;
			imgidx++;
			for (j = 1; j < width - 1; j++)
			{
				_COMPUTE_RANK_PIXEL_LASTROW
				pd[imgidx].cnt = 3;
				pd[imgidx].edgerank = compute_incidedge_queue(maxpixval, rankinfo[contidx - 1].alpha, rankinfo[contidx - 2].alpha, rankinfo[contidx + j - 1 - wstride].alpha);
				pd[imgidx].isVisited = 0;
				pd[imgidx].isAvailable = 0xe;
				imgidx++;
			}
			pd[imgidx].cnt = 2;
			pd[imgidx].edgerank = compute_incidedge_queue(maxpixval, maxpixval, rankinfo[contidx - 1].alpha, rankinfo[contidx + j - 1 - wstride].alpha);
			pd[imgidx].isVisited = 0;
			pd[imgidx].isAvailable = 0xc;
			imgidx++;
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
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width - 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width - 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width + 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width - 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width - 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
			}
		}

		for (i = 0; i < (Imgidx)(sizeof(Pixel) * 65536 >> 1); i++)
			hist[i] = 0;

		if (bit_depth < 16) // count sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankinfo[i].alpha;
				hist[hidx]++;
			}

			h = hist;

			hsum = 0;
			for (i = 0; i < (Imgidx)(hist_size); i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = rankinfo[i].alpha;
				j = --h[hidx];
				tmp[j] = rankinfo[i]; //slow!
			}
			r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
		}
		else //Radix sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankinfo[i].alpha;
				hist[hidx & mask]++;

				shamt = (Pixel)16;
				h_offset = (Pixel)65536;
				for (nbyte = 2; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
				{
					hidx = hidx >> shamt;
					hist[h_offset + (hidx & mask)]++;
					h_offset += (Pixel)65536;
				}
			}

			h = hist;
			shamt = 0;
			for (nbyte = 0; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
			{
				hsum = 0;
				for (i = 0; i < (Imgidx)(hist_size); i++)
				{
					hsum += h[i];
					h[i] = hsum;
				}
				for (i = nredges - 1; i >= 0; i--)
				{
					hidx = (rankinfo[i].alpha >> shamt) & mask;
					j = --h[hidx];
					tmp[j] = rankinfo[i]; //slow!
				}
				r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
				shamt += 16;
				h += hist_size;
			}
		}
		for (i = 0; i < nredges; i++)
		{
			//rank_to_alpha[i] = rankinfo[i].alpha;
			rank[rankinfo[i].dimgidx] = i;
		}



		//std::cout << "in_func2: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		//Free(rankinfo);
		Free(tmp);
		Free(hist);
		Free(incidents);
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
//			for (i = 0; i < ((imgsize + 1) >> 1); i++)
//				isAvailable[i] = 0xff;
			for (i = 0; i < imgsize; i++)
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

	void init_isAvailable_par(uint8* isAvailable, int npartitions)
	{
		int32 i, j, k;
		Imgidx imgsize = width * height;
		Imgidx wstride = width / npartitions;
		Imgidx hstride = height / npartitions;

		init_isAvailable(isAvailable);

		if (connectivity == 4)
		{
			//hor partitions
			j = (hstride - 1) * width;
			for (i = 0; i < npartitions - 1; i++)
			{
				k = j + width;
				//set_field(isAvailable, j, 0xa);
				//set_field(isAvailable, j + width, 0x3);
				for (; j < k; j++)
				{
					set_field(isAvailable, j, 0xe);
					set_field(isAvailable, j + width, 0x7);
				}
				//set_field(isAvailable, j, 0xc);
				//set_field(isAvailable, j + width, 0x5);

				j += (hstride - 1) * width;
			}

			//ver partitions

			for (i = 0; i < npartitions - 1; i++)
			{
				j = (i + 1) * wstride - 1;
				//k = j + height;
				//set_field(isAvailable, j, 0x5);
				//set_field(isAvailable, j + 1, 0x3);
				for (; j < imgsize; j += width)
				{
					set_field(isAvailable, j, 0xd);
					set_field(isAvailable, j + 1, 0xb);
				}
				//set_field(isAvailable, j, 0xc);
				//set_field(isAvailable, j + width, 0xa);
			}

			k = 0;
			for (i = 0; i < height; i++)
			{
				for (j = 0; j < width; j++)
				{
					uint16 p = get_field(isAvailable, k++);
					std::cout << p << '\t';
				}
				std::cout << std::endl;
			}

			/*		//partition crossings - no need
					k = (hstride - 1) * width;
					for(i = 0;i < npartitions - 1;i++)
					{
						k--;
						for(j = 0;j < npartitions - 1;j++)
						{
							k += wstride;
							set_field(isAvailable, k, 0xc);
							set_field(isAvailable, k + 1, 0xa);
							set_field(isAvailable, k + width, 0x5);
							set_field(isAvailable, k + width + 1, 0x3);
						}
						k += (hstride - 1) * width;
					}
					*/
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
		arr[idx] = in;
	}

	inline uint8 get_field(uint8* arr, Imgidx idx)
	{
		return arr[idx];
	}

	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		if (pNode->area) //possibly unnecessary branch..
			pNode->add(pix_val, pix_type);
		else
			pNode->set(1, level, (double)pix_val, pix_val, pix_val);
	}

	//#if DELAYED_NODE_ALLOC
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx *levelroot, int64 level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		Imgidx iNode = levelroot[level];
		if (iNode >= NODE_CANDIDATE)
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
			pNode->add(pix_val, pix_type);
		}
	}
	//#else
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode)
	{
		AlphaNode<Imgidx, Pixel> *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->add(pix_val, pix_type);
	}

	inline void connectPix2Node0(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		pNode->set(1, level, (double)pix_val, pix_val, pix_val, pix_type);
	}

	inline void connectNode2Node(Imgidx prev_top, Imgidx iPar, Pixel level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		pChild = node + prev_top;
		pPar = node + iPar;
		pChild->parentidx = iPar;
		if (pPar->area)
			pPar->add(pChild);
		else
		{
			pPar->alpha = level;
			pPar->copy(pChild);
		}
	}

	inline void connectNode2Node(Imgidx prev_top, Imgidx iPar)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		pChild = node + prev_top;
		pPar = node + iPar;
		pChild->parentidx = iPar;
		if (pPar->area)
			pPar->add(pChild);
		else
			pPar->copy(pChild);
	}

	inline void connectNode2Node(Imgidx* levelroot, Imgidx prev_top, int64 level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		Imgidx iPar = levelroot[level];
		if (iPar >= NODE_CANDIDATE)
		{
			iPar = NewAlphaNode();
			levelroot[level] = iPar;
			pPar = node + iPar;
			pChild = node + prev_top;
			pChild->parentidx = iPar;
			pPar->alpha = level;
			pPar->copy(pChild);
		}
		else
		{
			pPar = node + iPar;
			pChild = node + prev_top;
			pChild->parentidx = iPar;
			pPar->add(pChild);
		}
	}

	inline Imgidx NewAlphaNode()
	{
		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = min(2 * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		}
		return curSize++;
	}
	inline Imgidx NewAlphaNode(Pixel level, AlphaNode<Imgidx, Pixel> *pCopy)
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = min(2 * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			pNew = node + curSize;
		}
		pNew->alpha = level;
		pNew->copy(pCopy);
		return curSize++;
	}
	//#else
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
		pNew->alpha = level;
		pNew->minPix = (uint8)-1;
		pNew->maxPix = 0;
		pNew->sumPix = 0.0;
		//		pNew->parentidx = 0;
		pNew->area = 0;

		return curSize++;
	}
	//#endif
	inline uint8 is_visited(uint8* isVisited, Imgidx p)
	{
		return (isVisited[p >> 3] >> (p & 7)) & 1;
	}
	inline void visit(uint8* isVisited, Imgidx p)
	{
		isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
	}

#define _PUSH_NEIGHBOURS_4N(rettype) \
	{q = p << 1;\
	(is_available(isAv, 0) && !is_visited(isVisited, p + width)) ? queue->push(p + width, dimg[q]) : (rettype)0;\
	(is_available(isAv, 1) && !is_visited(isVisited, p + 1)) ? queue->push(p + 1, dimg[q + 1]) : (rettype)0;\
	(is_available(isAv, 2) && !is_visited(isVisited, p - 1)) ? queue->push(p - 1, dimg[q - 1]) : (rettype)0;\
	(is_available(isAv, 3) && !is_visited(isVisited, p - width)) ? queue->push(p - width, dimg[q - (width << 1)]) : (rettype)0;}

#define _PUSH_NEIGHBOURS_8N(rettype) \
	{q = p << 2;\
	(is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))?queue->push(p + wstride1, dimg[q]) : (rettype)0;\
	(is_available(isAv, 1) && !is_visited(isVisited, p + width))   ?queue->push(p + width, dimg[q + 1]) : (rettype)0;\
	(is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))?queue->push(p + wstride0, dimg[q + 2]) : (rettype)0;\
	(is_available(isAv, 3) && !is_visited(isVisited, p + 1))	   ?queue->push(p + 1, dimg[q + 3]) : (rettype)0;\
	(is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))?queue->push(p - wstride1, dimg[q - wstride_d + 4]) : (rettype)0;\
	(is_available(isAv, 5) && !is_visited(isVisited, p - width))   ?queue->push(p - width, dimg[q - wstride_d + 1]) : (rettype)0;\
	(is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))?queue->push(p - wstride0, dimg[q - wstride_d - 2]) : (rettype)0;\
	(is_available(isAv, 7) && !is_visited(isVisited, p - 1))	   ?queue->push(p - 1, dimg[q - 1]) : (rettype)0;}

#define _CREATE_NEW_NODE(minlev) \
	{Pixel pix_val = img[p];\
	current_level = minlev;\
	iNode = NewAlphaNode();\
	node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);\
	node[iNode].parentidx = stack_top;\
	stack_top = iNode;}

#define  _CREATE_SINGLETON_NODE_0 \
	{iNode = NewAlphaNode(0, node + stack_top);\
	node[iNode].parentidx = stack_top;\
	parentAry[p] = iNode;\
	prev_top = iNode;}\

#define  _CREATE_SINGLETON_NODE_1 \
	{Pixel pix_val = img[p];\
	iNode = NewAlphaNode();\
	node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);\
	node[stack_top].add(node + iNode);\
	node[iNode].parentidx = stack_top;\
	parentAry[p] = iNode;}

#define  _DECLARE_COMMON_LOCALS \
	Imgidx imgsize, dimgsize, nredges, x0;\
	int64 numlevels, max_level, current_level;\
	Imgidx *dhist;\
	Pixel *dimg;\
	Imgidx prev_top, stack_top, iNode;\
	uint8 *isVisited, *isAvailable, isAv;\
	Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;\
	imgsize = width * height;\
	nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);\
	dimgsize = (connectivity >> 1) * width * height;\
	numlevels = 1 << bit_depth;\
	max_level = numlevels - 1;
	
	Imgidx TreeSizeEstimation(Imgidx *dhist, int64 numlevels, Imgidx imgsize, Imgidx nredges)
	{
		nrmsd = 0;
		for (Imgidx p = 0; p < numlevels; p++)
			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		return min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
	}

#define _SET_COMMON_MEMORY	\
	Free(dhist);\
	isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));\
	isAvailable = (uint8*)Malloc((size_t)(imgsize));\
	memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));\
	init_isAvailable(isAvailable);\
	parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));\
	node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));


#define _INIT_FLOOD \
	stack_top = NewAlphaNode();/*dummy root*/\
	AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;\
	pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);\
	pNode->parentidx = stack_top;\
	current_level = max_level;\
	x0 = 0; /*arbitrary starting point*/\
	queue->push(x0, current_level);\
	prev_top = stack_top; /*to find redundant node*/\

	inline void remove_redundant_node(Imgidx& prev_top, Imgidx& stack_top)
	{
		if (node[prev_top].parentidx == stack_top && node[prev_top].area == node[stack_top].area)
		{
			node[prev_top].parentidx = node[stack_top].parentidx;
			stack_top = prev_top;
			curSize--;
		}
	}

#define _CREATE_NEW_STACKTOP(minlev) \
	{iNode = NewAlphaNode(minlev, node + stack_top);\
	node[iNode].parentidx = node[stack_top].parentidx;\
	node[stack_top].parentidx = iNode;}

	void Flood_HQueue(Pixel* img)
	{
		HierarQueue<Imgidx>* queue;

		_DECLARE_COMMON_LOCALS

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new HierarQueue<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
		curSize = 0;

		maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);

		_SET_COMMON_MEMORY
		_INIT_FLOOD
		queue->push(x0, current_level);
		while (1) //flooding
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (is_visited(isVisited, p))
				{
					queue->find_minlev();
					continue;
				}
				visit(isVisited, p);

				isAv = isAvailable[p];
				if (connectivity == 4) _PUSH_NEIGHBOURS_4N(void)
				else 				   _PUSH_NEIGHBOURS_8N(void)

				if (current_level > queue->min_level) //go to lower level
				{
					_CREATE_NEW_NODE(queue->min_level)
					if (current_level)
						_CREATE_SINGLETON_NODE_0
					else
						parentAry[p] = stack_top;
				}
				else
				{
					queue->find_minlev();
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
				}
			}
			
			remove_redundant_node(prev_top, stack_top);		

			//go to higher level
			iNode = node[stack_top].parentidx;
			if (queue->min_level < node[iNode].alpha) //new level from queue
				_CREATE_NEW_STACKTOP(queue->min_level)
			else //go to existing node
				node[iNode].add(node + stack_top);
			
			if (node[iNode].area == imgsize)	// root node found...done	
				break;

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = rootidx;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}
	void Flood_HQueue1(Pixel* img)
	{
		HQueue_l1idx<Imgidx>* queue;

		_DECLARE_COMMON_LOCALS

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new HQueue_l1idx<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
		curSize = 0;

		maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);

		_SET_COMMON_MEMORY
		_INIT_FLOOD
		queue->push(x0, current_level);
		while (1) //flooding
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (is_visited(isVisited, p))
				{
					queue->find_minlev();
					continue;
				}
				visit(isVisited, p);

				isAv = isAvailable[p];
				if (connectivity == 4) _PUSH_NEIGHBOURS_4N(int)
				else 				   _PUSH_NEIGHBOURS_8N(int)

				if (current_level > queue->min_level) //go to lower level
				{
					_CREATE_NEW_NODE(queue->min_level)
						if (current_level)
							_CREATE_SINGLETON_NODE_0
						else
							parentAry[p] = stack_top;
				}
				else
				{
					queue->find_minlev();
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
				}
			}

			remove_redundant_node(prev_top, stack_top);

			//go to higher level
			iNode = node[stack_top].parentidx;
			if (queue->min_level < node[iNode].alpha) //new level from queue
				_CREATE_NEW_STACKTOP(queue->min_level)
			else //go to existing node
				node[iNode].add(node + stack_top);

			if (node[iNode].area == imgsize)	// root node found...done	
				break;

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = rootidx;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}
	void Flood_HQueue2(Pixel* img)
	{
		HQueue_l2idx<Imgidx>* queue;

		_DECLARE_COMMON_LOCALS

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new HQueue_l2idx<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
		curSize = 0;

		maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);

		_SET_COMMON_MEMORY
		_INIT_FLOOD
		queue->push(x0, current_level);
		while (1) //flooding
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (is_visited(isVisited, p))
				{
					queue->find_minlev();
					continue;
				}
				visit(isVisited, p);

				isAv = isAvailable[p];
				if (connectivity == 4) _PUSH_NEIGHBOURS_4N(void)
				else 				   _PUSH_NEIGHBOURS_8N(void)

				if (current_level > queue->min_level) //go to lower level
				{
					_CREATE_NEW_NODE(queue->min_level)
						if (current_level)
							_CREATE_SINGLETON_NODE_0
						else
							parentAry[p] = stack_top;
				}
				else
				{
					queue->find_minlev();
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
				}
			}

			remove_redundant_node(prev_top, stack_top);

			//go to higher level
			iNode = node[stack_top].parentidx;
			if (queue->min_level < node[iNode].alpha) //new level from queue
				_CREATE_NEW_STACKTOP(queue->min_level)
			else //go to existing node
				node[iNode].add(node + stack_top);

			if (node[iNode].area == imgsize)	// root node found...done	
				break;

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = rootidx;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}
	void Flood_HeapQueue(Pixel* img)
	{
		HeapQueue<Imgidx, Pixel>* queue;

		_DECLARE_COMMON_LOCALS

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		//create heap-based priority queue from dhist
		queue = new HeapQueue<Imgidx, Pixel>(nredges);
		curSize = 0;

		maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);

		_SET_COMMON_MEMORY
		_INIT_FLOOD
		queue->push_run(x0, current_level);
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (is_visited(isVisited, p))
				{
					queue->pop_run();
					continue;
				}
				visit(isVisited, p);

				isAv = isAvailable[p];
				if (connectivity == 4) _PUSH_NEIGHBOURS_4N(void)
				else                   _PUSH_NEIGHBOURS_8N(void)

				queue->find_minlev();
				if (current_level > queue->get_minlev())
				{
					_CREATE_NEW_NODE(queue->get_minlev())
						if (current_level)
							_CREATE_SINGLETON_NODE_0
						else
							parentAry[p] = stack_top;
				}
				else
				{
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
				}

			}
			remove_redundant_node(prev_top, stack_top);

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
				_CREATE_NEW_STACKTOP(queue->get_minlev())
			else
				node[iNode].add(node + stack_top);

			if (node[iNode].area == imgsize)	// root node found...done	
				break;

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = rootidx;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

#define _RANK_SET_COMMON_LOCALS \
	Imgidx imgsize, dimgsize, nredges;																						\
	Imgidx current_rank, next_rank;																							\
	RankItem<Imgidx, Pixel> *rankinfo, *pRank;																				\
	AlphaNode<Imgidx, Pixel> *pNode;																						\
	Pixel maxpixval;																										\
	Imgidx *rank, top_rank;																									\
	int8 nbits;																												\
	Imgidx *dhist;																											\
	Imgidx prev_top;																										\
	uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;																\
	Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;										\
	imgsize = width * height;																								\
	nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);	\
	dimgsize = (connectivity >> 1) * width * height;																		\
	isVisited = (uint8*)Calloc((size_t)((imgsize)));																		\
	isAvailable = (uint8*)Malloc((size_t)(imgsize));																		\
	init_isAvailable(isAvailable);																							\
	rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));									\
	rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));																\
	maxSize = imgsize + nredges;																							\
	num_node = maxSize;																										\
	num_node_in = nredges;																									\
	parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));															\
	node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));							\
	node_in = node + imgsize;																								\
	nbits = ((sizeof(Pixel) << 3) - 1);																						\
	maxpixval = ~(1 << nbits);

#define _RANK_INIT_FLOOD(rettype)													\
	isVisited[0] = 1;																\
	(connectivity == 4)?queue->push(rank[0]):(rettype)0;							\
	(connectivity == 4)?queue->push(rank[1]):(rettype)0;							\
	(connectivity == 8)?queue->push(rank[1]):(rettype)0;							\
	(connectivity == 8)?queue->push(rank[2]):(rettype)0;							\
	(connectivity == 8)?queue->push(rank[3]):(rettype)0;							\
	current_rank = queue->top();													\
	node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize);		\
	prev_top = current_rank;

	void initialize_node(Pixel *img, RankItem<Imgidx, Pixel> *rankinfo, Pixel maxpixval)
	{
		AlphaNode<Imgidx, Pixel> *pNode;
		RankItem<Imgidx, Pixel> *pRank;
		Imgidx p;

		for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
			pNode->set(1, 0, (double)img[p], img[p], img[p]);		//init leaf nodes
		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0, pix_type); //init inner nodes
	}

#define	_RANK_PUSH_NEIGHBOURS_4N(rettype) {														\
	q = p << 1;																					\
	(is_available(isAv, 0) && !isVisited[p + width]) ? queue->push(rank[q]) : (rettype)0;		\
	(is_available(isAv, 1) && !isVisited[p + 1]) ? queue->push(rank[q + 1]) : (rettype)0;		\
	(is_available(isAv, 2) && !isVisited[p - 1]) ? queue->push(rank[q - 1]) : (rettype)0;		\
	(is_available(isAv, 3) && !isVisited[p - width]) ? queue->push(rank[q - (width << 1)]) : (rettype)0;}

#define	_RANK_PUSH_NEIGHBOURS_8N(rettype) {																	\
	q = p << 2;																								\
	(is_available(isAv, 0) && !isVisited[p + wstride1]) ? queue->push(rank[q]) : (rettype)0;				\
	(is_available(isAv, 1) && !isVisited[p + width]) ? queue->push(rank[q + 1]) : (rettype)0;				\
	(is_available(isAv, 2) && !isVisited[p + wstride0]) ? queue->push(rank[q + 2]) : (rettype)0;			\
	(is_available(isAv, 3) && !isVisited[p + 1]) ? queue->push(rank[q + 3]) : (rettype)0;					\
	(is_available(isAv, 4) && !isVisited[p - wstride1]) ? queue->push(rank[q - wstride_d + 4]) : (rettype)0;\
	(is_available(isAv, 5) && !isVisited[p - width]) ? queue->push(rank[q - wstride_d + 1]) : (rettype)0;	\
	(is_available(isAv, 6) && !isVisited[p - wstride0]) ? queue->push(rank[q - wstride_d - 2]) : (rettype)0;\
	(is_available(isAv, 7) && !isVisited[p - 1]) ? queue->push(rank[q - 1]) : (rettype)0;}

#define _RANK_REMOVE_REDUNDANT_NODE \
	(node_in[prev_top].parentidx == current_rank && node_in[prev_top].area == node_in[current_rank].area) ? current_rank = prev_top : (Imgidx)0;
	
	void Flood_HierarQueue_Rank(Pixel* img)
	{
		HybridQueue_HQueue_Rank1<Imgidx> *queue;

		_RANK_SET_COMMON_LOCALS

		queue = new HybridQueue_HQueue_Rank1<Imgidx>(nredges + 1, 32);

		compute_dimg(rank, rankinfo, img);
		initialize_node(img, rankinfo, maxpixval);

		_RANK_INIT_FLOOD(void)

			while (1)
			{
				while (1)
				{
					top_rank = queue->top();	//remove tmp variables later if possible
					pRank = rankinfo + top_rank;
					if (isVisited[pRank->p])
					{
						if (isVisited[pRank->q])
							break;
						p = pRank->q;
					}
					else
						p = pRank->p;
					isVisited[p] = 1;
					isAv = isAvailable[p];
					if (connectivity == 4)	_RANK_PUSH_NEIGHBOURS_4N(void)
					else					_RANK_PUSH_NEIGHBOURS_8N(void)

						next_rank = queue->top();
					node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);
					if (current_rank == next_rank)
						break;
					current_rank = next_rank;
				}

				queue->pop();
				next_rank = queue->top();

				_RANK_REMOVE_REDUNDANT_NODE

					node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
				if (node_in[next_rank].area == imgsize)
					break;

				prev_top = current_rank;
				current_rank = next_rank;
			}
		if (node_in[current_rank].area == imgsize)
			next_rank = current_rank;
		//std::cout << "Mean jump dist: " << queue->jumpdist / queue->jumpnum << std::endl;
//		std::cout << "Mean jump dist: " << queue->get_jumpdist() / queue->get_jumpnum() << std::endl;

		delete queue;
		Free(rank);
		Free(rankinfo);
		Free(isVisited);
		Free(isAvailable);
	}

	


	void Flood_HeapQueue_Rank(Pixel* img)
	{
		HeapQueue_rank<Imgidx> *queue;

		_RANK_SET_COMMON_LOCALS

		queue = new HeapQueue_rank<Imgidx>(nredges + 1);

		compute_dimg(rank, rankinfo, img);
		initialize_node(img, rankinfo, maxpixval);

		_RANK_INIT_FLOOD(void)

		while (1)
		{
			while (1)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
				pRank = rankinfo + top_rank;
				if (isVisited[pRank->p])
				{
					if (isVisited[pRank->q])
						break;
					p = pRank->q;
				}
				else
					p = pRank->p;
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)	_RANK_PUSH_NEIGHBOURS_4N(void)
				else					_RANK_PUSH_NEIGHBOURS_8N(void)

					next_rank = queue->top();
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			_RANK_REMOVE_REDUNDANT_NODE

				node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
				break;

			prev_top = current_rank;
			current_rank = next_rank;
		}
		if (node_in[current_rank].area == imgsize)
			next_rank = current_rank;


		delete queue;
		Free(rank);
		Free(rankinfo);
		Free(isVisited);
		Free(isAvailable);
	}


	void Flood_Trie(Pixel* img)
	{
		Trie<Imgidx, trieidx> *queue;

		_RANK_SET_COMMON_LOCALS

		queue = new Trie<Imgidx, trieidx>(nredges);

		compute_dimg(rank, rankinfo, img);
		initialize_node(img, rankinfo, maxpixval);
		
		_RANK_INIT_FLOOD(int8)

		while (1)
		{
			while (1)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
				pRank = rankinfo + top_rank;
				if (isVisited[pRank->p])
				{
					if (isVisited[pRank->q])
						break;
					p = pRank->q;
				}
				else
					p = pRank->p;
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)	_RANK_PUSH_NEIGHBOURS_4N(int8)
				else					_RANK_PUSH_NEIGHBOURS_8N(int8)

				next_rank = queue->top();
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			_RANK_REMOVE_REDUNDANT_NODE
			
			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
				break;

			prev_top = current_rank;
			current_rank = next_rank;
		}
		if (node_in[current_rank].area == imgsize)
			next_rank = current_rank;

		delete queue;
		Free(rank);
		Free(rankinfo);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_HybridQueue(Pixel* img)
	{
		HybridQueue_Trie<Imgidx, trieidx> *queue;

		_RANK_SET_COMMON_LOCALS

		queue = new HybridQueue_Trie<Imgidx, trieidx>(nredges);

		compute_dimg(rank, rankinfo, img);
		initialize_node(img, rankinfo, maxpixval);

		_RANK_INIT_FLOOD(void)

		while (1)
		{
			while (1)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
				pRank = rankinfo + top_rank;
				if (isVisited[pRank->p])
				{
					if (isVisited[pRank->q])
						break;
					p = pRank->q;
				}
				else
					p = pRank->p;
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)	_RANK_PUSH_NEIGHBOURS_4N(void)
				else					_RANK_PUSH_NEIGHBOURS_8N(void)

					next_rank = queue->top();
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			_RANK_REMOVE_REDUNDANT_NODE

				node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
				break;

			prev_top = current_rank;
			current_rank = next_rank;
		}
		if (node_in[current_rank].area == imgsize)
			next_rank = current_rank;

		delete queue;
		Free(rank);
		Free(rankinfo);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_Trie_sorted(Pixel* img) //experimental
	{
		//HierarQueue<Imgidx> *queue;
		Trie<Imgidx, int64> *queue;

		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank, next_rank;
		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank;
		int8 nbits;
		Imgidx *dhist;
		Imgidx prev_top;
		//		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;		
		uint8 isAv;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		//		isVisited = (uint8*)Malloc((size_t)((imgsize)));																		
		//		isAvailable = (uint8*)Malloc((size_t)(imgsize));																		
		//memset(isVisited, 0, (size_t)((imgsize)));
		//		init_isAvailable(isAvailable);																							
		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		maxSize = imgsize + nredges;
		num_node = maxSize;
		num_node_in = nredges;
		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;
		nbits = ((sizeof(Pixel) << 3) - 1);
		maxpixval = ~(1 << nbits);

		pixdata *pd;

		queue = new Trie<Imgidx, trieidx>(nredges);

		pd = (pixdata *)Malloc(width * height * sizeof(pixdata));


		compute_dimg(rank, rankinfo, pd, img); //yogi. performance penalty negligible.
		initialize_node(img, rankinfo, maxpixval);

		pd[0].isVisited = 1;
		(connectivity == 4) ? queue->push(rank[0]) : (int8)0;
		(connectivity == 4) ? queue->push(rank[1]) : (int8)0;
		(connectivity == 8) ? queue->push(rank[1]) : (int8)0;
		(connectivity == 8) ? queue->push(rank[2]) : (int8)0;
		(connectivity == 8) ? queue->push(rank[3]) : (int8)0;
		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize);
		prev_top = current_rank;

		Imgidx poffset[4] = { width, 1, -1, -width };
		Imgidx qoffset[4] = { 0, 1, -1, -(width << 1) };
		uint8 mask = 0x3, erank, eidx, cnt;
		while (1)
		{
			while (1)
			{
				//top_rank = queue->top();	//remove tmp variables later if possible
				pRank = rankinfo + queue->top();
				if (!pd[pRank->p].cnt)
				{
					if (!pd[pRank->q].cnt)
						break;
					p = pRank->q;
				}
				else
					p = pRank->p;
				// 				if (!pd[p].cnt)
				// 					break;


								//isAv = pd[p].isAvailable;
				if (connectivity == 4) {
					q = p << 1;

					erank = pd[p].edgerank;
					cnt = pd[p].cnt - 1;
					eidx = erank & mask;
					next_rank = rank[q + qoffset[eidx]];
					if (!pd[p].isVisited)
					{
						node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);
						pd[p].isVisited = 1;
					}
					//should break here...
					(pd[p + poffset[eidx]].isVisited || queue->push(next_rank)) &&
						cnt && cnt-- && (eidx = (erank >> 2) & mask) < connectivity && (pd[p + poffset[eidx]].isVisited || queue->push(rank[q + qoffset[eidx]])) &&
						cnt && cnt-- && (eidx = (erank >> 4) & mask) < connectivity && (pd[p + poffset[eidx]].isVisited || queue->push(rank[q + qoffset[eidx]])) &&
						cnt && cnt-- && (eidx = (erank >> 6) & mask) < connectivity && (pd[p + poffset[eidx]].isVisited || queue->push(rank[q + qoffset[eidx]]));
					pd[p].edgerank = erank >> ((pd[p].cnt - cnt) << 1);
					pd[p].cnt = cnt;
					//for ()
					//	;
// 						(is_available(isAv, 0) && !pd[p + width].isVisited) ? queue->push(rank[q], 1) : (int8)0; 
// 						(is_available(isAv, 1) && !pd[p + 1].isVisited) ? queue->push(rank[q + 1], 1) : (int8)0;
// 						(is_available(isAv, 2) && !pd[p - 1].isVisited) ? queue->push(rank[q - 1], 0) : (int8)0;
// 						(is_available(isAv, 3) && !pd[p - width].isVisited) ? queue->push(rank[q - (width << 1)], 0) : (int8)0;
				}
				//				else					_RANK_PUSH_NEIGHBOURS_8N

				// 				next_rank = queue->top();
				// 				if (current_rank == next_rank)
				// 					break;
				//				current_rank = next_rank;
				current_rank = queue->top();
			}

			queue->pop();
			next_rank = queue->top();

			_RANK_REMOVE_REDUNDANT_NODE

				node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
				break;

			prev_top = current_rank;
			current_rank = next_rank;
		}
		if (node_in[current_rank].area == imgsize)
			next_rank = current_rank;

		delete queue;
		Free(rank);
		Free(rankinfo);
		//		Free(isVisited);
		//		Free(isAvailable);
		Free(pd);
	}

public:
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Imgidx, Pixel> *node, *node_in;
	Imgidx num_node, num_node_in;
	Imgidx rootidx;
	Imgidx* parentAry;

	//double ratio_excessnodes;
	double nrmsd;
	int8 pix_type;
	int8 bit_depth;

	ATree(int8 ptype, int8 bit_depth) : maxSize(0), curSize(0), node(0), parentAry(0), pix_type(ptype), bit_depth(bit_depth) {}
	~ATree() { if (node) { Free(node); Free(parentAry); } }
	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }

	void BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, int connectivity_in, int algorithm)
	{
		this->height = (Imgidx)height_in;
		this->width = (Imgidx)width_in;
		this->channel = (Imgidx)channel_in;
		this->connectivity = (Imgidx)connectivity_in;
		curSize = 0;
		if (connectivity != 4 && connectivity != 8)
		{
			//std::cout << "connectivity should be 4 or 8\n" << std::endl;
			return;
		}

		switch (algorithm)
		{
		case(0): Flood_HQueue(img);	        break;
		case(1): Flood_HQueue1(img);        break;
		case(2): Flood_HQueue2(img);	    break;
		case(3): Flood_HeapQueue(img);      break;
		case(4): Flood_HeapQueue_Rank(img);	break;
		case(5): Flood_Trie(img);		    break;
		case(6): Flood_HybridQueue(img);	break;
		case(7): Flood_HierarQueue_Rank(img);	break;
		default: break;
		}
	}

	void AreaFilter(Pixel *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		AlphaNode<Imgidx, Pixel> *pNode;

		for (idx_lim = 0; idx_lim < num_node_in && node_in[idx_lim].alpha < alpha; idx_lim++)
			;

		for (i = 0; i < num_node; i++)
		{
			node[i].rootidx = 0;
			node[i].filter_val = 0;
		}

		imgsize = height * width;
		iarea = (Imgidx)(area * (double)imgsize);
		iarea = min(imgsize, max(0, iarea));
		//val = 1;
		for (i = 0; i < imgsize; i++)
		{
			pNode = &node[i];
			while (pNode->parentidx != -1 && pNode->alpha < alpha)
				pNode = &node[pNode->parentidx];
			outimg[i] = min(255, (Pixel)pNode->area / 200.0);
		}
	}

	void AreaFilter(double *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		AlphaNode<Imgidx, Pixel> *pNode;

		for (idx_lim = 0; idx_lim < num_node_in && node_in[idx_lim].alpha < alpha; idx_lim++)
			;

		for (i = 0; i < num_node; i++)
		{
			node[i].rootidx = 0;
			node[i].filter_val = 0;
		}

		imgsize = height * width;
		iarea = (Imgidx)(area * (double)imgsize);
		iarea = min(imgsize, max(0, iarea));
		//val = 1;
		for (i = 0; i < imgsize; i++)
		{
			pNode = &node[i];
			while (pNode->parentidx != -1 && pNode->alpha < alpha)
				pNode = &node[pNode->parentidx];
			outimg[i] = min(255, (double)pNode->area / 200.0);
		}
	}

	void print_tree()
	{
		for (int i = 0; i < maxSize; i++)
			node[i].print(node);
	}

};


class AlphaTree
{
	void *tree;
	int8 imgidx, pix_type, bit_depth;
public:
	AlphaTree(int8 bit_depth) : tree(0), bit_depth(bit_depth) {}
	~AlphaTree()
	{
		clear();
	}

	inline void clear()
	{
		if (tree)
		{
			if (imgidx == IMGIDX_32BITS)
			{
				if (pix_type == PIXEL_8BIT)			delete ((ATree<int32, uint8>*)tree);
				else if (pix_type == PIXEL_16BIT)	delete ((ATree<int32, uint16>*)tree);
				else if (pix_type == PIXEL_32BIT)	delete ((ATree<int32, uint32>*)tree);
				else if (pix_type == PIXEL_64BIT)	delete ((ATree<int32, uint64>*)tree);
				else if (pix_type == PIXEL_FLOAT)	delete ((ATree<int32, uint32>*)tree);
				else								delete ((ATree<int32, uint64>*)tree);
			}
			else
			{
				if (pix_type == PIXEL_8BIT)			delete ((ATree<int64, uint8>*)tree);
				else if (pix_type == PIXEL_16BIT)	delete ((ATree<int64, uint16>*)tree);
				else if (pix_type == PIXEL_32BIT)	delete ((ATree<int64, uint32>*)tree);
				else if (pix_type == PIXEL_64BIT)	delete ((ATree<int64, uint64>*)tree);
				else if (pix_type == PIXEL_FLOAT)	delete ((ATree<int64, uint32>*)tree);
				else								delete ((ATree<int64, uint64>*)tree);
			}
		}
		tree = 0;
	}

	void BuildAlphaTree(uint8 *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_8BIT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = new ATree<int32, uint8>(pix_type, bit_depth);
			((ATree<int32, uint8>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint8>(pix_type, bit_depth);
			((ATree<int64, uint8>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}

	void BuildAlphaTree(uint16 *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_16BIT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = new ATree<int32, uint16>(pix_type, bit_depth);
			((ATree<int32, uint16>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint16>(pix_type, bit_depth);
			((ATree<int64, uint16>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}
	void BuildAlphaTree(uint32 *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_32BIT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = new ATree<int32, uint32>(pix_type, bit_depth);
			((ATree<int32, uint32>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint32>(pix_type, bit_depth);
			((ATree<int64, uint32>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}
	void BuildAlphaTree(uint64 *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_64BIT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = new ATree<int32, uint64>(pix_type, bit_depth);
			((ATree<int32, uint64>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint64>(pix_type, bit_depth);
			((ATree<int64, uint64>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}
	void BuildAlphaTree(float *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_FLOAT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;

			tree = new ATree<int32, uint32>(pix_type, bit_depth);
			((ATree<int32, uint32>*)tree)->BuildAlphaTree((uint32*)img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint32>(pix_type, bit_depth);
			((ATree<int64, uint32>*)tree)->BuildAlphaTree((uint32*)img, height, width, channel, connectivity, algorithm);
		}
	}
	void BuildAlphaTree(double *img, int height, int width, int channel, int connectivity, int algorithm)
	{

		//		printf("Entered BuildAlphaTree\n");
		pix_type = PIXEL_DOUBLE;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = new ATree<int32, uint64>(pix_type, bit_depth);
			tree = Malloc(sizeof(ATree<int32, uint64>));
			((ATree<int32, uint64>*)tree)->BuildAlphaTree((uint64*)img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = new ATree<int64, uint64>(pix_type, bit_depth);
			((ATree<int64, uint64>*)tree)->BuildAlphaTree((uint64*)img, height, width, channel, connectivity, algorithm);
		}
	}

	void AreaFilter(double *outimg, double area, double alpha)
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			((ATree<int32, uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
			else if (pix_type == PIXEL_16BIT)	((ATree<int32, uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
			else if (pix_type == PIXEL_32BIT)	((ATree<int32, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else if (pix_type == PIXEL_64BIT)	((ATree<int32, uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
			else if (pix_type == PIXEL_FLOAT)	((ATree<int32, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else								((ATree<int32, uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			((ATree<int64, uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
			else if (pix_type == PIXEL_16BIT)	((ATree<int64, uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
			else if (pix_type == PIXEL_32BIT)	((ATree<int64, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else if (pix_type == PIXEL_64BIT)	((ATree<int64, uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
			else if (pix_type == PIXEL_FLOAT)	((ATree<int64, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else								((ATree<int64, uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);
		}

	}
	void print_tree()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			((ATree<int32, uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)	((ATree<int32, uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)	((ATree<int32, uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)	((ATree<int32, uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)	((ATree<int32, uint32>*)tree)->print_tree();
			else								((ATree<int32, uint64>*)tree)->print_tree();
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			((ATree<int64, uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)	((ATree<int64, uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)	((ATree<int64, uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)	((ATree<int64, uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)	((ATree<int64, uint32>*)tree)->print_tree();
			else								((ATree<int64, uint64>*)tree)->print_tree();
		}
	}

	int64 get_maxSize()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			return (int64)((ATree<int32, uint8>*)tree)->maxSize;
			else if (pix_type == PIXEL_16BIT)	return (int64)((ATree<int32, uint16>*)tree)->maxSize;
			else if (pix_type == PIXEL_32BIT)	return (int64)((ATree<int32, uint32>*)tree)->maxSize;
			else if (pix_type == PIXEL_64BIT)	return (int64)((ATree<int32, uint64>*)tree)->maxSize;
			else if (pix_type == PIXEL_FLOAT)	return (int64)((ATree<int32, uint32>*)tree)->maxSize;
			else								return (int64)((ATree<int32, uint64>*)tree)->maxSize;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			return (int64)((ATree<int64, uint8>*)tree)->maxSize;
			else if (pix_type == PIXEL_16BIT)	return (int64)((ATree<int64, uint16>*)tree)->maxSize;
			else if (pix_type == PIXEL_32BIT)	return (int64)((ATree<int64, uint32>*)tree)->maxSize;
			else if (pix_type == PIXEL_64BIT)	return (int64)((ATree<int64, uint64>*)tree)->maxSize;
			else if (pix_type == PIXEL_FLOAT)	return (int64)((ATree<int64, uint32>*)tree)->maxSize;
			else								return (int64)((ATree<int64, uint64>*)tree)->maxSize;
		}
	}
	int64 get_curSize()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			return (int64)((ATree<int32, uint8>*)tree)->curSize;
			else if (pix_type == PIXEL_16BIT)	return (int64)((ATree<int32, uint16>*)tree)->curSize;
			else if (pix_type == PIXEL_32BIT)	return (int64)((ATree<int32, uint32>*)tree)->curSize;
			else if (pix_type == PIXEL_64BIT)	return (int64)((ATree<int32, uint64>*)tree)->curSize;
			else if (pix_type == PIXEL_FLOAT)	return (int64)((ATree<int32, uint32>*)tree)->curSize;
			else								return (int64)((ATree<int32, uint64>*)tree)->curSize;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			return (int64)((ATree<int64, uint8>*)tree)->curSize;
			else if (pix_type == PIXEL_16BIT)	return (int64)((ATree<int64, uint16>*)tree)->curSize;
			else if (pix_type == PIXEL_32BIT)	return (int64)((ATree<int64, uint32>*)tree)->curSize;
			else if (pix_type == PIXEL_64BIT)	return (int64)((ATree<int64, uint64>*)tree)->curSize;
			else if (pix_type == PIXEL_FLOAT)	return (int64)((ATree<int64, uint32>*)tree)->curSize;
			else								return (int64)((ATree<int64, uint64>*)tree)->curSize;
		}
	}

	double get_nrmsd()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			return (double)((ATree<int32, uint8>*)tree)->nrmsd;
			else if (pix_type == PIXEL_16BIT)	return (double)((ATree<int32, uint16>*)tree)->nrmsd;
			else if (pix_type == PIXEL_32BIT)	return (double)((ATree<int32, uint32>*)tree)->nrmsd;
			else if (pix_type == PIXEL_64BIT)	return (double)((ATree<int32, uint64>*)tree)->nrmsd;
			else if (pix_type == PIXEL_FLOAT)	return (double)((ATree<int32, uint32>*)tree)->nrmsd;
			else								return (double)((ATree<int32, uint64>*)tree)->nrmsd;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			return (double)((ATree<int64, uint8>*)tree)->nrmsd;
			else if (pix_type == PIXEL_16BIT)	return (double)((ATree<int64, uint16>*)tree)->nrmsd;
			else if (pix_type == PIXEL_32BIT)	return (double)((ATree<int64, uint32>*)tree)->nrmsd;
			else if (pix_type == PIXEL_64BIT)	return (double)((ATree<int64, uint64>*)tree)->nrmsd;
			else if (pix_type == PIXEL_FLOAT)	return (double)((ATree<int64, uint32>*)tree)->nrmsd;
			else								return (double)((ATree<int64, uint64>*)tree)->nrmsd;
		}
	}
};
