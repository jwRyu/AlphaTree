#pragma once
#pragma once
#include<iostream>
#include<math.h>
#include "defines.h"
#include "HQueue.h"
#include "Trie.h"
#include "HybridQueue.h"
#include "PriorityQueue.h"

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


using namespace std;
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
		/*
		if(pix_type < 4)
			this->sumPix += q->sumPix;
		else if (pix_type == PIXEL_DOUBLE)
			this->sumPix += *((double*)(&q->sumPix));
		else
			this->sumPix += *((float*)(&q->sumPix));
			*/

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

	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>*& rankinfo, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;
		int8 p64 = (pix_type == PIXEL_64BIT) || (pix_type == PIXEL_DOUBLE);

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
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
						: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 2;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				rankinfo[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))
					: (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
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

	/*
	inline void set_field(uint8* arr, Imgidx idx, uint8 in)
	{
		uint8 shamt = (idx & 1) << 2;
		arr[idx >> 1] &= (in << shamt) | ((uint8)(0x0f) << (4 - shamt));
	}

	inline uint8 get_field(uint8* arr, Imgidx idx)
	{
		return (arr[idx >> 1] >> ((idx & 1) << 2)) & 0x0f;
	}
	*/

	inline void push_neighbour(HQueue<Imgidx> *hqueue, Imgidx* levelroot, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		hqueue->push(idx, (int64)dissim);
		// 		if (levelroot[dissim] == NULL_LEVELROOT)
		// #if DELAYED_NODE_ALLOC
		// 			levelroot[dissim] = NODE_CANDIDATE;
		// #else
		// 			levelroot[dissim] = NewAlphaNode1(dissim);
		// #endif
	}

	inline void push_neighbour(HQueue<Imgidx> *hqueue, Imgidx* levelroot, int64* levelroot_flag, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		hqueue->push(idx, (int64)dissim);
		if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_NODE_ALLOC
			levelroot[dissim] = NODE_CANDIDATE;
		levelroot_flag[dissim >> 6] |= (int64)1 << (dissim & 63);
#else
			levelroot[dissim] = NewAlphaNode1(dissim);
#endif
	}

	inline void push_neighbour(HQueue<Imgidx> *hqueue, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		Imgidx idx_new;
		hqueue->push(idx, (int64)dissim);
	}

	inline void push_neighbour(HQueue_hdr<Imgidx> *hqueue, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		Imgidx idx_new;
		hqueue->push(idx, (int64)dissim);
	}

	inline void push_neighbour(HQueue_hdr2<Imgidx> *hqueue, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		Imgidx idx_new;
		hqueue->push(idx, (int64)dissim);
	}

	inline void push_neighbour(PriorityQueue<Imgidx, uint64> *pqueue, Imgidx* levelroot, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		pqueue->push(idx, dissim);
		if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_NODE_ALLOC
			levelroot[dissim] = NODE_CANDIDATE;
#else
			levelroot[dissim] = NewAlphaNode1(dissim);
#endif
	}

	inline void push_neighbour(PriorityQueue<Imgidx, Pixel> *pqueue, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		pqueue->push(idx, dissim);
	}

	inline void push_neighbour1(PriorityQueue<Imgidx, Pixel> *pqueue, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		pqueue->push_run(idx, dissim);
	}


	inline void push_neighbour(PriorityQueue<Imgidx, uint64> *pqueue, Imgidx* levelroot, int64* levelroot_flag, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
	{
		Pixel dissim = dimg[dimgidx];
		pqueue->push(idx, dissim);
		if (levelroot[dissim] == NULL_LEVELROOT)
#if DELAYED_NODE_ALLOC
			levelroot[dissim] = NODE_CANDIDATE;
		levelroot_flag[dissim >> 6] |= (int64)1 << (dissim & 63);
#else
			levelroot[dissim] = NewAlphaNode1(dissim);
#endif
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
	//#endif


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

	inline void connectNode2Node0(Imgidx prev_top, Imgidx iPar)
	{
		node[iPar].add(node + prev_top);
	}

	// 	inline void connectNode2Node0(Imgidx prev_top, Imgidx iPar, Pixel level)
	// 	{
	// 		node[iPar].alpha = level;
	// 		node[iPar].copy(node + prev_top);
	// 		node[prev_top].parentidx = iPar;
	// 	}

#if DELAYED_NODE_ALLOC
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
#else
	inline void connectNode2Node(AlphaNode<Imgidx, Pixel>* pPar, Imgidx iPar, AlphaNode<Imgidx, Pixel>* pNode)
	{
		pNode->parentidx = iPar;
		pPar->add(pNode);
	}
#endif
	//#if DELAYED_NODE_ALLOC
	inline Imgidx NewAlphaNode()
	{
		//		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = min(2 * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			//pNew = node + curSize;
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

	// 	void Flood_HQueue(Pixel* img)
	// 	{
	// 		Imgidx imgsize, dimgsize, nredges, x0;
	// 		int64 numlevels, max_level, current_level, next_level;
	// 		HQueue<Imgidx>* queue;
	// 		Imgidx *dhist;
	// 		Pixel *dimg;
	// 		Imgidx prev_top, *levelroot;
	// 		uint8 *isVisited, *isAvailable, isAv;;
	// 		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
	// 
	// 
	// 		imgsize = width * height;
	// 		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
	// 		dimgsize = (connectivity >> 1) * width * height;
	// 		numlevels = 1 << bit_depth;
	// 
	// 		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
	// 		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
	// 		levelroot = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	// 		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
	// 		if (connectivity == 4)
	// 			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
	// 		else
	// 			isAvailable = (uint8*)Malloc((size_t)(imgsize));
	// 		for (p = 0; p < numlevels; p++)
	// 			levelroot[p] = NULL_LEVELROOT;
	// 		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
	// 		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
	// 		init_isAvailable(isAvailable);
	// 
	// 		max_level = numlevels - 1;
	// 
	// 		compute_dimg(dimg, dhist, img);
	// 
	// 		dhist[max_level]++;
	// 		queue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels);
	// 		curSize = 0;
	// 
	// 		/////////////////////////////////////////
	// 		//tree size estimation (TSE)
	// 		nrmsd = 0;
	// 		for (p = 0; p < numlevels; p++)
	// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	// 		maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		/////////////////////////////////////////
	// 	//	printf("NRMSD: %f\tEst. NTS: %f\tEst. Tree size: %d\n", nrmsd, ((A * exp(SIGMA * nrmsd) + B) + M), tree->maxSize);
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		Free(dhist);
	// 
	// 		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
	// 		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
	// 
	// #if DELAYED_NODE_ALLOC
	// 		//		levelroot[max_level + 1] = NODE_CANDIDATE;
	// 		levelroot[max_level + 1] = NewAlphaNode();
	// 		AlphaNode<Imgidx, Pixel> *pNode = node + levelroot[max_level + 1];
	// 		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
	// 		//(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	// 
	// #else
	// 		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	// 		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
	// #endif
	// 
	// 		current_level = max_level;
	// 		x0 = imgsize >> 1;
	// 		queue->push(x0, current_level);
	// 
	// 		prev_top = levelroot[max_level + 1];
	// 		while (1)
	// 		{
	// 			while (queue->min_level <= current_level)
	// 			{
	// 				p = queue->pop();
	// 				if (is_visited(isVisited, p))
	// 				{
	// 					queue->find_min_level();
	// 					continue;
	// 				}
	// 				visit(isVisited, p);
	// 
	// #if !HQUEUE_COST_AMORTIZE
	// 				queue->find_min_level();
	// #endif
	// 				if (connectivity == 4)
	// 				{
	// 					isAv = get_field(isAvailable, p);
	// 					q = p << 1;
	// 
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(queue, levelroot, dimg, p + width, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))				push_neighbour(queue, levelroot, dimg, p + 1, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))				push_neighbour(queue, levelroot, dimg, p - 1, q - 1);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(queue, levelroot, dimg, p - width, q - (width << 1));
	// 				}
	// 				else
	// 				{
	// 					isAv = isAvailable[p];
	// 					q = p << 2;
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, levelroot, dimg, p + wstride1, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, levelroot, dimg, p + width, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, levelroot, dimg, p + wstride0, q + 2);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, levelroot, dimg, p + 1, q + 3);
	// 					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
	// 					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, levelroot, dimg, p - width, q - wstride_d + 1);
	// 					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
	// 					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, levelroot, dimg, p - 1, q - 1);
	// 				}
	// 
	// 				if (current_level > queue->min_level)
	// 					current_level = queue->min_level;
	// #if HQUEUE_COST_AMORTIZE
	// 				else
	// 					queue->find_min_level();
	// #endif
	// 
	// #if DELAYED_NODE_ALLOC
	// 				connectPix2Node(p, img[p], levelroot, current_level);
	// #else
	// 				connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
	// #endif
	// 
	// 			}
	// 			//std::cout << "checking redundancy..." << std::endl;
	// 			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
	// 			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
	// 			if (node[prev_top].parentidx == levelroot[current_level] &&
	// 				node[levelroot[current_level]].area == node[prev_top].area)
	// 			{
	// 				levelroot[current_level] = prev_top;
	// #if DELAYED_NODE_ALLOC
	// 				curSize--;
	// #endif
	// 			}
	// 
	// 			next_level = current_level + 1;
	// 			while((levelroot[next_level] == NULL_LEVELROOT) && (next_level != queue->get_minlev()))
	// 				next_level++;
	// 
	// 			//next_level = queue->get_minlev();
	// //			next_level = current_level + 1;
	// //			while ((levelroot[next_level] == NULL_LEVELROOT))
	// //				next_level++;
	// 
	// // 			double  sumpix = 0;
	// // 
	// // 			if (next_level == 256)
	// // 			{
	// // 				for (q = 0; q < curSize; q++)
	// // 				{
	// // 					if (node[q].area && node[q].parentidx == NULL_LEVELROOT )
	// // 					{
	// // 						sumpix += node[q].sumPix;
	// // 						q = q;
	// // 					}
	// // 				}
	// // 			}
	// // 			q = q;
	// 				
	// 
	// #if DELAYED_NODE_ALLOC
	// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
	// #else
	// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
	// #endif
	// 			if (node[levelroot[next_level]].area == imgsize)
	// 			{
	// 				if (node[levelroot[current_level]].area == imgsize)
	// 					rootidx = levelroot[current_level];
	// 				else
	// 					rootidx = levelroot[next_level];
	// 				//				node_in[next_rank].parentidx = next_rank;
	// 				// node_in[next_rank].parentidx = -1;
	// 				node[rootidx].parentidx = rootidx;
	// 				break;
	// 			}
	// 
	// 			prev_top = levelroot[current_level];
	// 			levelroot[current_level] = NULL_LEVELROOT;
	// 			current_level = next_level;
	// 
	// 		}
	// 		node[prev_top].parentidx = prev_top;
	// 		rootidx = prev_top;
	// 		curSize--;
	// 
	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}
	// 
	// 		delete queue;
	// 		Free(dimg);
	// 		Free(levelroot);
	// 		Free(isVisited);
	// 		Free(isAvailable);
	// 	}



	// 	int verify(Imgidx *levelroot, int64 *levelroot_flag, Imgidx size)
	// 	{
	// 		for (Imgidx i = 0; i < size; i++)
	// 		{
	// 			int lvlroot = levelroot[i] != NULL_LEVELROOT;
	// 			int flag = (levelroot_flag[i >> 6] >> (i & 63)) & 1;
	// 			if (lvlroot != flag)
	// 				return 0;
	// 		}
	// 		return 1;
	// 	}
	// 
	// 	void Flood_HQueue_doubleQ(Pixel* img)
	// 	{
	// 		Imgidx imgsize, dimgsize, nredges, x0;
	// 		int64 numlevels, max_level, current_level, next_level;
	// 		HQueue<Imgidx>* queue;
	// 		Imgidx *dhist;
	// 		Pixel *dimg;
	// 		Imgidx prev_top, *levelroot;
	// 		int64 *levelroot_flag;
	// 		uint8 *isVisited, *isAvailable, isAv;;
	// 		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
	// 
	// 
	// 		imgsize = width * height;
	// 		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
	// 		dimgsize = (connectivity >> 1) * width * height;
	// 		numlevels = 1 << bit_depth;
	// 
	// 		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
	// 		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
	// 		levelroot = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	// 		levelroot_flag = (int64*)Malloc((size_t)((numlevels + 1 + 63) >> 6) * sizeof(int64));
	// 		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
	// 		if (connectivity == 4)
	// 			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
	// 		else
	// 			isAvailable = (uint8*)Malloc((size_t)(imgsize));
	// 		for (p = 0; p < numlevels; p++)
	// 			levelroot[p] = NULL_LEVELROOT;
	// 		for (p = 0; p < (numlevels + 1 + 63) >> 6; p++)
	// 			levelroot_flag[p] = 0;
	// 		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
	// 		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
	// 		init_isAvailable(isAvailable);
	// 
	// 		max_level = numlevels - 1;
	// 
	// 		compute_dimg(dimg, dhist, img);
	// 
	// 		dhist[max_level]++;
	// 		queue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels);
	// 		curSize = 0;
	// 
	// 		/////////////////////////////////////////
	// 		//tree size estimation (TSE)
	// 		nrmsd = 0;
	// 		for (p = 0; p < numlevels; p++)
	// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	// 		maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		/////////////////////////////////////////
	// 	//	printf("NRMSD: %f\tEst. NTS: %f\tEst. Tree size: %d\n", nrmsd, ((A * exp(SIGMA * nrmsd) + B) + M), tree->maxSize);
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		Free(dhist);
	// 
	// 		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
	// 		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
	// 
	// #if DELAYED_NODE_ALLOC
	// 		//		levelroot[max_level + 1] = NODE_CANDIDATE;
	// 		levelroot[max_level + 1] = NewAlphaNode();
	// 		levelroot_flag[(max_level + 1) >> 6] |= (int64)1 << ((max_level + 1) & 63);
	// 		AlphaNode<Imgidx, Pixel> *pNode = node + levelroot[max_level + 1];
	// 		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
	// 		//(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	// 
	// #else
	// 		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	// 		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
	// #endif
	// 
	// 		current_level = max_level;
	// 		x0 = imgsize >> 1;
	// 		queue->push(x0, current_level);
	// 
	// 		prev_top = levelroot[max_level + 1];
	// 		while (current_level <= max_level)
	// 		{
	// 			while (queue->min_level <= current_level)
	// 			{
	// 				p = queue->pop();
	// 				if (is_visited(isVisited, p))
	// 				{
	// 					queue->find_min_level();
	// 					continue;
	// 				}
	// 				visit(isVisited, p);
	// #if !HQUEUE_COST_AMORTIZE
	// 				queue->find_min_level();
	// #endif
	// 				if (connectivity == 4)
	// 				{
	// 					isAv = get_field(isAvailable, p);
	// 					q = p << 1;
	// 
	// 					if (q == 0x0081c6aa)
	// 						q = q;
	// 
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p + width, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))				push_neighbour(queue, levelroot, levelroot_flag, dimg, p + 1, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))				push_neighbour(queue, levelroot, levelroot_flag, dimg, p - 1, q - 1);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p - width, q - (width << 1));
	// 				}
	// 				else
	// 				{
	// 					isAv = isAvailable[p];
	// 					q = p << 2;
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p + wstride1, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, levelroot, levelroot_flag, dimg, p + width, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p + wstride0, q + 2);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, levelroot, levelroot_flag, dimg, p + 1, q + 3);
	// 					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p - wstride1, q - wstride_d + 4);
	// 					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, levelroot, levelroot_flag, dimg, p - width, q - wstride_d + 1);
	// 					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, levelroot, levelroot_flag, dimg, p - wstride0, q - wstride_d - 2);
	// 					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, levelroot, levelroot_flag, dimg, p - 1, q - 1);
	// 				}
	// 
	// 				if (current_level > queue->min_level)
	// 					current_level = queue->min_level;
	// #if HQUEUE_COST_AMORTIZE
	// 				else
	// 					queue->find_min_level();
	// #endif
	// 
	// #if DELAYED_NODE_ALLOC
	// 				connectPix2Node(p, img[p], levelroot, current_level);
	// #else
	// 				connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
	// #endif
	// 
	// 			}
	// 			//std::cout << "checking redundancy..." << std::endl;
	// 			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
	// 			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
	// 			if (node[prev_top].parentidx == levelroot[current_level] &&
	// 				node[levelroot[current_level]].area == node[prev_top].area)
	// 			{
	// 				levelroot[current_level] = prev_top;
	// #if DELAYED_NODE_ALLOC
	// 				curSize--;
	// #endif
	// 			}
	// 			
	// 			Imgidx nlvl;
	// 			if (current_level == 42)
	// 				current_level = 42;
	// 			prev_top = levelroot[current_level];
	// 			levelroot[current_level] = NULL_LEVELROOT;
	// 			levelroot_flag[current_level >> 6] &= ~((int64)1 << (current_level & 63));
	// 	
	// 			next_level = (current_level + 1) >> 6;
	// 			while (!levelroot_flag[next_level])
	// 				next_level++;
	// 			next_level = next_level << 6;
	// 			while (levelroot[next_level] == NULL_LEVELROOT)
	// 				next_level++;
	// //			nlvl = next_level;
	// //			next_level = current_level + 1;
	// //			while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
	// //				next_level++;
	// 
	// //			if (nlvl != next_level)
	// //				nlvl = nlvl;
	// 
	// //			verify(levelroot, levelroot_flag, max_level + 1);
	// 
	// #if DELAYED_NODE_ALLOC
	// 			connectNode2Node(levelroot, prev_top, next_level);
	// #else
	// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
	// #endif
	// 					   			
	// 			current_level = next_level;
	// 
	// 		}
	// 		node[prev_top].parentidx = prev_top;
	// 		rootidx = prev_top;
	// 		curSize--;
	// 
	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}
	// 
	// 		delete queue;
	// 		Free(dimg);
	// 		Free(levelroot);
	// 		Free(levelroot_flag);
	// 		Free(isVisited);
	// 		Free(isAvailable);
	// 	}

	void Flood_HQueue(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		HQueue<Imgidx>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create hierarchical queue from dhist
		queue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
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
		maxSize = 2 * imgsize;
		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		queue->push(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;
				if (is_visited(isVisited, p))
				{
					queue->find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width)) push_neighbour(queue, dimg, p + width, q);//bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		push_neighbour(queue, dimg, p + 1, q + 1);//right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		push_neighbour(queue, dimg, p - 1, q - 1);//left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width)) push_neighbour(queue, dimg, p - width, q - (width << 1));//top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE
					queue->find_min_level();
#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_PQueue(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		PriorityQueue<Imgidx, Pixel>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create heap-based priority queue from dhist
		queue = new PriorityQueue<Imgidx, Pixel>(nredges);
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
		maxSize = 2 * imgsize;

		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		//queue->push(x0, current_level);
		queue->push_run(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (p == 7)
					p = p;
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;
				if (is_visited(isVisited, p))
				{
					queue->pop_run();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width)) push_neighbour(queue, dimg, p + width, q);//bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		push_neighbour(queue, dimg, p + 1, q + 1);//right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		push_neighbour(queue, dimg, p - 1, q - 1);//left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width)) push_neighbour(queue, dimg, p - width, q - (width << 1));//top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				queue->find_min_level();

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE

#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}


	void Flood_HQueue1(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		HQueue_hdr<Imgidx>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create hierarchical queue from dhist
		queue = new HQueue_hdr<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
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
		maxSize = 2 * imgsize;
		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		queue->push(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;
				if (is_visited(isVisited, p))
				{
					queue->find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width)) push_neighbour(queue, dimg, p + width, q);//bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		push_neighbour(queue, dimg, p + 1, q + 1);//right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		push_neighbour(queue, dimg, p - 1, q - 1);//left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width)) push_neighbour(queue, dimg, p - width, q - (width << 1));//top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE
					queue->find_min_level();
#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}


	void Flood_HQueue2(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		HQueue_hdr2<Imgidx>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create hierarchical queue from dhist
		queue = new HQueue_hdr2<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
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
		maxSize = 2 * imgsize;
		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		queue->push(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;
				if (is_visited(isVisited, p))
				{
					queue->find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width)) push_neighbour(queue, dimg, p + width, q);//bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		push_neighbour(queue, dimg, p + 1, q + 1);//right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		push_neighbour(queue, dimg, p - 1, q - 1);//left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width)) push_neighbour(queue, dimg, p - width, q - (width << 1));//top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE
					queue->find_min_level();
#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_HQueue3(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		HybridQueue_HQueue<Imgidx>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create hierarchical queue from dhist
		queue = new HybridQueue_HQueue<Imgidx>(nredges + 1, dhist, numlevels); // +1 for the dummy node
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
		maxSize = 2 * imgsize;
		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		queue->push(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;

				if (p == 18 && current_level == 7)
					p = p;

				if (is_visited(isVisited, p))
				{
					queue->find_min_level();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))   queue->push(p + width, (int64)dimg[q]);                //bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		  queue->push(p + 1, (int64)dimg[q + 1]);                //right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		  queue->push(p - 1, (int64)dimg[q - 1]);                //left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))   queue->push(p - width, (int64)dimg[q - (width << 1)]); //top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
// 					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
// 					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
// 					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
// 					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
// 					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
// 					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE
					queue->find_min_level();
#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}



	void Flood_PQueue1(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		int64 numlevels, max_level, current_level;// , next_level;
		PriorityQueue<Imgidx, Pixel>* queue;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx prev_top, stack_top, iNode;// , *levelroot;
		uint8 *isVisited, *isAvailable, isAv;;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = 1 << bit_depth;
		max_level = numlevels - 1;

		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
		init_isAvailable(isAvailable);

		//calculate pixel differences and make histogram
		compute_dimg(dimg, dhist, img);

		//create heap-based priority queue from dhist
		queue = new PriorityQueue<Imgidx, Pixel>(nredges);
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
		maxSize = 2 * imgsize;

		Free(dhist); //free some memory for node array

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32)); //pointers to leaf nodes for each pixel
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

#if DELAYED_NODE_ALLOC
		stack_top = NewAlphaNode(); //dummy root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
#else
		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
#endif

		current_level = max_level;
		x0 = 0; //first pixel to visit (can be any pixel)
		//queue->push(x0, current_level);
		queue->push_run(x0, current_level);

		prev_top = stack_top; //to detect redundant nodes

		//flooding start
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->pop_run();
				//cout << "visiting " << (int)p << " at level " << (int)current_level << endl;
				if (is_visited(isVisited, p))
				{
					//queue->pop_run();
					continue;
				}
				visit(isVisited, p);
#if !HQUEUE_COST_AMORTIZE
				queue->find_min_level();
#endif
				//push neighbours on queue
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !is_visited(isVisited, p + width)) push_neighbour1(queue, dimg, p + width, q);//bottom
					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))		push_neighbour1(queue, dimg, p + 1, q + 1);//right
					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))		push_neighbour1(queue, dimg, p - 1, q - 1);//left
					if (is_available(isAv, 3) && !is_visited(isVisited, p - width)) push_neighbour1(queue, dimg, p - width, q - (width << 1));//top

				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;
					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, dimg, p + wstride1, q);
					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(queue, dimg, p + width, q + 1);
					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, dimg, p + wstride0, q + 2);
					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(queue, dimg, p + 1, q + 3);
					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, dimg, p - wstride1, q - wstride_d + 4);
					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(queue, dimg, p - width, q - wstride_d + 1);
					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, dimg, p - wstride0, q - wstride_d - 2);
					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(queue, dimg, p - 1, q - 1);
				}

				//queue->find_min_level();

				if (current_level > queue->get_minlev())
				{
					Pixel pix_val = img[p];
					current_level = queue->get_minlev();
					iNode = NewAlphaNode();

					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
#if HQUEUE_COST_AMORTIZE

#endif
					if (current_level)
					{
						Pixel pix_val = img[p];
						//current_level = queue->get_minlev();
						iNode = NewAlphaNode();

						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode);
						node[iNode].parentidx = stack_top;
						parentAry[p] = iNode;
					}
					else
						//#if DELAYED_NODE_ALLOC
						connectPix2Node(p, img[p], stack_top);
					//#else
					//					connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
					//#endif
				}

			}
			//std::cout << "checking redundancy..." << std::endl;
			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
			if (node[prev_top].parentidx == stack_top &&
				node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
#if DELAYED_NODE_ALLOC
				curSize--;
#endif
			}

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				//connectnode2node
				iNode = NewAlphaNode(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else
				connectNode2Node0(stack_top, iNode);

			// #if DELAYED_NODE_ALLOC
			// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
			// #else
			// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			// #endif


			if (node[iNode].area == imgsize)
			{
				if (node[stack_top].area == imgsize)
					rootidx = stack_top;
				else
					rootidx = iNode;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;
				node[rootidx].parentidx = rootidx;
				break;
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}
		//	node[prev_top].parentidx = prev_top;
		//	rootidx = prev_top;
		//	curSize--;

	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}

		delete queue;
		Free(dimg);
		//Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}


	void Flood_PQueue_UnionbyRank(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank, next_rank;
		PriorityQueue_ubr<Imgidx> *queue;
		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank, top_rank;
		int8 nbits;

		Imgidx *dhist;

		Imgidx prev_top;
		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;

		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize)));
		//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; //
		memset(isVisited, 0, (size_t)((imgsize)));
		//		memset(isVisited_edges, 0, (size_t)((nredges)));
		init_isAvailable(isAvailable);

		//	std::cout << "before:" << rankinfo->alpha << std::endl;
	//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
		compute_dimg(rank, rankinfo, img);

		//std::cout << "after:" << rankinfo->alpha << std::endl;
//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		queue = new PriorityQueue_ubr<Imgidx>(nredges + 1);
		//queue = new Trie<Imgidx, trieidx>(nredges);
		//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

				//trie = new Trie<Imgidx, trieidx>(nredges << 1);
				//trie->push(nredges, 0);



		// 		//tree size estimation (TSE)
		//

		Pixel r = rankinfo[0].alpha, s;
		Imgidx count = 1;
		nrmsd = 0;
		for (p = 1; p < nredges; p++)
		{
			s = rankinfo[p].alpha;
			if (r != s)
			{
				nrmsd += ((double)count) * ((double)count);
				r = s;
				count = 1;
			}
			else
				count++;
		}

		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		//maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		// 		//maxSize = imgsize;
		// 		Free(dhist);
		maxSize = imgsize + nredges;
		//maxSize = nredges;
		num_node = maxSize;
		num_node_in = nredges;

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
			pNode->set(1, 0, (double)img[p], img[p], img[p]);

		nbits = ((sizeof(Pixel) << 3) - 1);
		maxpixval = ~(1 << nbits);
		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0, pix_type);
		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

		isVisited[0] = 1;
		if (connectivity == 4)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
		}
		else
		{
			queue->push(rank[1]);
			queue->push(rank[2]);
			queue->push(rank[3]);
		}

		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize);
		//isVisited_edges[current_rank] = 1;
		prev_top = current_rank;

		while (1)//(current_rank <= nredges)
		{
			while (1)//((trie->top() >> 1) <= current_rank)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
			//	if (isVisited_edges[top_rank])
			//		break;
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
#if !HQUEUE_COST_AMORTIZE
				//find_min_level();
#endif
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !isVisited[p + width])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])
						queue->push(rank[q - (width << 1)]);
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;

					if (is_available(isAv, 0) && !isVisited[p + wstride1])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + width])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + wstride0])
						queue->push(rank[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p + 1])
						queue->push(rank[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - wstride1])
						queue->push(rank[q - wstride_d + 4]);
					if (is_available(isAv, 5) && !isVisited[p - width])
						queue->push(rank[q - wstride_d + 1]);
					if (is_available(isAv, 6) && !isVisited[p - wstride0])
						queue->push(rank[q - wstride_d - 2]);
					if (is_available(isAv, 7) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
				}

				next_rank = queue->top();
				/*				if(p == 27579)
								{
									node[p].print(node);
									node[5773602].print(node);
									getchar();

									disp = 0;
								}
								*/
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);

				/*				disp = 0;
								if(p == 27579)
								{
									node[p].print(node);
									node_in[node[p].parentidx].print(node);
									getchar();
								}
				*/
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			//Redundant node removal
			if (node_in[prev_top].parentidx == current_rank &&
				node_in[prev_top].area == node_in[current_rank].area)
				current_rank = prev_top;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
					next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;

				break;
			}

			prev_top = current_rank;
			current_rank = next_rank;
		}

		Imgidx num_excessnodes = 0;
		for (pNode = node_in; pNode < node + maxSize; pNode++)
		{
			pNode->parentidx += imgsize;
			if (pNode->area == 0)
				num_excessnodes++;
		}
		//num_excessnodes += num_danglingnodes;
		node_in[next_rank].parentidx = -1;

		//ratio_danglingnodes = (double)num_danglingnodes / (double)maxSize;
		//ratio_excessnodes = (double)num_excessnodes / (double)nredges;


		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
//		Free(levelroot);
		Free(isVisited);
		//Free(isVisited_edges);
		Free(isAvailable);
	}



	// 	void Flood_PQueue(Pixel* img)
	// 	{
	// 		Imgidx imgsize, dimgsize, nredges, x0;
	// 		int32 numlevels, max_level, current_level, next_level;
	// 		PriorityQueue<Imgidx, Pixel>* queue;
	// 		Imgidx *dhist;
	// 		Pixel *dimg;
	// 		Imgidx prev_top, *levelroot;
	// 		uint8 *isVisited, *isAvailable, isAv;;
	// 		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
	// 
	// 		imgsize = width * height;
	// 		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
	// 		dimgsize = (connectivity >> 1) * width * height;
	// 		numlevels = 1 << bit_depth;
	// 
	// 		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
	// 		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
	// 		levelroot = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	// 		isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
	// 		if (connectivity == 4)
	// 			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
	// 		else
	// 			isAvailable = (uint8*)Malloc((size_t)(imgsize));
	// 		for (p = 0; p < numlevels; p++)
	// 			levelroot[p] = NULL_LEVELROOT;
	// 		memset(dhist, 0, (size_t)numlevels * sizeof(int32));
	// 		memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
	// 		init_isAvailable(isAvailable);
	// 
	// 		max_level = numlevels - 1;
	// 
	// 		compute_dimg(dimg, dhist, img);
	// 
	// 		dhist[max_level]++;
	// 		queue = new PriorityQueue<Imgidx, Pixel>(nredges);
	// 		curSize = 0;
	// 
	// 		/*
	// 		char str[10];
	// 		fstream fs;
	// 		fs.open("D:/RUG/2019/TTMA_ISMM/qstat0.dat", std::fstream::in);
	// 		int pushidx = 0;
	// 		int popidx = 0;
	// 		int testsize = 10240;
	// 		Pixel popalpha;
	// 		while (1)
	// 		{
	// 			fs.getline(str, 10);
	// 			if (str[0] == 0)
	// 				break;
	// 			int op = std::atoi(str);
	// 			fs.getline(str, 10);
	// 			int numer = std::atoi(str);
	// 
	// 			if (op == 0)
	// 				queue->push(pushidx++, numer);
	// 			else
	// 			{
	// 				//popidx = (int)queue->pop();
	// 				popalpha = queue->min_level;
	// 				queue->pop();
	// 				if (popalpha != numer)
	// 					popalpha = popalpha;
	// 			}
	// 		}
	// 		fs.close();
	// 		*/
	// 
	// 		/////////////////////////////////////////
	// 		//tree size estimation (TSE)
	// 		nrmsd = 0;
	// 		for (p = 0; p < numlevels; p++)
	// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	// 		maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		/////////////////////////////////////////
	// 	//	printf("NRMSD: %f\tEst. NTS: %f\tEst. Tree size: %d\n", nrmsd, ((A * exp(SIGMA * nrmsd) + B) + M), tree->maxSize);
	// 		//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	// 		Free(dhist);
	// 
	// 		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
	// 		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
	// 
	// #if DELAYED_NODE_ALLOC
	// 		//		levelroot[max_level + 1] = NODE_CANDIDATE;
	// 		levelroot[max_level + 1] = NewAlphaNode();
	// 		AlphaNode<Imgidx, Pixel> *pNode = node + levelroot[max_level + 1];
	// 		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
	// 		//(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	// 
	// #else
	// 		levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	// 		node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
	// #endif
	// 
	// 		current_level = max_level;
	// 		x0 = imgsize >> 1;
	// 		queue->push(x0, (Pixel)current_level);
	// 
	// 		prev_top = levelroot[max_level + 1];
	// 		while (current_level <= max_level)
	// 		{
	// 			while (queue->min_level <= current_level)
	// 			{
	// 				p = queue->pop();
	// 				if (is_visited(isVisited, p))
	// 				{
	// 					queue->find_min_level();
	// 					continue;
	// 				}
	// 				visit(isVisited, p);
	// #if !HQUEUE_COST_AMORTIZE
	// 				queue->find_min_level();
	// #endif
	// 				if (connectivity == 4)
	// 				{
	// 					isAv = get_field(isAvailable, p);
	// 					q = p << 1;
	// 
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(queue, levelroot, dimg, p + width, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))			push_neighbour(queue, levelroot, dimg, p + 1, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))			push_neighbour(queue, levelroot, dimg, p - 1, q - 1);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(queue, levelroot, dimg, p - width, q - (width << 1));
	// 				}
	// 				else
	// 				{
	// 					isAv = isAvailable[p];
	// 					q = p << 2;
	// 					if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(queue, levelroot, dimg, p + wstride1, q);
	// 					if (is_available(isAv, 1) && !is_visited(isVisited, p + width))			push_neighbour(queue, levelroot, dimg, p + width, q + 1);
	// 					if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(queue, levelroot, dimg, p + wstride0, q + 2);
	// 					if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))				push_neighbour(queue, levelroot, dimg, p + 1, q + 3);
	// 					if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(queue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
	// 					if (is_available(isAv, 5) && !is_visited(isVisited, p - width))			push_neighbour(queue, levelroot, dimg, p - width, q - wstride_d + 1);
	// 					if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(queue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
	// 					if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))				push_neighbour(queue, levelroot, dimg, p - 1, q - 1);
	// 				}
	// 
	// 				if (current_level > queue->min_level)
	// 					current_level = queue->min_level;
	// #if HQUEUE_COST_AMORTIZE
	// 				else
	// 					queue->find_min_level();
	// #endif
	// 
	// #if DELAYED_NODE_ALLOC
	// 				connectPix2Node(p, img[p], levelroot, current_level);
	// #else
	// 				connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
	// #endif
	// 
	// 			}
	// 			//std::cout << "checking redundancy..." << std::endl;
	// 			//std::cout << "node[prev_top].parentidx: " << node[prev_top].parentidx << std::endl;
	// 			//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
	// 			if (node[prev_top].parentidx == levelroot[current_level] &&
	// 				node[levelroot[current_level]].area == node[prev_top].area)
	// 			{
	// 				levelroot[current_level] = prev_top;
	// #if DELAYED_NODE_ALLOC
	// 				curSize--;
	// #endif
	// 			}
	// 
	// 			next_level = current_level + 1;
	// 			while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
	// 				next_level++;
	// 
	// #if DELAYED_NODE_ALLOC
	// 			connectNode2Node(levelroot, levelroot[current_level], next_level);
	// #else
	// 			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
	// #endif
	// 
	// 
	// 			prev_top = levelroot[current_level];
	// 			levelroot[current_level] = NULL_LEVELROOT;
	// 			current_level = next_level;
	// 
	// 		}
	// 		node[prev_top].parentidx = prev_top;
	// 		rootidx = prev_top;
	// 		curSize--;
	// 
	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
	// 			{
	// 				x0 = NewAlphaNode();
	// 				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	// 				node[x0].parentidx = parentAry[p];
	// 				parentAry[p] = x0;
	// 			}
	// 		}
	// 
	// 		delete queue;
	// 		Free(dimg);
	// 		Free(levelroot);
	// 		Free(isVisited);
	// 		Free(isAvailable);
	// 	}

	void Flood_Trie(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank, next_rank;
		Trie<Imgidx, trieidx> *queue;
		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank, top_rank;
		int8 nbits;

		Imgidx *dhist;

		Imgidx prev_top;
		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;

		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize)));
		//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; //
		memset(isVisited, 0, (size_t)((imgsize)));
		//		memset(isVisited_edges, 0, (size_t)((nredges)));
		init_isAvailable(isAvailable);

		//	std::cout << "before:" << rankinfo->alpha << std::endl;
	//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
		compute_dimg(rank, rankinfo, img);

		//std::cout << "after:" << rankinfo->alpha << std::endl;
//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		queue = new Trie<Imgidx, trieidx>(nredges);
		//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

				//trie = new Trie<Imgidx, trieidx>(nredges << 1);
				//trie->push(nredges, 0);



		// 		//tree size estimation (TSE)
		//

		Pixel r = rankinfo[0].alpha, s;
		Imgidx count = 1;
		nrmsd = 0;
		for (p = 1; p < nredges; p++)
		{
			s = rankinfo[p].alpha;
			if (r != s)
			{
				nrmsd += ((double)count) * ((double)count);
				r = s;
				count = 1;
			}
			else
				count++;
		}

		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		//maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		// 		//maxSize = imgsize;
		// 		Free(dhist);
		maxSize = imgsize + nredges;
		//maxSize = nredges;
		num_node = maxSize;
		num_node_in = nredges;

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
			pNode->set(1, 0, (double)img[p], img[p], img[p]);

		nbits = ((sizeof(Pixel) << 3) - 1);
		maxpixval = ~(1 << nbits);
		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0, pix_type);
		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

		isVisited[0] = 1;
		if (connectivity == 4)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
		}
		else
		{
			queue->push(rank[1]);
			queue->push(rank[2]);
			queue->push(rank[3]);
		}

		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize);
		//isVisited_edges[current_rank] = 1;
		prev_top = current_rank;

		while (1)//(current_rank <= nredges)
		{
			while (1)//((trie->top() >> 1) <= current_rank)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
			//	if (isVisited_edges[top_rank])
			//		break;
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
#if !HQUEUE_COST_AMORTIZE
				//find_min_level();
#endif
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !isVisited[p + width])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])
						queue->push(rank[q - (width << 1)]);
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;

					if (is_available(isAv, 0) && !isVisited[p + wstride1])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + width])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + wstride0])
						queue->push(rank[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p + 1])
						queue->push(rank[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - wstride1])
						queue->push(rank[q - wstride_d + 4]);
					if (is_available(isAv, 5) && !isVisited[p - width])
						queue->push(rank[q - wstride_d + 1]);
					if (is_available(isAv, 6) && !isVisited[p - wstride0])
						queue->push(rank[q - wstride_d - 2]);
					if (is_available(isAv, 7) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
				}

				next_rank = queue->top();
				/*				if(p == 27579)
								{
									node[p].print(node);
									node[5773602].print(node);
									getchar();

									disp = 0;
								}
								*/
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);

				/*				disp = 0;
								if(p == 27579)
								{
									node[p].print(node);
									node_in[node[p].parentidx].print(node);
									getchar();
								}
				*/
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			//Redundant node removal
			if (node_in[prev_top].parentidx == current_rank &&
				node_in[prev_top].area == node_in[current_rank].area)
				current_rank = prev_top;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
					next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;

				break;
			}

			prev_top = current_rank;
			current_rank = next_rank;
		}

		Imgidx num_excessnodes = 0;
		for (pNode = node_in; pNode < node + maxSize; pNode++)
		{
			pNode->parentidx += imgsize;
			if (pNode->area == 0)
				num_excessnodes++;
		}
		//num_excessnodes += num_danglingnodes;
		node_in[next_rank].parentidx = -1;

		//ratio_danglingnodes = (double)num_danglingnodes / (double)maxSize;
		//ratio_excessnodes = (double)num_excessnodes / (double)nredges;


		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
//		Free(levelroot);
		Free(isVisited);
		//Free(isVisited_edges);
		Free(isAvailable);
	}

	void Flood_HybridQueue(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank, next_rank;
		HybridQueue_Trie<Imgidx, trieidx> *queue;
		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank, top_rank;
		int8 nbits;

		Imgidx *dhist;

		Imgidx prev_top;
		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;

		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize)));
		//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; //
		memset(isVisited, 0, (size_t)((imgsize)));
		//		memset(isVisited_edges, 0, (size_t)((nredges)));
		init_isAvailable(isAvailable);

		//	std::cout << "before:" << rankinfo->alpha << std::endl;
	//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
		compute_dimg(rank, rankinfo, img);

		//std::cout << "after:" << rankinfo->alpha << std::endl;
//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		queue = new HybridQueue_Trie<Imgidx, trieidx>(nredges);
		//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

				//trie = new Trie<Imgidx, trieidx>(nredges << 1);
				//trie->push(nredges, 0);



		// 		//tree size estimation (TSE)
		//

		Pixel r = rankinfo[0].alpha, s;
		Imgidx count = 1;
		nrmsd = 0;
		for (p = 1; p < nredges; p++)
		{
			s = rankinfo[p].alpha;
			if (r != s)
			{
				nrmsd += ((double)count) * ((double)count);
				r = s;
				count = 1;
			}
			else
				count++;
		}

		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		//maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		// 		//maxSize = imgsize;
		// 		Free(dhist);
		maxSize = imgsize + nredges;
		//maxSize = nredges;
		num_node = maxSize;
		num_node_in = nredges;

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
			pNode->set(1, 0, (double)img[p], img[p], img[p]);

		nbits = ((sizeof(Pixel) << 3) - 1);
		maxpixval = ~(1 << nbits);
		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0, pix_type);
		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

		isVisited[0] = 1;
		if (connectivity == 4)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
		}
		else
		{
			queue->push(rank[1]);
			queue->push(rank[2]);
			queue->push(rank[3]);
		}

		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize);
		//isVisited_edges[current_rank] = 1;
		prev_top = current_rank;

		while (1)//(current_rank <= nredges)
		{
			while (1)//((trie->top() >> 1) <= current_rank)
			{
				top_rank = queue->top();	//remove tmp variables later if possible
			//	if (isVisited_edges[top_rank])
			//		break;
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
#if !HQUEUE_COST_AMORTIZE
				//find_min_level();
#endif
				if (connectivity == 4)
				{
					isAv = isAvailable[p];
					q = p << 1;

					if (is_available(isAv, 0) && !isVisited[p + width])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])
						queue->push(rank[q - (width << 1)]);
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;

					if (is_available(isAv, 0) && !isVisited[p + wstride1])
						queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + width])
						queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + wstride0])
						queue->push(rank[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p + 1])
						queue->push(rank[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - wstride1])
						queue->push(rank[q - wstride_d + 4]);
					if (is_available(isAv, 5) && !isVisited[p - width])
						queue->push(rank[q - wstride_d + 1]);
					if (is_available(isAv, 6) && !isVisited[p - wstride0])
						queue->push(rank[q - wstride_d - 2]);
					if (is_available(isAv, 7) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
				}

				next_rank = queue->top();
				/*				if(p == 27579)
								{
									node[p].print(node);
									node[5773602].print(node);
									getchar();

									disp = 0;
								}
								*/
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize);

				/*				disp = 0;
								if(p == 27579)
								{
									node[p].print(node);
									node_in[node[p].parentidx].print(node);
									getchar();
								}
				*/
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			//Redundant node removal
			if (node_in[prev_top].parentidx == current_rank &&
				node_in[prev_top].area == node_in[current_rank].area)
				current_rank = prev_top;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
					next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;

				break;
			}

			prev_top = current_rank;
			current_rank = next_rank;
		}

		Imgidx num_excessnodes = 0;
		for (pNode = node_in; pNode < node + maxSize; pNode++)
		{
			pNode->parentidx += imgsize;
			if (pNode->area == 0)
				num_excessnodes++;
		}
		//num_excessnodes += num_danglingnodes;
		node_in[next_rank].parentidx = -1;

		//ratio_danglingnodes = (double)num_danglingnodes / (double)maxSize;
		//ratio_excessnodes = (double)num_excessnodes / (double)nredges;


		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
//		Free(levelroot);
		Free(isVisited);
		//Free(isVisited_edges);
		Free(isAvailable);
	}
	// 
	// 	void par_Flood_HybridQueue(Pixel* img, int8 listsize)
	// 	{
	// 		Imgidx imgsize, dimgsize, nredges;
	// 		Imgidx current_rank, next_rank;
	// 		HybridQueue<Imgidx, trieidx> *queue;
	// 		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
	// 		AlphaNode<Imgidx, Pixel> *pNode;
	// 		Pixel maxpixval;
	// 		Imgidx *rank, top_rank;
	// 		int8 nbits;
	// 
	// 		Imgidx prev_top;
	// 		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
	// 		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
	// 
	// 		imgsize = width * height;
	// 		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
	// 		dimgsize = (connectivity >> 1) * width * height;
	// 
	// 		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
	// 		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
	// 		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
	// 		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
	// 		isVisited = (uint8*)Malloc((size_t)((imgsize)));
	// 		//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
	// 		if (connectivity == 4)
	// 			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
	// 		else
	// 			isAvailable = (uint8*)Malloc((size_t)(imgsize));
	// 
	// 		//	levelroot[p] = NULL_LEVELROOT;
	// 		//levelroot[nredges] = NODE_CANDIDATE; //
	// 		memset(isVisited, 0, (size_t)((imgsize)));
	// 		//		memset(isVisited_edges, 0, (size_t)((nredges)));
	// 		init_isAvailable_par(isAvailable, 2); //yogi
	// 
	// 	//	std::cout << "before:" << rankinfo->alpha << std::endl;
	// 	//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
	// 		compute_dimg(rank, rankinfo, img);
	// 		//std::cout << "after:" << rankinfo->alpha << std::endl;
	// 	//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
	// 
	// 		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);
	// 		//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);
	// 
	// 				//trie = new Trie<Imgidx, trieidx>(nredges << 1);
	// 				//trie->push(nredges, 0);
	// 
	// 
	// 
	// 		// 		//tree size estimation (TSE)
	// 		//
	// 		// 		nrmsd = 0;
	// 		// 		for (p = 0; p < numlevels; p++)
	// 		// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	// 		// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	// 		// 		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
	// 		// 		//maxSize = imgsize;
	// 		// 		Free(dhist);
	// 		maxSize = imgsize + nredges;
	// 		num_node = maxSize;
	// 		num_node_in = nredges;
	// 
	// 		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
	// 		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
	// 		node_in = node + imgsize;
	// 
	// 		std::cout << "initiating nodes" << std::endl;
	// 		//		pNode = node;pNode++
	// #pragma omp parallel for
	// 		for (p = 0; p < imgsize; p++)
	// 		{
	// 			node[p].set(1, 0, (double)img[p], img[p], img[p]);
	// 			//node[p].print(node);
	// 		}
	// 
	// 
	// 		pNode = node_in;
	// 		nbits = ((sizeof(Pixel) << 3) - 1);
	// 		maxpixval = ~(1 << nbits);
	// 		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
	// 			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0);
	// 		//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	// 
	// 		isVisited[0] = 1;
	// 		if (connectivity == 4)
	// 		{
	// 			queue->push(rank[0]);
	// 			queue->push(rank[1]);
	// 		}
	// 		else
	// 		{
	// 			queue->push(rank[1]);
	// 			queue->push(rank[2]);
	// 			queue->push(rank[3]);
	// 		}
	// 
	// 		current_rank = queue->top();
	// 		node[0].connect_to_parent(&node_in[current_rank], current_rank);
	// 		//isVisited_edges[current_rank] = 1;
	// 		prev_top = current_rank;
	// 
	// 		while (1)//(current_rank <= nredges)
	// 		{
	// 			while (1)//((trie->top() >> 1) <= current_rank)
	// 			{
	// 				top_rank = queue->top();	//remove tmp variables later if possible
	// 			//	if (isVisited_edges[top_rank])
	// 			//		break;
	// 				pRank = rankinfo + top_rank;
	// 				if (isVisited[pRank->p])
	// 				{
	// 					if (isVisited[pRank->q])
	// 						break;
	// 					p = pRank->q;
	// 				}
	// 				else
	// 					p = pRank->p;
	// 
	// 				isVisited[p] = 1;
	// #if !HQUEUE_COST_AMORTIZE
	// 				//find_min_level();
	// #endif
	// 				if (connectivity == 4)
	// 				{
	// 					isAv = get_field(isAvailable, p);
	// 					q = p << 1;
	// 
	// 					if (is_available(isAv, 0) && !isVisited[p + width])
	// 						queue->push(rank[q]);
	// 					if (is_available(isAv, 1) && !isVisited[p + 1])
	// 						queue->push(rank[q + 1]);
	// 					if (is_available(isAv, 2) && !isVisited[p - 1])
	// 						queue->push(rank[q - 1]);
	// 					if (is_available(isAv, 3) && !isVisited[p - width])
	// 						queue->push(rank[q - (width << 1)]);
	// 				}
	// 				else
	// 				{
	// 					isAv = isAvailable[p];
	// 					q = p << 2;
	// 
	// 					if (is_available(isAv, 0) && !isVisited[p + wstride1])
	// 						queue->push(rank[q]);
	// 					if (is_available(isAv, 1) && !isVisited[p + width])
	// 						queue->push(rank[q + 1]);
	// 					if (is_available(isAv, 2) && !isVisited[p + wstride0])
	// 						queue->push(rank[q + 2]);
	// 					if (is_available(isAv, 3) && !isVisited[p + 1])
	// 						queue->push(rank[q + 3]);
	// 					if (is_available(isAv, 4) && !isVisited[p - wstride1])
	// 						queue->push(rank[q - wstride_d + 4]);
	// 					if (is_available(isAv, 5) && !isVisited[p - width])
	// 						queue->push(rank[q - wstride_d + 1]);
	// 					if (is_available(isAv, 6) && !isVisited[p - wstride0])
	// 						queue->push(rank[q - wstride_d - 2]);
	// 					if (is_available(isAv, 7) && !isVisited[p - 1])
	// 						queue->push(rank[q - 1]);
	// 				}
	// 
	// 				next_rank = queue->top();
	// 				/*				if(p == 27579)
	// 							{
	// 								node[p].print(node);
	// 								node[5773602].print(node);
	// 								getchar();
	// 
	// 								disp = 0;
	// 							}
	// 							*/
	// 				node[p].connect_to_parent(&node_in[next_rank], next_rank);
	// 
	// 				/*				disp = 0;
	// 							if(p == 27579)
	// 							{
	// 								node[p].print(node);
	// 								node_in[node[p].parentidx].print(node);
	// 								getchar();
	// 							}
	// 				*/
	// 				if (current_rank == next_rank)
	// 					break;
	// 				current_rank = next_rank;
	// 			}
	// 
	// 			queue->pop();
	// 			next_rank = queue->top();
	// 
	// 			//Redundant node removal
	// 			if (node_in[prev_top].parentidx == current_rank &&
	// 				node_in[prev_top].area == node_in[current_rank].area)
	// 				current_rank = prev_top;
	// 
	// 			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
	// 			if (node_in[next_rank].area == imgsize)
	// 			{
	// 				if (node_in[current_rank].area == imgsize)
	// 					next_rank = current_rank;
	// 				//				node_in[next_rank].parentidx = next_rank;
	// 				// node_in[next_rank].parentidx = -1;
	// 
	// 				break;
	// 			}
	// 
	// 			prev_top = current_rank;
	// 			current_rank = next_rank;
	// 		}
	// 
	// 		for (pNode = node; pNode < node + maxSize; pNode++)
	// 			pNode->parentidx += imgsize;
	// 		node_in[next_rank].parentidx = -1;
	// 
	// 		delete queue;
	// 		Free(rank);
	// 		Free(rankinfo);
	// 		//Free(rank2alpha);
	// 	//		Free(levelroot);
	// 		Free(isVisited);
	// 		//Free(isVisited_edges);
	// 		Free(isAvailable);
	// 	}
	Imgidx node_at_minlev;
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
		case(0): //HQueue
			Flood_HQueue(img);
			break;
		case(1): //PQueue
			Flood_PQueue(img);
			break;
		case(2): //Trie
			Flood_Trie(img);
			break;
		case(3): //HybridQueue
			Flood_HybridQueue(img);
			break;
		case(4):
			Flood_PQueue_UnionbyRank(img);
			break;
		case(5):
			Flood_HQueue1(img);
			break;
		case(6):
			Flood_HQueue2(img);
			break;
		default:
			break;
		}
		/*
		if (sizeof(Pixel) == 1)
		{
			//tmp
			//std::cout << "buildalphatree - hqueue" << std::endl;
			Flood_HQueue(img);
			//Flood_PQueue(img);
			//Flood_Trie(img);
			//Flood_HybridQueue(img);
		}
		else
		{
			Flood_HQueue(img);
			//Flood_PQueue(img);
			//Flood_HybridQueue(img);
			//Flood_Trie(img);
			//Flood_HybridQueue(img);
		}
		*/
	}

	void AreaFilter(Pixel *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		//	Pixel val;
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
			while (pNode->parentidx != -1 &&
				pNode->alpha < alpha
				//&& pNode->alpha < alpha
				)
			{
				//					if(val == 9260)
				{
					//						pNode->print(node);
				}
				pNode = &node[pNode->parentidx];
			}
			// 				if(pNode->filter_val == 0 && val < 200)
			//					pNode->filter_val = val++;

					//	printf("%lld\n",(uint64)val);
			// 				outimg[i] = ((double)pNode->area / (double)imgsize) * 255;
			outimg[i] = min(255, (Pixel)pNode->area / 200.0);
		}
		//	printf("%lld\n",(uint64)val);
	}

	void AreaFilter(double *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		//	Pixel val;
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
			while (pNode->parentidx != -1 &&
				pNode->alpha < alpha
				//&& pNode->alpha < alpha
				)
			{
				//					if(val == 9260)
				{
					//						pNode->print(node);
				}
				pNode = &node[pNode->parentidx];
			}
			// 				if(pNode->filter_val == 0 && val < 200)
			//					pNode->filter_val = val++;

						//	printf("%lld\n",(uint64)val);
			// 				outimg[i] = ((double)pNode->area / (double)imgsize) * 255;
			outimg[i] = min(255, (double)pNode->area / 200.0);
		}
		//	printf("%lld\n",(uint64)val);
	}

	void print_tree()
	{
		for (int i = 0; i < maxSize; i++)
		{
			node[i].print(node);
		}
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
		if (!tree)
			return;
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)
				delete ((ATree<int32, uint8>*)tree);
			else if (pix_type == PIXEL_16BIT)
				delete ((ATree<int32, uint16>*)tree);
			else if (pix_type == PIXEL_32BIT)
				delete ((ATree<int32, uint32>*)tree);
			else if (pix_type == PIXEL_64BIT)
				delete ((ATree<int32, uint64>*)tree);
			else if (pix_type == PIXEL_FLOAT)
				delete ((ATree<int32, uint32>*)tree);
			else
				delete ((ATree<int32, uint64>*)tree);
		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				delete ((ATree<int64, uint8>*)tree);
			else if (pix_type == PIXEL_16BIT)
				delete ((ATree<int64, uint16>*)tree);
			else if (pix_type == PIXEL_32BIT)
				delete ((ATree<int64, uint32>*)tree);
			else if (pix_type == PIXEL_64BIT)
				delete ((ATree<int64, uint64>*)tree);
			else if (pix_type == PIXEL_FLOAT)
				delete ((ATree<int64, uint32>*)tree);
			else
				delete ((ATree<int64, uint64>*)tree);
		}
	}

	inline void clear()
	{
		/*		if(tree)
				{
					switch(imgidx)
					{
						case (IMGIDX_32BIT):
						case (IMGIDX_64BIT):
					}
					}
					tree->\
				}
			*/
	}

	void BuildAlphaTree(uint8 *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pix_type = PIXEL_8BIT;
		if ((int64)height * (int64)width < (int64)0xefffffff)
		{
			//tmp
			//std::cout << "Entered wrapper build" << std::endl;

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
			if (pix_type == PIXEL_8BIT)
				((ATree<int32, uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
			else if (pix_type == PIXEL_16BIT)
				((ATree<int32, uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
			else if (pix_type == PIXEL_32BIT)
				((ATree<int32, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else if (pix_type == PIXEL_64BIT)
				((ATree<int32, uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
			else if (pix_type == PIXEL_FLOAT)
				((ATree<int32, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else
				((ATree<int32, uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);

		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				((ATree<int64, uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
			else if (pix_type == PIXEL_16BIT)
				((ATree<int64, uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
			else if (pix_type == PIXEL_32BIT)
				((ATree<int64, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else if (pix_type == PIXEL_64BIT)
				((ATree<int64, uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
			else if (pix_type == PIXEL_FLOAT)
				((ATree<int64, uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
			else
				((ATree<int64, uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);
		}

	}
	void print_tree()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)
				((ATree<int32, uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)
				((ATree<int32, uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)
				((ATree<int32, uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)
				((ATree<int32, uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)
				((ATree<int32, uint32>*)tree)->print_tree();
			else
				((ATree<int32, uint64>*)tree)->print_tree();
		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				((ATree<int64, uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)
				((ATree<int64, uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)
				((ATree<int64, uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)
				((ATree<int64, uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)
				((ATree<int64, uint32>*)tree)->print_tree();
			else
				((ATree<int64, uint64>*)tree)->print_tree();
		}
	}

	int64 get_maxSize()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)
				return (int64)((ATree<int32, uint8>*)tree)->maxSize;
			else if (pix_type == PIXEL_16BIT)
				return (int64)((ATree<int32, uint16>*)tree)->maxSize;
			else if (pix_type == PIXEL_32BIT)
				return (int64)((ATree<int32, uint32>*)tree)->maxSize;
			else if (pix_type == PIXEL_64BIT)
				return (int64)((ATree<int32, uint64>*)tree)->maxSize;
			else if (pix_type == PIXEL_FLOAT)
				return (int64)((ATree<int32, uint32>*)tree)->maxSize;
			else
				return (int64)((ATree<int32, uint64>*)tree)->maxSize;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				return (int64)((ATree<int64, uint8>*)tree)->maxSize;
			else if (pix_type == PIXEL_16BIT)
				return (int64)((ATree<int64, uint16>*)tree)->maxSize;
			else if (pix_type == PIXEL_32BIT)
				return (int64)((ATree<int64, uint32>*)tree)->maxSize;
			else if (pix_type == PIXEL_64BIT)
				return (int64)((ATree<int64, uint64>*)tree)->maxSize;
			else if (pix_type == PIXEL_FLOAT)
				return (int64)((ATree<int64, uint32>*)tree)->maxSize;
			else
				return (int64)((ATree<int64, uint64>*)tree)->maxSize;
		}
	}
	int64 get_curSize()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)
				return (int64)((ATree<int32, uint8>*)tree)->curSize;
			else if (pix_type == PIXEL_16BIT)
				return (int64)((ATree<int32, uint16>*)tree)->curSize;
			else if (pix_type == PIXEL_32BIT)
				return (int64)((ATree<int32, uint32>*)tree)->curSize;
			else if (pix_type == PIXEL_64BIT)
				return (int64)((ATree<int32, uint64>*)tree)->curSize;
			else if (pix_type == PIXEL_FLOAT)
				return (int64)((ATree<int32, uint32>*)tree)->curSize;
			else
				return (int64)((ATree<int32, uint64>*)tree)->curSize;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				return (int64)((ATree<int64, uint8>*)tree)->curSize;
			else if (pix_type == PIXEL_16BIT)
				return (int64)((ATree<int64, uint16>*)tree)->curSize;
			else if (pix_type == PIXEL_32BIT)
				return (int64)((ATree<int64, uint32>*)tree)->curSize;
			else if (pix_type == PIXEL_64BIT)
				return (int64)((ATree<int64, uint64>*)tree)->curSize;
			else if (pix_type == PIXEL_FLOAT)
				return (int64)((ATree<int64, uint32>*)tree)->curSize;
			else
				return (int64)((ATree<int64, uint64>*)tree)->curSize;
		}
	}

	double get_nrmsd()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)
				return (double)((ATree<int32, uint8>*)tree)->nrmsd;
			else if (pix_type == PIXEL_16BIT)
				return (double)((ATree<int32, uint16>*)tree)->nrmsd;
			else if (pix_type == PIXEL_32BIT)
				return (double)((ATree<int32, uint32>*)tree)->nrmsd;
			else if (pix_type == PIXEL_64BIT)
				return (double)((ATree<int32, uint64>*)tree)->nrmsd;
			else if (pix_type == PIXEL_FLOAT)
				return (double)((ATree<int32, uint32>*)tree)->nrmsd;
			else
				return (double)((ATree<int32, uint64>*)tree)->nrmsd;
		}
		else
		{
			if (pix_type == PIXEL_8BIT)
				return (double)((ATree<int64, uint8>*)tree)->nrmsd;
			else if (pix_type == PIXEL_16BIT)
				return (double)((ATree<int64, uint16>*)tree)->nrmsd;
			else if (pix_type == PIXEL_32BIT)
				return (double)((ATree<int64, uint32>*)tree)->nrmsd;
			else if (pix_type == PIXEL_64BIT)
				return (double)((ATree<int64, uint64>*)tree)->nrmsd;
			else if (pix_type == PIXEL_FLOAT)
				return (double)((ATree<int64, uint32>*)tree)->nrmsd;
			else
				return (double)((ATree<int64, uint64>*)tree)->nrmsd;
		}
	}
	/*
	void BuildAlphaTree(float *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pixel = PIXEL_8BIT;
		if((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = Malloc(sizeof(ATree<int32,float>));
			((ATree<int32,float>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = Malloc(sizeof(ATree<int64,float>));
			((ATree<int64,float>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}
	void BuildAlphaTree(double *img, int height, int width, int channel, int connectivity, int algorithm)
	{
		pixel = PIXEL_8BIT;
		if((int64)height * (int64)width < (int64)0xefffffff)
		{
			imgidx = IMGIDX_32BITS;
			tree = Malloc(sizeof(ATree<int32,double>));
			((ATree<int32,double>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
		else
		{
			imgidx = IMGIDX_64BITS;
			tree = Malloc(sizeof(ATree<int64,double>));
			((ATree<int64,double>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, algorithm);
		}
	}
	*/
};
