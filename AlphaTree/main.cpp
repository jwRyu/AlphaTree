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

#include "AlphaTree.hpp"
#include "HQueue.hpp"
#include "allocator.h"
#include "Trie.hpp"
using namespace std;

#define OUTPUT_FNAME "C:/Users/jwryu/RUG/2018/AlphaTree/test.dat"

#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Grey"
#define INPUTIMAGE_DIR_COLOUR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Colour" //colour images are used after rgb2grey conversion
#define REPEAT 10
#define RUN_TSE_ONLY 0

#define DEBUG 0

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

void RandomizedHDRimage(uint64* hdrimg, uint8* ldrimg, int64 imgsize)
{
	uint64 pix;

	for (int64 i = 0; i < imgsize; i++)
	{
		pix = ((uint64)ldrimg[i]) << 56;
		pix |= ((uint64)(rand() & 0xff) << 48);
		pix |= ((uint64)(rand() & 0xff) << 40);
		pix |= ((uint64)(rand() & 0xff) << 32);
		pix |= ((uint64)(rand() & 0xff) << 24);
		pix |= ((uint64)(rand() & 0xff) << 16);
		pix |= ((uint64)(rand() & 0xff) << 8);
		pix |= ((uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
}

void RandomizedHDRimage(uint32* hdrimg, uint8* ldrimg, int64 imgsize)
{
	uint32 pix;

	for (int64 i = 0; i < imgsize; i++)
	{
		pix = ((uint64)ldrimg[i]) << 24;
		pix |= ((uint64)(rand() & 0xff) << 16);
		pix |= ((uint64)(rand() & 0xff) << 8);
		pix |= ((uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
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
	AlphaTree<int32, uint8> *tree;
	AlphaTree<int32, uint64> *hdrtree;
	int32 width, height, channel;
	int32 cnt = 0;
	ofstream f;
	ifstream fcheck;
	char in;
	int32 i, contidx;
	std::string path;
	double time_elapsed = 0;
	double pixels_processed = 0;
	uint8 testimg[4*4] = {4, 7, 4, 0, 4, 1, 1, 0, 0, 3, 0, 0, 2, 2, 0, 5};

	srand(time(NULL));
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

			//uint64 *hdrimg = new uint64[width * height * channel];
			//RandomizedHDRimage(hdrimg, testimg, 16);

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

				//hdrtree = (AlphaTree<int32, uint64>*)Malloc(sizeof(AlphaTree<int32, uint64>));
				//hdrtree->BuildAlphaTree(hdrimg, 4, 4, channel, 4);

				tree = (AlphaTree<int32, uint8>*)Malloc(sizeof(AlphaTree<int32, uint8>));
				tree->BuildAlphaTree(cvimg.data, height, width, channel, 4);
				//tree->BuildAlphaTree(testimg, 4, 4, channel, 4);
				std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
				runtime = wctduration.count();
				minruntime = testrep == 0 ? runtime : min(runtime, minruntime);

				if (testrep < (REPEAT - 1))
					tree->clear();
			}
			f << p.path().string().c_str() << '\t' << height << '\t' << width << '\t' << max_memuse << '\t' << tree->nrmsd << '\t' << tree->maxSize << '\t' << tree->curSize << '\t' << minruntime << i << endl;

			pixels_processed += width * height;
			time_elapsed += minruntime;
			cout << "Time Elapsed: " << minruntime << "# Nodes: " << tree->curSize << " mean processing speed(Mpix/s): " << pixels_processed / (time_elapsed * 1000000) << endl;
			cvimg.release();
			str1.clear();
			tree->clear();
			//return 0;
		}
	}

	f.close();
	Free(tree);
	return 0;
}