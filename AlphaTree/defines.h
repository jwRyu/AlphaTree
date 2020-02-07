#pragma once

#define DELAYED_NODE_ALLOC		1
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define NODE_CANDIDATE		0xfffffffe
// 
// #define dimg_idx_v(pidx) ((pidx)<<1)
// #define dimg_idx_h(pidx) ((pidx)<<1)+1
// 
// #define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
// #define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
// #define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
// #define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

#define HIERARCHICAL_QUEUE			0
#define HIERARCHICAL_L1IDX_QUEUE	1
#define HIERARCHICAL_L2IDX_QUEUE	2
#define HEAP_QUEUE					3
#define HEAP_RANK_QUEUE				4
#define TRIE_QUEUE					5
#define TRIE_HYBRID_QUEUE			6

#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)>(b)?(b):(a)

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long uint32;
typedef unsigned long long uint64;
typedef char int8;
typedef short int16;
typedef long int32;
typedef long long int64;

typedef int64 trieidx;