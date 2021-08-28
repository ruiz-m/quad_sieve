#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include"heap.h"

heap *new_heap(uint32_t size)
{
	heap *h = (heap*)malloc(sizeof(heap));
	h->count = 0;
	h->queue = (hnode*)malloc(size*sizeof(hnode));
	return h;
}

void reheapify(heap *h, uint32_t index)
{
	uint32_t c1 = 2*index + 1;
	uint32_t c2 = 2*index + 2;

	uint32_t smallest = index;
	int64_t x = h->queue[index].x;
	uint64_t absx = abs(x);
	uint64_t absc1x = abs(h->queue[c1].x);
	uint64_t absc2x = abs(h->queue[c2].x);

	if(c1 < h->count && absx > absc1x)
	{
		smallest = c1;
		absx = absc1x;
	}
	if(c2 < h->count && absx > absc2x)
	{
		smallest = c2;
		absx = absc2x;
	}
	if(smallest != index)
	{
		int64_t x = h->queue[smallest].x;
		uint32_t qInd = h->queue[smallest].qInd;

		h->queue[smallest].x = h->queue[index].x;
		h->queue[smallest].qInd = h->queue[index].qInd;

		h->queue[index].x = x;
		h->queue[index].qInd = qInd;
		
		reheapify(h, smallest);
	}
}

void hInsert(heap *h, int64_t x, uint32_t qInd)
{
	int32_t index = h->count;
	int32_t parent = (h->count-1) >> 1;

	h->queue[index].x = x;
	h->queue[index].qInd = qInd;
	
	while(parent < h->count && abs(h->queue[parent].x) > abs(x))
	{
		int64_t x = h->queue[parent].x;
		uint32_t qInd = h->queue[parent].qInd;

		h->queue[parent].x = h->queue[index].x;
		h->queue[parent].qInd = h->queue[index].qInd;

		h->queue[index].x = x;
		h->queue[index].qInd = qInd;

		index = parent;
		parent = (parent-1) >> 1;
	}
	++h->count;
}

void hDelete(heap *h)
{
	h->queue[0].x = h->queue[h->count-1].x;
	h->queue[0].qInd = h->queue[h->count-1].qInd;

	--h->count;

	reheapify(h, 0);
}

int64_t peekX(heap *h)
{
	return h->queue[0].x;
}

uint32_t peekQind(heap *h)
{
	return h->queue[0].qInd;
}

void hPrint(heap *h)
{
	for(int i=0; i<h->count; ++i)
	{
		printf("[%lld, %u], ", h->queue[i].x, h->queue[i].qInd);
	}
	printf("\n");
}

void hFree(heap *h)
{
	free(h->queue);
	free(h);
}
