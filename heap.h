#ifndef __HEAP_H__
#define __HEAP_H__
#include<stdint.h>

typedef struct node
{
	int64_t x;
	uint32_t qInd;
}hnode;

typedef struct heap
{
	uint32_t count;
	hnode *queue;
}heap;

heap *new_heap(uint32_t size);
void hInsert(heap *h, int64_t x, uint32_t qInd);
void hDelete(heap *h);
int64_t peekX(heap *h);
uint32_t peekQind(heap *h);
void hPrint();
void hFree(heap *h);


#endif //__HEAP_H__
