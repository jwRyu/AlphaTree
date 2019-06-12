#include <cstdlib>
#include "defines.h"

size_t memuse, max_memuse, memary_size;

void *memary;

void Preallocate(size_t size)
{
	memary_size = size;
	memary = calloc(size, 1);
	memuse = 0;
	//max_memuse = 0;
}

void* Malloc(size_t size)
{
	void* cur = (void*)((char*)memary + memuse);
	memuse += size;
	//max_memuse = max(memuse, max_memuse);
	return cur;
}

void Free()
{
	free(memary);
}


/*

void* Malloc(size_t size)
{
	void* pNew = malloc(size + sizeof(size_t));

	memuse += size;
	max_memuse = max(memuse, max_memuse);

	*((size_t*)pNew) = size;
	return (void*)((size_t*)pNew + 1);
}

void* Realloc(void* ptr, size_t size)
{
	void* pOld = (void*)((size_t*)ptr - 1);
	size_t oldsize = *((size_t*)pOld);
	void* pNew = realloc(pOld, size + sizeof(size_t));

	if (pOld != pNew)
		max_memuse = max(memuse + size, max_memuse);
	else
		max_memuse = max(memuse + size - oldsize, max_memuse);
	memuse += size - oldsize;

	*((size_t*)pNew) = size;
	return (void*)((size_t*)pNew + 1);
}

void Free(void* ptr)
{
	size_t size = *((size_t*)ptr - 1);
	memuse -= size;
	free((void*)((size_t*)ptr - 1));
}
*/