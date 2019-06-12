#pragma once

extern size_t memuse, max_memuse;

void Preallocate(size_t size);
void* Malloc(size_t size);
void Free();
//void* Realloc(void* ptr, size_t size);
//void Free(void* ptr);