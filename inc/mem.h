#ifdef TARGET_EXECUTION
#ifndef MEM_H_
#define MEM_H_

#include <stddef.h>

#include <pmsis.h>

#define malloc(size) my_malloc(size)
#define free(ptr) my_free(ptr)

void *my_malloc(size_t size);
void my_free(void *ptr);

#endif  /* MEM_H_ */
#endif  /* TARGET_EXECUTION */