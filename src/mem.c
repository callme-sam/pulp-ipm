#ifdef TARGET_EXECUTION

#include "mem.h"

void *my_malloc(size_t size)
{
    return pmsis_l2_malloc(size);
}

void my_free(void *ptr)
{
    return;
}

#endif  /* TARGET_EXECUTION */