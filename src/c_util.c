// Utilities for interfacing with C from Rust, primarily variadics


#include "metislib.h"

void 
gk_free_one(void **ptr)
{
  gk_free(ptr, LTERM);
}
