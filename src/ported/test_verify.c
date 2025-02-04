#include <stdint.h>
#include <stddef.h>
#include "ifunc.h"

// as rename.h does
#define TestVerify			libmetis__TestVerify

uint32_t TestVerify(void);

IFUNC(uint32_t, TestVerify, (void));
uint32_t c__libmetis__TestVerify(void) {
  return 0xabad1dea;
}


uint32_t call_test_verify_from_c(void) {
  return TestVerify();
}
