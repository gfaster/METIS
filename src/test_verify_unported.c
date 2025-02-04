#include <stdint.h>
#include <stddef.h>

// as rename.h does
#define TestVerify			libmetis__TestVerify

uint32_t TestVerify(void);

uint32_t call_test_verify_from_c(void) {
  return TestVerify();
}
