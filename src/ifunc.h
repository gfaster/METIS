// Ifunc implementations

#ifndef RS_LIBMETIS_IFUNC_
#define RS_LIBMETIS_IFUNC_

#define IFUNC(ret, name, ...) \
extern void * resolve_##name(void); \
static void * resolve_inner_##name(void) { \
        return resolve_##name(); \
} \
ret __attribute__((ifunc("resolve_inner_"#name))) name __VA_ARGS__;
// ret c__libmetis__##name __VA_ARGS__


#endif // RS_LIBMETIS_IFUNC_
