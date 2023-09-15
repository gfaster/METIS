# Appendix
There is a bunch of wacky stuff in this library, so here is what I found out.

## functions

### `?set` (eg `iset`)

just a simple memset. Basically equiv to:
```c
int* iset(size_t n, int v, int *p) {
    size_t i;
    for (i = 0; i < n; i++) {
        p[i] = v;
    }
    return p;
}
```


### `?copy` (eg `icopy`)

just a simple memmove, but with the operands swapped
I might be able to replace these with memcpy which may be faster
```c
    Tcopy(n, src, dst)

    // equals

    memmove(dst, src, n * sizeof(T))
```
