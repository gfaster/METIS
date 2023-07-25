# Appendix
There is a bunch of wacky stuff in this library, so here is what I found out.

## functions

### `?set` (eg `iset`)

just a simple memset. Basically equiv to:
```c
void iset(size_t n, int v, int *p) {
    size_t i;
    for (i = 0; i < n; i++) {
        p[i] = v;
    }
    return p;
}
```
where different prefixes replace uses of `int` with something else.
