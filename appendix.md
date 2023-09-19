# Appendix
There is a bunch of wacky stuff in this library, so here is what I found out.

Something I didn't realize before is that the BLAS functions actually follow
the conventions of standard BLAS routines.

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

### `?Insert` (eg `ipqInsert`)

Key-value priority queue insertion. Created by `gk_mkpqueue.h` (nvim doesn't
like the file).

```c
    void ipqInsert(queue_t *queue, node_t node, value_t val);
```

for `ipq`, (index priority queue, as opposed to real priority queue), `node_t`
and `value_t` are both `idx_t`.

In addition to inserting into the priority queue, it seems to also keep a map
of the values. This basically implies that the max value is always less than or
equal to the number of items in the queue. It's also a fixed size - there's no
reallocing that goes on.

In rust terms, here is approximately what the construction is:
```rust
struct MetisQueue<K, V>
where
K: Ord,
[K]: Index<V>
{ /* ... */ }
```

### `?Update` (eg `ipqUpdate`)

Similar to `?Insert`

```c
ipqUpdate(queue_t *queue, value_t node, key_t newkey);
```

It updates the key (which only determines the pop order) and moves it if
applicable.
