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

## Direct Access List
the format that the boundary list is stored in, changed by the `BNDInsert` and
`BNDDelete` macros.

This direct access list is really a way of storing a bounded set of integers
with `O(1)` insert, delete, and access, but `O(U)` memory. It's stored as two
slices, `lptr` (`bndptr`) and `lind` (`bndind`), along with an integer to keep
track of the size.

`lptr` keeps track of the location of elements by index, whereas `lind`
stores the set contiguously. 

Inserting an integer `x` into the set is easy: we push `x` to the end of `lind`
and set `lptr[x]` to the size variable, which is then incremented.

ex:
```
lind: [ 4,  0,  1,  0,  0]; // since size == 3 last 2 elements not in set
lptr: [ 1,  2, -1, -1,  0]; // -1 means element isn't in set
size = 3;
// contents = {4, 0, 1}

// insert 2
lind: [ 4,  0,  1,  2,  0]; // appended to the end
lptr: [ 1,  2,  3, -1,  0]; // lptr[2] set to position where 2 was appended
size = 4;
// contents = {4, 0, 1, 2}
```

Deletion is a swap-remove (see `Vec::swap_remove`), but we need to update the
swapped value's `lptr` entry to reflect its new position.

ex:
```
lind: [ 4,  0,  1,  0,  0]; // since size == 3 last 2 elements not in set
lptr: [ 1,  2, -1, -1,  0]; // -1 means element isn't in set
size = 3;
// contents = {4, 0, 1}

// remove 4
// 4's spot in lind was taken by the one at the end
lind: [ 1,  0,  1,  0,  0]; // now only first two elements are in set
lptr: [ 1,  0, -1, -1, -1]; // set lptr[4] = -1 and lptr[1] = 0 (4's old pos)
size = 2;
```

### `gk_malloc` and `?wspacemalloc`

`wspacemalloc` has scope defined by `WCOREPUSH` and `WCOREPOP`, where
`WCOREPOP` frees all `wspacemalloc` calls made since the last `WCOREPUSH`.
Between these calls, `wspacemalloc` is a bump allocator as long as the core has
space, otherwise it defers to `gk_malloc`. When `AllocateWorkspace` is called,
the core used by workspace allocations is created with a best-estimate size.

`gk_malloc` somewhat works similarly but with scope defined by `gk_malloc_init`
and `gk_malloc_cleanup`. However, unlike the workspace allocations, they can be
freed individually. It's also worth noting that `gk_free` cannot free any
allocation made outside of it's init/cleanup scope.
